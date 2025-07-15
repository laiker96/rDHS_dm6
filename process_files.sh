#! /bin/bash
WORKING_DIRECTORY=$1
METADATA_FILE=$2

mkdir -p "$WORKING_DIRECTORY"/macs_peaks

macs3_run() {
    local sample_id="$1"
    local workdir="$2"

    echo "Running MACS3 for ${sample_id} in ${workdir}"

    macs3 callpeak \
        -t "${workdir}/bams/${sample_id}_out.bam" \
        -f BAM \
        --nomodel \
        -n "${sample_id}" \
        --outdir "${workdir}/macs_peaks/${sample_id}" \
        --shift -75 \
        --extsize 150 \
        --keep-dup all \
        -B
}

macs3_coverage() {
    local sample_id="$1"
    local workdir="$2"

    macs3 bdgcmp -t "${workdir}/macs_peaks/${sample_id}/${sample_id}_treat_pileup.bdg" \
	    -c "${workdir}/macs_peaks/${sample_id}/${sample_id}_control_lambda.bdg" \
	    -m qpois -o "${workdir}/macs_peaks/${sample_id}/${sample_id}_qpois.bdg"
}

process_qpois_file() {
    local infile="$1"
    local outfile="$2"

    if [[ -z "$infile" || -z "$outfile" ]]; then
        echo "Usage: process_qpois_file <input_file> <output_file>"
        return 1
    fi

LC_NUMERIC=C awk 'BEGIN { OFS = "\t"; log10 = log(10) }
     {
         $4 = exp(-$4 * log10);
         print
     }' "$infile" > "$outfile"
}

intersect_intervals() {
    local sample_id="$1"
    local workdir="$2"

    echo "Intersecting intervals for ${sample_id}"

    bedtools intersect -u -wa -f 1.0 \
        -a "${workdir}/macs_peaks/${sample_id}/${sample_id}_qpois_proc.bdg" \
        -b "${workdir}/macs_peaks/${sample_id}/${sample_id}_peaks.narrowPeak" \
        > "${workdir}/macs_peaks/${sample_id}/${sample_id}_qpois.bed"
}

process_enrichment_file() {
    local enrich="$1"
    local workdir="$2"
    local dataDir="${workdir}/macs_peaks/${enrich}"
    local minP=325
    local scriptDir=~/data_dm6_PhD/rDHS_dm6/pipeline_rATAC_peaks/Toolkit

    echo "Step 1 ..."
    cp "$dataDir/${enrich}_qpois.bed" $dataDir/tmp.1
    mkdir -p ${dataDir}/tmp_files

    for j in $(seq 2 1 $minP); do
        cutoff=$(awk "BEGIN { printf \"1E-%d\", $j }")
        echo "$cutoff"
        
        python3 "$scriptDir/filter-long-double.py" "${dataDir}/tmp.1" "$cutoff" > ${dataDir}/tmp.2

        bedtools merge -d 1 -c 4 -o min -i ${dataDir}/tmp.2 | \
            awk '{if ($3-$2 >= 50) print $0}' > "${dataDir}/tmp_files/$enrich.$cutoff.bed"
        
        mv ${dataDir}/tmp.2 ${dataDir}/tmp.1

        num=$(wc -l < "${dataDir}/tmp_files/$enrich.$cutoff.bed")
        echo "$cutoff $num"
    done

    echo "Step 2 ..."
    cutoff="1E-2"
    awk '{if ($3 - $2 <= 400) print $0}' "${dataDir}/tmp_files/$enrich.$cutoff.bed" > ${dataDir}/peaks

    for j in $(seq 3 1 $minP); do
        cutoff=$(awk "BEGIN { printf \"1E-%d\", $j }")
        
        bedtools intersect -v -a "${dataDir}/tmp_files/$enrich.$cutoff.bed" -b ${dataDir}/peaks > ${dataDir}/tmp
        awk '{if ($3 - $2 <= 400) print $0}' ${dataDir}/tmp >> ${dataDir}/peaks
    done

    mv ${dataDir}/peaks "${dataDir}/${enrich}.DHSs.bed"

    bedtools intersect -v -a "${dataDir}/tmp_files/$enrich.1E-$minP.bed" -b "${dataDir}/${enrich}.DHSs.bed" > "${dataDir}/${enrich}.Excluded.bed"

    mkdir -p "$workdir/DHSs"
    mv "${dataDir}/${enrich}.Excluded.bed" "${dataDir}/${enrich}.DHSs.bed" "$workdir/DHSs/"

    rm -r ${dataDir}/tmp*
}

process_dhs_signal() {
    local id="$1"
    local workdir="$2"

    local dataDir="${workdir}/DHSs"
    local scriptDir="${workdir}/rDHS_dm6/pipeline_rATAC_peaks/Toolkit"
    local dhs="${dataDir}/${id}.DHSs.bed"
    local signal="${workdir}/signal_files/S3norm_rc_bedgraph/${id}_S3.bw"

    if [[ ! -f "$dhs" || ! -f "$signal" ]]; then
        echo "Error: Missing input file(s) for ID '$id'"
        echo "Expected DHS: $dhs"
        echo "Expected Signal: $signal"
        return 1
    fi

    echo "Processing ID: $id"

    local tmp_bed="${dataDir}/${id}_tmp.bed"
    local new="${dataDir}/${id}.new"
    local new_bed="${dataDir}/${id}.new.bed"
    local out_tab="${dataDir}/${id}.out.tab"
    local tmp1="${dataDir}/${id}.tmp.1"
    local output="${dataDir}/Processed-DHSs/output.${id}"

    mkdir -p "${dataDir}/Processed-DHSs/"

    cp "$dhs" "$tmp_bed"

    awk -v id="$id" '{print $1 "\t" $2 "\t" $3 "\t" id "-" NR "\t" $4}' "$tmp_bed" | sort -k4,4 > "$new"
    awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' "$new" > "$new_bed"

    bigWigAverageOverBed "$signal" "$new_bed" "$out_tab"
    sort -k1,1 "$out_tab" > "$tmp1"

    paste "$new" "$tmp1" | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $10}' | \
        awk '$5 != 1' > "$output"

    # Clean up
    rm -f "$tmp_bed" "$new" "$new_bed" "$out_tab" "$tmp1"
}


export -f macs3_run
export -f macs3_coverage
export -f process_qpois_file
export -f intersect_intervals
export -f process_enrichment_file
export -f process_dhs_signal

cut -f 5 -d , "$METADATA_FILE" | tail -n +2 \
    | parallel macs3_run {} "$WORKING_DIRECTORY"

cut -f 5 -d , "$METADATA_FILE" | tail -n +2 \
    | parallel macs3_coverage {} "$WORKING_DIRECTORY"

cut -f 5 -d , "$METADATA_FILE" | tail -n +2 \
    | parallel process_qpois_file "$WORKING_DIRECTORY/macs_peaks/{}/{}_qpois.bdg" "$WORKING_DIRECTORY/macs_peaks/{}/{}_qpois_proc.bdg"

cut -f 5 -d , "$METADATA_FILE" | tail -n +2 \
    | parallel intersect_intervals {} "$WORKING_DIRECTORY"

cut -f 5 -d , "$METADATA_FILE" | tail -n +2 \
    | parallel process_enrichment_file {} "$WORKING_DIRECTORY"

cut -f 5 -d , "$METADATA_FILE" | tail -n +2 \
    | parallel process_dhs_signal {} "$WORKING_DIRECTORY"
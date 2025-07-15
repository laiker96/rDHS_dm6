#! /bin/bash
WORKING_DIRECTORY=$1
METADATA_FILE=$2

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

export -f process_dhs_signal

cut -f 5 -d , "$METADATA_FILE" | tail -n +2 \
    | parallel process_dhs_signal {} "$WORKING_DIRECTORY"


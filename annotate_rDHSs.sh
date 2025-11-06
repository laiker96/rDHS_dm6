#!/bin/bash
set -euo pipefail

# ------------------ HELP & ARGUMENT PARSING ------------------ #
usage() {
    echo "Usage: $0 -w working_directory -m metadata_file -g chrom_sizes -a h3k27ac_dir -e metadata_h3k27ac"
    exit 1
}

while getopts ":w:m:g:a:e:" opt; do
  case $opt in
    w) working_directory="$OPTARG" ;;
    m) metadata_file="$OPTARG" ;;
    g) chrom_sizes="$OPTARG" ;;
    a) h3k27ac_dir="$OPTARG" ;;
    e) metadata_h3k27ac="$OPTARG" ;;
    *) usage ;;
  esac
done

if [[ -z "${working_directory:-}" || -z "${metadata_file:-}" || -z "${chrom_sizes:-}" || -z "${h3k27ac_dir:-}" || -z "${metadata_h3k27ac:-}" ]]; then
    usage
fi

# ------------------ FUNCTIONS ------------------ #

sort_dhs_outputs() {
    local workdir="$1"

    for bed in "${workdir}/DHSs/Processed-DHSs/output"*; do
        sort-bed --max-mem 8G "$bed" > tmp && mv tmp "$bed"
    done

    sort-bed --max-mem 8G "${workdir}/DHSs/rPeaks.bed" > "${workdir}/DHSs/tmp_sort"
    mv "${workdir}/DHSs/tmp_sort" "${workdir}/DHSs/rPeaks.bed"
}

build_reference_bed() {
    local workdir="$1"
    local metadata="$2"

    cut -f 1-4 "${workdir}/DHSs/rPeaks.bed" > "${workdir}/reference.bed"

    cat "$metadata" | while read -r sample_id; do
        bedmap --echo \
               --indicator \
               --fraction-either 0.6 \
               --delim '\t' \
               "${workdir}/reference.bed" \
               "${workdir}/DHSs/Processed-DHSs/output.${sample_id}" \
               > "${workdir}/tmp"
        mv "${workdir}/tmp" "${workdir}/reference.bed"
    done
}

generate_intervals() {
    local workdir="$1"
    local chrom_sizes="$2"

    local peaks="${workdir}/DHSs/rPeaks.bed"

    cut -f 1-4 "$peaks" | bedtools flank -i - -g "$chrom_sizes" -l 500 -r 0 > "${workdir}/left_intervals.bed"
    cut -f 1-4 "$peaks" | bedtools flank -i - -g "$chrom_sizes" -r 500 -l 0 > "${workdir}/right_intervals.bed"
    cut -f 1-4 "$peaks" > "${workdir}/central_intervals.bed"

    cp "${workdir}/central_intervals.bed" "${workdir}/dataset_right.bed"
    cp "${workdir}/central_intervals.bed" "${workdir}/dataset_left.bed"
    cp "${workdir}/central_intervals.bed" "${workdir}/dataset_central.bed"
}

process_signal() {
    local region="$1"
    local sample_id="$2"
    local workdir="$3"
    local bw_dir="$4"

    local interval_file="${workdir}/${region}_intervals.bed"
    local dataset_file="${workdir}/dataset_${region}.bed"
    local bw="${bw_dir}/${sample_id}.bw"
    local out_bed="${workdir}/out.bed"
    local out_tab="${workdir}/out.tab"

    bigWigAverageOverBed -bedOut="$out_bed" "$bw" "$interval_file" "$out_tab"
    cut -f 5 "$out_bed" | paste "$dataset_file" - > "${workdir}/tmp"
    mv "${workdir}/tmp" "$dataset_file"
    rm -f "$out_bed" "$out_tab"
}

process_all_regions() {
    local workdir="$1"
    local bw_dir="$2"
    local metadata="$3"

    cat "$metadata" | while read -r sample_id; do
        for region in left right central; do
            process_signal "$region" "$sample_id" "$workdir" "$bw_dir"
        done
    done
}

cleanup_intervals() {
    local workdir="$1"
    rm -f "${workdir}/left_intervals.bed" \
          "${workdir}/right_intervals.bed" \
          "${workdir}/central_intervals.bed"
}

# ------------------ MAIN EXECUTION ------------------ #

sort_dhs_outputs "$working_directory"
build_reference_bed "$working_directory" "$metadata_file"
generate_intervals "$working_directory" "$chrom_sizes"
process_all_regions "$working_directory" "$h3k27ac_dir" "$metadata_h3k27ac"
cleanup_intervals "$working_directory"

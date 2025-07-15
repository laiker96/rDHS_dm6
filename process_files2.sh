#! /bin/bash
WORKING_DIRECTORY=$1

#!/bin/bash

finalize_dhs() {
    local workdir="$1"

    local dhs_dir="${workdir}/DHSs"
    local processed_dir="${dhs_dir}/Processed-DHSs"
    local script_dir="${workdir}/rDHS_dm6/pipeline_rATAC_peaks/Toolkit"
    local all_bed="${dhs_dir}/DHS-All.bed"
    local filtered="${dhs_dir}/DHS-Filtered.bed"
    local tmp_sorted="${dhs_dir}/sorted"
    local rpeaks="${dhs_dir}/rPeaks"
    local cutoff=10

    echo "Combining DHSs from: $processed_dir"
    cat "${processed_dir}"/output.* > "$all_bed"

    echo "Filtering DHSs..."
    awk -v c="$cutoff" '{
        if ($1 !~ /_/ && $3 - $2 >= 150 && $6 > c && $5 <= 0.001)
            print $0
    }' "$all_bed" | grep -v -E 'chrM|chrY' > "$filtered"

    echo "Sorting DHSs..."
    sort -k1,1 -k2,2n "$filtered" > "$tmp_sorted"
    cp "$filtered" "${dhs_dir}/tmp.bed"
    > "$rpeaks"

    echo "Merging DHSs..."
    while [[ $(wc -l < "$tmp_sorted") -gt 0 ]]; do
        echo -e "\t$(wc -l < "$tmp_sorted") remaining..."

        bedtools merge -i "$tmp_sorted" -c 4,6 -o collapse,collapse > "${dhs_dir}/merge"
        python3 "${script_dir}/pick-best-peak.py" "${dhs_dir}/merge" > "${dhs_dir}/peak-list"

        awk 'FNR==NR {x[$1]; next} ($4 in x)' "${dhs_dir}/peak-list" "$tmp_sorted" >> "$rpeaks"
        bedtools intersect -v -a "$tmp_sorted" -b "$rpeaks" > "${dhs_dir}/remaining"

        mv "${dhs_dir}/remaining" "$tmp_sorted"
    done

    mv "$rpeaks" "${dhs_dir}/rPeaks.bed"

    echo "Cleaning up..."
    rm -f "${dhs_dir}/merge" "${dhs_dir}/peak-list" "$tmp_sorted" "${dhs_dir}/tmp.bed"
    echo "Final DHS set written to: ${dhs_dir}/rPeaks.bed"
}

# Export if needed for GNU parallel or external sourcing
export -f finalize_dhs
finalize_dhs $WORKING_DIRECTORY

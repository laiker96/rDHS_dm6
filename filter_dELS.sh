#!/bin/bash
set -euo pipefail

# ------------------ USAGE FUNCTION ------------------ #

usage() {
    echo "Usage: $0 -c CREs.bed -t TSS.bed -g genome_file -m metadata.csv"
    exit 1
}

# ------------------ ARGUMENT PARSING WITH getopts ------------------ #

while getopts ":c:t:g:m:" opt; do
  case ${opt} in
    c) CRES="$OPTARG" ;;
    t) TSS_BED="$OPTARG" ;;
    g) GENOME_FILE="$OPTARG" ;;
    m) METADATA_FILE="$OPTARG" ;;
    *) usage ;;
  esac
done

# Check if required arguments are set
if [[ -z "${CRES:-}" || -z "${TSS_BED:-}" || -z "${GENOME_FILE:-}" || -z "${METADATA_FILE:-}" ]]; then
    echo "Error: Missing required argument(s)."
    usage
fi

# ------------------ FUNCTIONS ------------------ #

define_output_file() {
    local cres="$1"
    local cre_dir
    cre_dir=$(dirname "$cres")
    echo "${cre_dir}/dELS.bed"
}

generate_dels() {
    local cres="$1"
    local tss="$2"
    local genome="$3"
    local output="$4"

    bedtools slop -b 500 -i "$tss" -g "$genome" \
        | bedtools intersect -wa -v -a "$cres" -b - > "$output"
}

split_dels_by_sample() {
    local metadata="$1"
    local dels_file="$2"
    local output_dir
    output_dir=$(dirname "$dels_file")

    local column_index=5
    cat "$metadata" | while read -r sample_id; do
        awk -v col="$column_index" '$col == 1' "$dels_file" | cut -f1-4 > "${output_dir}/dELS_${sample_id}.bed"
        ((column_index++))
    done
}

# ------------------ MAIN ------------------ #

main() {
    local cres="$1"
    local tss="$2"
    local genome="$3"
    local metadata="$4"

    local output_file
    output_file=$(define_output_file "$cres")

    generate_dels "$cres" "$tss" "$genome" "$output_file"
    split_dels_by_sample "$metadata" "$output_file"
}

main "$CRES" "$TSS_BED" "$GENOME_FILE" "$METADATA_FILE"

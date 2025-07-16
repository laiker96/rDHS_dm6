#!/bin/bash

CREs=$1
TSS_BED=$2
GENOME_FILE=$3
METADATA_FILE=$4

CREs_DIR=$(dirname "${CREs}")

OUTFILE="${CREs_DIR}/dELS.bed"

bedtools slop -b 500 -i "${TSS_BED}" -g "${GENOME_FILE}" \
    | bedtools intersect -wa -v -a "${CREs}" -b - > "${OUTFILE}"

COLUMN_INDEX=5  
tail -n +2 "${METADATA_FILE}" | cut -f 1 -d , | while read LINE; do
    awk -v col="$COLUMN_INDEX" '$col == 1' "${CREs}" | cut -f1-4 > "${CREs_DIR}/dELS_${LINE}.bed"
    ((COLUMN_INDEX++))
done
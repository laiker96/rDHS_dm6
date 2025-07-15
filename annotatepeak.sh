#! /bin/bash

WORKING_DIRECTORY=$1
METADATA_FILE=$2
GENOME_DIR=$3
H3K27ac_DIR=$4

for BED in ${WORKING_DIRECTORY}/DHSs/Processed-DHSs/output*; do            
sort-bed --max-mem 32G $BED > tmp && mv tmp $BED
done

cut -f 1-4 ${WORKING_DIRECTORY}/DHSs/rPeaks.bed > ${WORKING_DIRECTORY}/reference.bed

while read LINE; do  
    bedmap --echo \
           --indicator \
           --fraction-both 0.5 \
           --delim '\t' \
           "${WORKING_DIRECTORY}/reference.bed" \
           "${WORKING_DIRECTORY}/DHSs/Processed-DHSs/output.${LINE}" \
           > "${WORKING_DIRECTORY}/tmp" && \
    mv "${WORKING_DIRECTORY}/tmp" "${WORKING_DIRECTORY}/reference.bed"
done < <(cut -f 5 -d , "${METADATA_FILE}" | tail -n +2)


cut -f 1-4 ${WORKING_DIRECTORY}/DHSs/rPeaks.bed | bedtools flank -i - -g ${GENOME_DIR}/dm6.chrom.sizes -l 500 -r 0 > ${WORKING_DIRECTORY}/left_intervals.bed
cut -f 1-4 ${WORKING_DIRECTORY}/DHSs/rPeaks.bed | bedtools flank -i - -g ${GENOME_DIR}/dm6.chrom.sizes -r 500 -l 0 > ${WORKING_DIRECTORY}/right_intervals.bed

cut -f 1-4 ${WORKING_DIRECTORY}/DHSs/rPeaks.bed > ${WORKING_DIRECTORY}/dataset_right.bed
cut -f 1-4 ${WORKING_DIRECTORY}/DHSs/rPeaks.bed > ${WORKING_DIRECTORY}/dataset_left.bed
cut -f 1-4 ${WORKING_DIRECTORY}/DHSs/rPeaks.bed > ${WORKING_DIRECTORY}/dataset_central.bed

cp ${WORKING_DIRECTORY}/dataset_central.bed ${WORKING_DIRECTORY}/central_intervals.bed

for FILE in ${H3K27ac_DIR}/*.bw; do
    bigWigAverageOverBed -bedOut=${WORKING_DIRECTORY}/out.bed "$FILE" "${WORKING_DIRECTORY}/right_intervals.bed" ${WORKING_DIRECTORY}/out.tab && \
    rm ${WORKING_DIRECTORY}/out.tab && \
    cut -f 5 ${WORKING_DIRECTORY}/out.bed | paste ${WORKING_DIRECTORY}/dataset_right.bed - > ${WORKING_DIRECTORY}/tmp && \
    mv ${WORKING_DIRECTORY}/tmp ${WORKING_DIRECTORY}/dataset_right.bed && \
    rm ${WORKING_DIRECTORY}/out.bed
done

for FILE in ${H3K27ac_DIR}/*.bw; do
    bigWigAverageOverBed -bedOut=${WORKING_DIRECTORY}/out.bed "$FILE" "${WORKING_DIRECTORY}/left_intervals.bed" ${WORKING_DIRECTORY}/out.tab && \
    rm ${WORKING_DIRECTORY}/out.tab && \
    cut -f 5 ${WORKING_DIRECTORY}/out.bed | paste ${WORKING_DIRECTORY}/dataset_left.bed - > ${WORKING_DIRECTORY}/tmp && \
    mv ${WORKING_DIRECTORY}/tmp ${WORKING_DIRECTORY}/dataset_left.bed && \
    rm ${WORKING_DIRECTORY}/out.bed
done

for FILE in ${H3K27ac_DIR}/*.bw; do
    bigWigAverageOverBed -bedOut=${WORKING_DIRECTORY}/out.bed "$FILE" "${WORKING_DIRECTORY}/central_intervals.bed" ${WORKING_DIRECTORY}/out.tab && \
    rm ${WORKING_DIRECTORY}/out.tab && \
    cut -f 5 ${WORKING_DIRECTORY}/out.bed | paste ${WORKING_DIRECTORY}/dataset_central.bed - > ${WORKING_DIRECTORY}/tmp && \
    mv ${WORKING_DIRECTORY}/tmp ${WORKING_DIRECTORY}/dataset_central.bed && \
    rm ${WORKING_DIRECTORY}/out.bed
done

rm ${WORKING_DIRECTORY}/left_intervals.bed ${WORKING_DIRECTORY}/right_intervals.bed ${WORKING_DIRECTORY}/central_intervals.bed

#!/bin/bash

dir=~/data_dm6_PhD/DHSs
scriptDir=~/data_dm6_PhD/pipeline_rATAC_peaks/Toolkit/

cd $dir/Processed-DHSs

echo -e "Combining DHSs..."
cat output.* > tmp
mv tmp $dir/DHS-All.bed

cutoff=10 #remove zero signal peaks

echo -e "Filtering DHSs..."
cd $dir
awk '{if ($1 !~ /_/ && $3-$2 >= 150 && $6 > '$cutoff' && $5 <= 0.001) print $0}' \
    DHS-All.bed | grep -v "chrM" | grep -v "chrY" > DHS-Filtered.bed

cp DHS-Filtered.bed tmp.bed

echo -e "Sorting DHSs..."
sort -k1,1 -k2,2n tmp.bed > sorted
num=$(wc -l sorted | awk '{print $1}')

echo -e "Merging DHSs..."
while [ $num -gt 0 ]
do
    echo -e "\t" $num
    bedtools merge -i sorted -c 4,6 -o collapse,collapse > merge
    python3 $scriptDir/pick-best-peak.py merge > peak-list
    awk 'FNR==NR {x[$1];next} ($4 in x)' peak-list sorted >> rPeaks
    bedtools intersect -v -a sorted -b rPeaks > remaining
    mv remaining sorted
    num=$(wc -l sorted | awk '{print $1}')
done

rm merge peak-list sorted tmp.bed 
mv rPeaks rPeaks.bed

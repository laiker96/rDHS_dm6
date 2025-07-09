#!/bin/bash

enrich=$1
dataDir=~/data_dm6_PhD
minP=325
scriptDir=~/data_dm6_PhD/pipeline_rATAC_peaks/Toolkit
    
echo "Step 1 ..."
cp $enrich.bed tmp.1
mkdir -p tmp_files

for j in `seq 2 1 $minP`
do
	cutoff=$(awk 'BEGIN{print "1E-'$j'"}')
	echo $cutoff 
	python3 $scriptDir/filter-long-double.py tmp.1 $cutoff > tmp.2
	bedtools merge -d 1 -c 4 -o min -i tmp.2 | \
	    awk '{if ($3-$2 >= 50) print $0}' > tmp_files/$enrich.$cutoff.bed
	mv tmp.2 tmp.1
	num=$(wc -l tmp_files/$enrich.$cutoff.bed | awk '{print $1}')
	echo $cutoff $num
done

echo "Step 2 ..." 

cutoff=1E-2
awk '{if ($3-$2 <= 400) print $0}' tmp_files/$enrich.$cutoff.bed > peaks

for j in `seq 3 1 $minP`
	do
	cutoff=$(awk 'BEGIN{print "1E-'$j'"}')
	bedtools intersect -v -a tmp_files/$enrich.$cutoff.bed -b peaks > tmp
	awk '{if ($3-$2 <= 400) print $0}' tmp >> peaks
done

mv peaks $enrich.DHSs.bed

bedtools intersect -v -a tmp_files/$enrich.1E-$minP.bed -b $enrich.DHSs.bed > $enrich.Excluded.bed

mkdir -p $dataDir/DHSs

mv $enrich.Excluded.bed $enrich.DHSs.bed $dataDir/DHSs/
rm -r tmp*

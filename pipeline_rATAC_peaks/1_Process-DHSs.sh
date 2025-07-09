#!/bin/bash

dataDir=~/data_dm6_PhD/DHSs
scriptDir=~/data_dm6_PhD/pipeline_rATAC_peaks/Toolkit/


dhs=$1
signal=$2
id=$3

cp $dhs bed
awk '{print $1 "\t" $2 "\t" $3 "\t" "'$id'-"NR "\t" $4}' bed | sort -k4,4 > new
awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' new > new.bed
bigWigAverageOverBed $signal new.bed out.tab
sort -k1,1 out.tab > tmp.1
paste new tmp.1 | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $10}' | awk '$5 != 1' >  output.$id

mkdir -p $dataDir/Processed-DHSs/
mv output.$id $dataDir/Processed-DHSs/
rm new new.bed out.tab tmp.1 bed

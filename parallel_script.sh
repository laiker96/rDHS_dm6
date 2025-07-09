for BAM in *.bam; do
alignmentSieve --bam $BAM --outFile $(basename $BAM .bam)_size_filter.bam -p 24 --maxFragmentLength 150 --minFragmentLength 40 --ATACshift
samtools sort --threads 24 $(basename $BAM .bam)_size_filter.bam > tmp && mv tmp $(basename $BAM .bam)_size_filter.bam
samtools index -@ 24 $(basename $BAM .bam)_size_filter.bam
done


ls bams/*.bam | parallel 'macs3 callpeak -t {} -f BAM --nomodel -n {/.} --outdir macs_peaks/{/.} --shift -75 --extsize 150 --keep-dup all -B'

parallel 'macs3 bdgcmp -t peak_files/peaks_all_contexts_size_filter/{}/{}_treat_pileup.bdg -c peak_files/peaks_all_contexts_size_filter/{}/{}_control_lambda.bdg -m qpois -o peak_files/peaks_all_contexts_size_filter/{}/{}_qpois.bdg' < <(cut -f 5 -d , metadata/contexts_complete_data.csv | tail -n +2)

parallel 'Rscript pipeline_rATAC_peaks/process_qpois.R --file peak_files/peaks_all_contexts_size_filter/{}/{}_qpois.bdg --outname peak_files/peaks_all_contexts_size_filter/{}/{}_qpois_proc.bdg' < <(cut -f 5 -d , metadata/contexts_complete_data.csv | tail -n +2)

parallel 'bedtools intersect -u -wa -f 1.0 -a peak_files/peaks_all_contexts_size_filter/{}/{}_qpois_proc.bdg -b peak_files/peaks_all_contexts/{}/{}_peaks.narrowPeak > peak_files/peaks_all_contexts_size_filter/{}/{}_qpois.bed' < <(cut -f 5 -d , metadata/contexts_complete_data.csv | tail -n +2)

parallel 'cd peak_files/peaks_all_contexts_size_filter/{} && bash $HOME/data_dm6_PhD/pipeline_rATAC_peaks/0_Call-DHSs.sh {}_qpois && cd ../../..' < <(cut -f 5 -d , metadata/contexts_complete_data.csv | tail -n +2)

parallel 'cd peak_files/peaks_all_contexts_size_filter/{} && bash $HOME/data_dm6_PhD/pipeline_rATAC_peaks/1_Process-DHSs.sh $HOME/data_dm6_PhD/DHSs/{}_qpois.DHSs.bed $HOME/data_dm6_PhD/signal_files/size_filtered/standard_all_contexts/S3norm_rc_bedgraph/{}_S3.bw {} && cd ../../..' < <(cut -f 5 -d , metadata/contexts_complete_data.csv | tail -n +2)

bash pipeline_rATAC_peaks/2_Create-rDHSs.sh


for BED in DHSs/*.bed; do            
sort-bed --max-mem 32G $BED > tmp && mv tmp $BED
done

cut -f 1-4 DHSs/rPeaks.bed > reference.bed

while read LINE; do  
bedmap --echo --indicator --fraction-both 0.5 --delim '\t' reference.bed DHSs/"$LINE"_qpois.DHS.bed > tmp && mv tmp reference.bed 
done < <(cut -f 5 -d , contexts_complete_data.csv | tail -n +2)

cut -f 1-4 DHSs/rPeaks.bed | bedtools flank -i - -g genome_files/dm6.chrom.sizes -l 500 -r 0 > left_intervals.bed
cut -f 1-4 DHSs/rPeaks.bed | bedtools flank -i - -g genome_files/dm6.chrom.sizes -r 500 -l 0 > right_intervals.bed

cut -f 1-4 DHSs/rPeaks.bed > dataset_right.bed
cut -f 1-4 DHSs/rPeaks.bed > dataset_left.bed
cut -f 1-4 DHSs/rPeaks.bed > dataset_central.bed
cp dataset_central.bed central_intervals.bed

for FILE in signal_files/H3K27ac/*.bw; do           
bigWigAverageOverBed -bedOut=out.bed $FILE right_intervals.bed out.tab && rm out.tab && cut -f 5 out.bed | paste dataset_right.bed - > tmp && mv tmp dataset_right.bed && rm out.bed
done

for FILE in signal_files/H3K27ac/*.bw; do                                                             
bigWigAverageOverBed -bedOut=out.bed $FILE left_intervals.bed out.tab && rm out.tab && cut -f 5 out.bed | paste dataset_left.bed - > tmp && mv tmp dataset_left.bed && rm out.bed   
done

for FILE in signal_files/H3K27ac/*.bw; do                                                             
bigWigAverageOverBed -bedOut=out.bed $FILE central_intervals.bed out.tab && rm out.tab && cut -f 5 out.bed | paste dataset_central.bed - > tmp && mv tmp dataset_central.bed && rm out.bed   
done

rm left_intervals.bed right_intervals.bed central_intervals.bed


bedtools slop -b 500 -i genome_files/dm6_TSS_RefSeqNCBI.bed -g genome_files/dm6.chrom.sizes | bedtools intersect -wa -v -a cCREs.bed -b - > dELS.bed

awk '$5 == 1' dELS.bed | cut -f 1-4 > dELS_AB.bed
awk '$6 == 1' dELS.bed | cut -f 1-4 > dELS_E11.bed 
awk '$7 == 1' dELS.bed| cut -f 1-4 > dELS_E13.bed 
awk '$8 == 1' dELS.bed | cut -f 1-4 > dELS_E5.bed
awk '$9 == 1' dELS.bed | cut -f 1-4 > dELS_EAD.bed 
awk '$10 == 1' dELS.bed | cut -f 1-4 > dELS_HID.bed
awk '$11 == 1' dELS.bed | cut -f 1-4 > dELS_LB.bed
awk '$12 == 1' dELS.bed | cut -f 1-4 > dELS_O.bed
awk '$13 == 1' dELS.bed | cut -f 1-4 > dELS_WID.bed

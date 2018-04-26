#!/bin/bash

# sorting raw (no trimming) BAMs
module load samtools

alignment_path='/home/CBSB_UofK.RNA_seq_assesment/Alignment/align_2pass_raw'
output_path='/home/CBSB_UofK.RNA_seq_assesment/DE/sorted-alignments/raw'

for i in `seq 1 6`;
do
	samtools sort -n -o $output_path/sample_${i}Aligned.sorted.out.bam $alignment_path/sample_${i}Aligned.out.bam 
done  

#!/bin/bash

module load stringtie

alignment_path='/home/CBSB_UofK.RNA_seq_assesment/Alignment/align_2pass_trimmed'
ref_annotation='/home/CBSB_UofK.RNA_seq_assesment/Genome/gencode.v27.chr_patch_hapl_scaff.annotation.gtf'
output_dir='/home/CBSB_UofK.RNA_seq_assesment/DE/stringtie-output/partial'

for i in `seq 1 6`;
do
	stringtie -p 8 -B -G $ref_annotation -o $output_dir/sample${i}/sample${i}.gtf $alignment_path/sample_${i}Aligned.sortedByCoord.out.bam
done

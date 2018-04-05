#!/bin/bash

tm_path='/home/users_current/assesuser7/tablemaker-2.1.1.Linux_x86_64'
alignment_path='/home/CBSB_UofK.RNA_seq_assesment/Alignment/align.clean.1pass'
ref_annotation='/home/CBSB_UofK.RNA_seq_assesment/Genome/gencode.v27.chr_patch_hapl_scaff.annotation.gtf'
output_dir='/home/CBSB_UofK.RNA_seq_assesment/DE/tablemaker-output/1pass/partial-trim'

for i in `seq 1 6`;
do
	$tm_path/tablemaker -p 4 -q -W -G $ref_annotation -o $output_dir/sample${i}_output $alignment_path/sample_${i}Aligned.sortedByCoord.out.bam
done

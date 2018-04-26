#!/bin/bash

set -x
# Defining path
project_dir=/home/aeahmed/rna_assement/
transcriptome_dir=${project_dir}/refernce/transcriptome
transcriptome_file=Homo_sapiens.GRCh38.cdna.all.fa.gz
genome_dir=${project_dir}/genome
reads=${project_dir}/reads/
quant_result=${project_dir}/quant

set +x
# Loading needed modules
module load kallisto/0.44.0

set -x

echo "Indexing start: `date`"
# Building and index:
kallisto index -i ${transcriptome_dir}/transcripts.idx ${transcriptome_dir}/${transcriptome_file}

# Quantification (per sample):
for i in `seq 1 6`;
do
        echo "**************************** Running on sample ${i} ******************"
	echo "Sample ${i} start: `date`"
	mkdir -p ${quant_result}/sample_${i}
	kallisto quant -i ${transcriptome_dir}/transcripts.idx \
		 -o ${quant_result}/sample_${i} \
		 -b 100  \
		  <(zcat $reads/sample${i}_R1.fastq.gz) <(zcat $reads/sample${i}_R2.fastq.gz)


done

echo "Analysis seems done!" | mail -s "RNA seq exercise" azzaea@gmail.com

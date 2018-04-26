#!/bin/bash


module load  nextflow/0.27.6.4775
module load R/3.2.3
module load fastqc/0.11.5 
module load TrimGalore/1.3.4c 
module load STAR/2.5.3a 
module load subread/1.6.0 
module load picard-tools/2.6.0
module load samtools/1.3.1 
module load stringtie/1.3.4c 
module load RSeQC/2.6.4 
module load preseq/2.0.0 
module load bedops/2.4.20
module load hisat2/2.1.0
module load multiqc/1.5

rm  $PWD/results/NGI-RNAseq_t_*

NXF_WORK="/home/aeahmed/tmp_nextflow"
export PICARD_HOME="/usr/src/picard-tools/picard-tools-2.6.0/"

nextflow run main.nf \
	--reads '/home/aeahmed/rna_assement/reads/*_R{1,2}.fastq.gz'\
	--genome 'GRCh37' \
	--project 'test' \
	-profile standard\
	-with-dag execution_dag.html -resume

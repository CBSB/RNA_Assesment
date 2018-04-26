#!/bin/bash


module load  nextflow/0.27.6.4775
module load kallisto/0.44.0
module load R/3.4.3

rm -rf execution_*

NXF_WORK="/home/aeahmed/tmp_nextflow"

nextflow run main.nf -with-dag execution_dag.html -resume

echo "done slueth pipeline" | mail -s "slueth DGE" "azzaea@gmail.com"

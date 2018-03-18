#!/bin/bash


module load  nextflow/0.27.6.4775
module load kallisto/0.44.0

nextflow run kallisto.nf --resume -with-dag execution_dag.html

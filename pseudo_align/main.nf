/*
 * Copyright (c) 2015-2018, Centre for Genomic Regulation (CRG) and the authors.
 * Main Kallisto-NF pipeline script
 * @authors
 * Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * Evan Floden <evanfloden@gmail.com> 
 */

params.bootstrap = 100

log.info "K A L L I S T O - N F  ~  version 0.9"
log.info "====================================="
log.info "reads                  : ${params.reads}"
log.info "transcriptome          : ${params.reftranscriptome}"
log.info "experimental design    : ${params.metadata}"
log.info "\n"


/*
 * Input parameters validation
 */
transcriptome_file     = file(params.reftranscriptome)
exp_file               = file(params.metadata) 


/*
 * validate input files
 */
if( !transcriptome_file.exists() ) exit 1, "Missing transcriptome file: ${transcriptome_file}"
if( !exp_file.exists() ) exit 1, "Missing experimental design file: ${exp_file}"


/*
 * Create a channel for read files 
 */
Channel
    .fromFilePairs( params.reads, size: -1 )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_files } 


process index {
    publishDir "${params.pseudo}", mode: 'copy'
    
    input:
    file transcriptome_file
    
    output:
    file "transcriptome.index" into transcriptome_index
      
    script:

    """
    kallisto index -i transcriptome.index ${transcriptome_file}
    """
}


process kallisto {
    tag "reads: $name"

    publishDir "${params.pseudo}", mode: 'copy'

    input:
    file index from transcriptome_index
    set val(name), file(reads) from read_files

    output:
    file "kallisto_${name}" into kallisto_out_dirs 

    script:
     
    """
    mkdir kallisto_${name}
    kallisto quant -b ${params.bootstrap} -i ${index} -t 10 -o kallisto_${name} <(zcat ${reads[0]}) <(zcat ${reads[1]})
    """
}


process sleuth_vs_limma {
    publishDir "${params.diff}", mode: 'copy'

    input:
    file 'kallisto/*' from kallisto_out_dirs.collect()   
    file exp_file
    file transcriptome_file

    output: 
    
    file "*.{csv,pdf,so,Rda}" into sleuth_limma__results

    script:
 
    """
    sleuth_vs_limma.R kallisto ${exp_file} ${transcriptome_file}
    """
}


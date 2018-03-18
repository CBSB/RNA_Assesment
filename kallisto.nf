/*
 * Copyright (c) 2015-2018, Centre for Genomic Regulation (CRG) and the authors.
 * Main Kallisto-NF pipeline script
 * @authors
 * Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * Evan Floden <evanfloden@gmail.com> 
 */


//params.transcriptome = "$baseDir/test_data/transcriptome/transcriptome.fa"
params.transcriptome = "/home/aeahmed/rna_assement//refernce/transcriptome/Homo_sapiens.GRCh38.cdna.all.fa.gz"
params.name          = "RNA-Seq Abundance Analysis"
params.reads         = "/home/CBSB_UofK.RNA_seq_assesment/dataset1/*.fastq.gz"
params.fragment_len  = '180'
params.fragment_sd   = '20'
params.bootstrap     = '100'
params.experiment    = "/home/CBSB_UofK.RNA_seq_assesment/dataset1/README.txt"
params.output        = "results/"


log.info "K A L L I S T O - N F  ~  version 0.9"
log.info "====================================="
log.info "name                   : ${params.name}"
log.info "reads                  : ${params.reads}"
log.info "transcriptome          : ${params.transcriptome}"
log.info "fragment length        : ${params.fragment_len} nt"
log.info "fragment SD            : ${params.fragment_sd} nt"
log.info "bootstraps             : ${params.bootstrap}"
log.info "experimental design    : ${params.experiment}"
log.info "output                 : ${params.output}"
log.info "\n"


/*
 * Input parameters validation
 */
transcriptome_file     = file(params.transcriptome)
exp_file               = file(params.experiment) 


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
    input:
    file transcriptome_file
    
    output:
    file "transcriptome.index" into transcriptome_index
      
    script:
    //
    // Kallisto tools mapper index
    //
    module load 'kallisto/0.44.0'

    """
    kallisto index -i transcriptome.index ${transcriptome_file}
    """
}


process mapping {
    tag "reads: $name"

    input:
    file index from transcriptome_index
    set val(name), file(reads) from read_files

    output:
    file "kallisto_${name}" into kallisto_out_dirs 

    script:
    //
    // Kallisto tools mapper
    //

    module load 'kallisto/0.44.0'

    def single = reads instanceof Path
    if( !single ) {
        """
        mkdir kallisto_${name}
        kallisto quant -b ${params.bootstrap} -i ${index} -t ${task.cpus} -o kallisto_${name} ${reads}
        """
    }  
    else {
        """
        mkdir kallisto_${name}
        kallisto quant --single -l ${params.fragment_len} -s ${params.fragment_sd} -b ${params.bootstrap} -i ${index} -t ${task.cpus} -o kallisto_${name} ${reads}
        """
    }

}


process sleuth {
    input:
    file 'kallisto/*' from kallisto_out_dirs.collect()   
    file exp_file

    output: 
    file 'sleuth_object.so'
    file 'gene_table_results.txt'

    script:
    //
    // Setup sleuth R dependancies and environment
    //
 
    """
    sleuth.R kallisto ${exp_file}
    """
}




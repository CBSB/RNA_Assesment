
execute_dge_analysis <- function(bg) {
  ## DGE analysis
  # Stat test
  genes_test_results = stattest(bg, feature='gene', meas='FPKM', covariate='group')
  
  # Filter for results with qval less than 0.05
  diff_genes_results = genes_test_results[!is.nan(genes_test_results$qval) & genes_test_results$qval < 0.05,]
  
  # Genes IDs and Genes Names
  diff_genes_results_w_names = diff_genes_results
  diff_genes_results_w_names$geneName = geneNames(bg)[rownames(diff_genes_results)]
  
  return(diff_genes_results_w_names)
}


execute_dte_analysis <- function(bg) {
  ## DTE analysis
  # Stat test
  # Without libadjust=FALSE, no result gets past the 0.05 threshold
  trans_test_results = stattest(bg, feature='transcript', meas='FPKM', covariate='group', libadjust=FALSE)
  
  # Filter for results with qval less than 0.05
  diff_trans_results = trans_test_results[!is.nan(trans_test_results$qval) & trans_test_results$qval < 0.05,]
  
  # Join differentiated transcripts and gene IDs and Names
  
  transcript_gene_table = indexes(bg)$t2g
  
  # Dataframe containing transcripts IDs, genes IDs, and genes Names
  trans_gens_names = transcript_gene_table
  trans_gens_names$transcriptName = transcriptNames(bg)
  trans_gens_names$geneName = geneNames(bg)
  
  diff_trans_results_w_genes_names = merge(x = diff_trans_results, y = trans_gens_names, by.x = "id", by.y = "t_id")
  
  return(diff_trans_results_w_genes_names)
}


library(ballgown)

## Load data
## ===========
raw_bg = readRDS('raw_bg.rd')
partial_bg = readRDS('partial_bg.rd')
full_bg = readRDS('full_bg.rd')


## DTE analysis
## ===========
raw_dte_results = execute_dte_analysis(raw_bg)
partial_dte_results = execute_dte_analysis(partial_bg)
full_dte_results = execute_dte_analysis(full_bg)

View(head(raw_dte_results))
View(head(partial_dte_results))
View(head(full_dte_results))
#nrow(diff_trans_results_w_genes_names)


## DGE analysis
## ===========
raw_dge_results = execute_dge_analysis(raw_bg)
partial_dge_results = execute_dge_analysis(partial_bg)
full_dge_results = execute_dge_analysis(full_bg)

View(head(raw_dge_results))
View(head(partial_dge_results))
View(head(full_dge_results))
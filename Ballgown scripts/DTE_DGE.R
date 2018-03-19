
load('bg.rda')

library(ballgown)
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

View(head(diff_trans_results_w_genes_names))
nrow(diff_trans_results_w_genes_names)

## DGE analysis

# Stat test
genes_test_results = stattest(bg, feature='gene', meas='FPKM', covariate='group')

# Filter for results with qval less than 0.05
diff_genes_results = genes_test_results[!is.nan(genes_test_results$qval) & genes_test_results$qval < 0.05,]

# Genes IDs and Genes Names
diff_genes_results_w_names = diff_genes_results
diff_genes_results_w_names$geneName = geneNames(bg)[rownames(diff_genes_results)]

View(diff_genes_results_w_names)




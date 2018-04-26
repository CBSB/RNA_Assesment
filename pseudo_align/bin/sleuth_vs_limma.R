#!/usr/bin/env Rscript

# Load packages ----
library(dplyr)
library(readr) #general use package from Hadley Wickham for reading in tables
library(Biostrings) #we'll use this to edit our reference fasta file
library(tximport) #package for getting Kallisto results into R
library(ensembldb) ##used together with your organism-specific database package to get annotation info
library(EnsDb.Hsapiens.v86) #replace with your organism-specific database package
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(trelliscopejs)
library(genefilter)
library(limma) #powerful package for differential gene expression using linear modeling
library(edgeR)
library(Biobase)
library(sleuth)
library(tibble)
library(DT)
library(gplots) 
library(VennDiagram)
library(UpSetR)

# library(GSEABase)
# library(GSVA)
# library(Biobase)

# Read in the inputs:
args1="/home/aeahmed/rna_assement/results/pseudo_align/raw"
args2="/home/aeahmed/rna_assement/reads/README.txt"
args3="/home/aeahmed/rna_assement/refernce/transcriptome/Homo_sapiens.GRCh38.cdna.all.fa"

print("/n/n")
args <- commandArgs(TRUE)
print("First argument is:")
kallisto_output <- args[1]
kallisto_output
print("Second argument is:")
metadata_file <- args[2]
metadata_file
print("Third argument is:")
ref_transcriptome <- args[3]
ref_transcriptome

# Get annotations ----

#use the 'transcripts' function from the EnsemblDB package to get annotation info
Tx <- transcripts(EnsDb.Hsapiens.v86, 
                  columns=c(listColumns(EnsDb.Hsapiens.v86,
                                        "tx"), "gene_name"))

Tx <- as.data.frame(Tx)

#need to change first column name to 'target_id'
Tx <- dplyr::rename(Tx, target_id = tx_id)
row.names(Tx) <- NULL
head(Tx)

#transcript ID needs to be the first column in the dataframe
Tx <- Tx[,c(6,12)]

# Import Kallisto transcript counts into R using Tximport ----
sample_id <- dir(kallisto_output, pattern = "sample") 
print("samples:")
sample_id

kal_dirs <-  sapply(sample_id, function(id) file.path(kallisto_output, id))
print("kal_dirs is")
kal_dirs

targets <- read.table(metadata_file, row.names=NULL, header = T, as.is = T) %>% 
  dplyr::rename(sample = SampleID, condition = Condition) %>%
  mutate(path = kal_dirs)

files <- file.path(targets$path, "abundance.h5")

Txi_gene <- tximport(files, type = "kallisto", tx2gene = Tx,  ignoreTxVersion = TRUE,
                     txOut = TRUE, #Do DTE, to change to DGE, make it FALSE
                     countsFromAbundance = "lengthScaledTPM") 

#take a look at the object you just created
head(Txi_gene$counts)
head(Txi_gene$abundance)

# =============================================================================

# Normalizing the data (for limma): ----
targets
groups <- targets$condition
groups <- factor(groups)
sampleLabels <- targets$sample

# Examine your data up to this point ----
myTPM <- Txi_gene$abundance
myCounts <- Txi_gene$counts
# Do you think counts or TPM are appropriate for multivariate statistical analysis....why?

# Take a look at the heteroskedasticity of the data ----
# first, calculate row means and standard deviations for each transcript or gene and add these to your data frame
myTPM.stats <- transform(myTPM, SD=rowSds(myTPM), AVG=rowMeans(myTPM), MED=rowMedians(myTPM))
myCounts.stats <- transform(myCounts, SD=rowSds(myCounts), AVG=rowMeans(myCounts), MED=rowMedians(myCounts))

print("TPM stats:")
head(myTPM.stats)
tpm_plot <- ggplot(myTPM.stats, aes(x=SD, y=MED)) + geom_point()
ggsave(plot = tpm_plot, filename ="0.heteroskedasticity_tpmplot.pdf")

print("Count stats")
head(myCounts.stats)
count_plot <- ggplot(myCounts.stats, aes(x=SD, y=MED)) + geom_point()
ggsave(count_plot, filename = "0.heteroskedasticity_countplot.pdf")
# how might you expect that counts and TPM compare if used as input for PCA analysis?
# how would these graphs change if you log2 converted the data?

# Make a DGElist from your counts ----
DGEList <- DGEList(Txi_gene$counts)
# use the 'cpm' function from EdgeR to get counts per million
cpm <- cpm(DGEList) 
log2.cpm <- cpm(DGEList, log=TRUE)

# Take a look at the distribution of the Log2 CPM

nsamples <- ncol(log2.cpm)
myColors <- brewer.pal(nsamples, "Paired")
pdf("1.explorations.pdf")
boxplot(log2.cpm, las=2, cex=1, col = myColors, names = sampleLabels, main="non-normalized log2 cpm")
# Filter your data ----
#first, take a look at how many genes or transcripts have no read counts at all
print("Distribution of 0 count transcripts across samples")
table(rowSums(DGEList$counts==0))
data.frame(sum = rowSums(cpm>0)) %>%  group_by(sum) %>%  summarise(n = n()) %>%
  ggplot() + geom_col(aes(x = sum, y = n), fill = brewer.pal(nsamples+1, "Paired")) + 
  geom_text(aes(x = sum, y = n + 0.05*n , label = n)) + 
  xlab("Transcripts not expressed across samples") + ylab("Count") +
  ggtitle("Distribution of 0 count transcripts across samples")

print("Count of Transcripts with no counts (cpm also = 0) across all samples")
table(rowSums(DGEList$counts==0)==nsamples)
table(rowSums(cpm == 0)==nsamples)

print("Therefore, fraction of transcripts with no counts to all other transcripts")
table(rowSums(cpm == 0)==nsamples)[2]/sum(table(rowSums(cpm == 0)==nsamples))

# now set some cut-off to get rid of genes/transcripts with low counts
keepers <- filterByExpr(DGEList)
print(paste0("Transcripts left after filtering:", sum(keepers), 
             " which is ", round(sum(keepers)/length(keepers)*100), "% of the total ", length(keepers)))

DGEList.filtered <- DGEList[keepers,, keep.lib.sizes=F]

log2.cpm.filtered <- cpm(DGEList.filtered, log=TRUE)
boxplot(log2.cpm.filtered, las=2, cex=1, col = myColors, names = sampleLabels, 
        main="Filtered log2 cpm distribution")

# Normalize your data: the TMM normalization method has been found to perform well in comparative studies.
DGEList.filtered.norm <- calcNormFactors(DGEList.filtered, method = "TMM")

# use the 'cpm' function from EdgeR to get counts per million from you 

log2.cpm.filtered.norm <- cpm(DGEList.filtered.norm, log = TRUE, prior.count = 3)

boxplot(log2.cpm.filtered.norm, las=2, cex=1, col = myColors,  names = sampleLabels, 
        main="Filtered and normalized log2 cpm distribution")

# Hierarchical clustering ---------------
#hierarchical clustering can only work on a data matrix, not a data frame
#try using filtered and unfiltered data...how does this change the results

distance <- dist(t(log2.cpm), method="maximum") #other dist methods are "maximum", "manhattan", 
                                                # "canberra", "binary" or "minkowski"
clusters <- hclust(distance, method = "complete") #other methods are ward.D, ward.D2, single, complete, average
plot(clusters, labels=sampleLabels)

distance <- dist(t(log2.cpm.filtered), method="maximum") #other dist methods are "maximum", "manhattan", 
                                                          # "canberra", "binary" or "minkowski"
clusters <- hclust(distance, method = "complete") #other methods are ward.D, ward.D2, single, complete, average
plot(clusters, labels=sampleLabels)

distance <- dist(t(log2.cpm.filtered.norm), method="maximum") #other dist methods are "maximum", "manhattan", 
                                                              # "canberra", "binary" or "minkowski"
clusters <- hclust(distance, method = "complete") #other methods are ward.D, ward.D2, single, complete, average
plot(clusters, labels=sampleLabels)

# Pricipal component analysis (PCA) -------------

pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
save(pca.res, file = "pca.filtered.norm.Rda")
summary(pca.res) # Prints variance summary for all principal components.
pca.res$x #$x shows you how much each SAMPLE influenced each PC (called 'scores')

#lets first plot any two PCs against each other
data.frame <- as.data.frame(pca.res$x)  %>% rownames_to_column("samples")
ggplot(data.frame, aes(x=PC1, y=PC2, color=groups)) +
  geom_text(aes(x=PC1, y=PC2+15, label = samples)) +
  geom_point(size=5) +  theme(legend.position="right")

# Create a PCA 'small multiples' chart ----
# this is another way to view PCA to understand impact of each variable on each pricipal component
melted <- cbind(groups, melt(pca.res$x[,1:4]))
head(melted)
#look at your 'melted' data
ggplot(melted) +
  geom_bar(aes(x=Var1, y=value, fill=groups), stat="identity") + 
  theme(axis.text.x = element_text(angle = 45)) +
  facet_wrap(~Var2)
dev.off()

# =============================================================================
# Identify differentially expressed genes (DEGs) or transcripts (DETs)
# whether you look for differential expression at the gene or transcript level depends on how you read the
# Kallisto output into R using TxImport

# DE analysis using Limma/VOOM (alternatively, EdgeR or DESeq2) -----
# first create a DGEList object from your original count data using the DGEList function from EdgeR
DGEList.filtered.norm

#set up your design matrix
groups 

design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)

# We may look at the PCA or MDS of the normalized data:
pdf("2.Limma_voom.pdf")
{plotMDS(DGEList.filtered.norm, col = as.numeric(groups), main = "MDS plot of filtered normalized data",
        plot = T)
# normalize your data using the mean-variance relationship using the VOOM function from Limma
v.DEGList.filtered.norm <- voom(DGEList.filtered.norm, design, plot = TRUE)

# fit a linear model to your data
fit <- lmFit(v.DEGList.filtered.norm, design)
# Contrast matrix ----
contrast.matrix <- makeContrasts(status = disease - normal, levels=design)

# extract the linear model fit -----
fits <- contrasts.fit(fit, contrast.matrix)
#get bayesian stats for your linear model fit
ebFit <- eBayes(fits)
# TopTable to view DEGs -----
myTopHits_voom <- topTable(ebFit, adjust ="BH", number=7000, sort.by="logFC")
head(myTopHits_voom)
sum(myTopHits_voom$adj.P.Val<0.05)
}


{v.DEGList.filtered.norm <- voomWithQualityWeights(DGEList.filtered.norm, design, normalization = "none",
                                                   plot = TRUE)

# fit a linear model to your data
fit <- lmFit(v.DEGList.filtered.norm, design)
# Contrast matrix ----
contrast.matrix <- makeContrasts(status = disease - normal, levels=design)
# extract the linear model fit -----
fits <- contrasts.fit(fit, contrast.matrix)
#get bayesian stats for your linear model fit
ebFit <- eBayes(fits)
# TopTable to view DEGs -----
myTopHits <- topTable(ebFit, adjust ="BH", number=7000, sort.by="logFC")
head(myTopHits)
sum(myTopHits$adj.P.Val<0.05)

}


# decideTests to pull out the DEGs and make Venn Diagram ----
results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.05, lfc=1)

# take a look at what the results of decideTests looks like
tail(results)
summary(results)
vennDiagram(results, include="up")

# retrieve expression data for your DEGs ----
head(v.DEGList.filtered.norm$E) #numeric matrix of normalized expression values on the log2 scale
colnames(v.DEGList.filtered.norm$E) <- sampleLabels
diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,]
head(diffGenes)
dim(diffGenes)
#write your DEGs to a file
write_csv(as_tibble(diffGenes, rownames = "targetid"),
	 "DEA_transcript_level_LIMMA_significance.csv") 


# Volcano Plots ----------

myTopHits %>%
  mutate(significance = (adj.P.Val<0.05 & abs(logFC)>1) ) %>%
  ggplot() + geom_point(aes(x = logFC, y = -log10(adj.P.Val), color = significance)) + 
  ggtitle("Volcano plot - Limma")
dev.off()

# DE analysis using Sleuth ----

# Now you're ready to construct a sleuth object from Kallisto data:
mySleuth <- sleuth_prep(targets, ~ condition, target_mapping = Tx, read_bootstrap_tpm=TRUE, 
                        extra_bootstrap_summary=TRUE) 
models(mySleuth)
# Estimate parameters for the full and reduced models, to carry the likelhood ratio test (to tell transcripts 
# that are really affected by the condition)
mySleuth <- sleuth_fit(mySleuth); mySleuth <- sleuth_fit(mySleuth, ~1, 'reduced')
mySleuth <- sleuth_lrt(mySleuth, 'reduced', 'full')
results_lrt <- sleuth_results(mySleuth, 'reduced:full', 'lrt')
sig_results_lrt <- results_lrt %>% dplyr::arrange(qval) %>% dplyr::filter(qval < 0.05)

write_csv(sig_results_lrt, "DEA_transcript_level_sleuth_lrt_significance.csv")

models(mySleuth)
# use a wald test (WT) to also yield a beta statistic that approximates the fold change in expression
mySleuth <- sleuth_wt(mySleuth, which_beta = "conditionnormal")

# Now, export Sleuth results:
results_wt <- sleuth::sleuth_results(mySleuth, "conditionnormal")
sig_results_wt <- results_wt %>% dplyr::arrange(qval) %>% dplyr::filter(qval < 0.05) 

print("Minimum significant* effect size (Beta) in data:")
min(abs(sig_results_wt$b))

write_csv(sig_results_wt, "DEA_transcript_level_sleuth_wt_significance.csv") 

save(mySleuth, file = "mySlueth.so")

pdf("3.Sleuth_plots.pdf")
sleuth::plot_volcano(mySleuth, "conditionnormal")
sleuth::plot_ma(mySleuth, "conditionnormal")
plot_transcript_heatmap(mySleuth, sig_results_wt$target_id)
#plot_transcript_heatmap(mySleuth, sig_results_lrt$target_id)
plot_transcript_heatmap(mySleuth, sig_results_wt$target_id[1:10])
plot_vars(mySleuth, test_type = "wt")

p1 <- plot_fld(mySleuth, "sample1") + coord_cartesian(xlim = c(0,500), ylim = c(0, 0.035))   
p2 <- plot_fld(mySleuth, "sample2") + coord_cartesian(xlim = c(0,500), ylim = c(0, 0.035)) 
p3 <- plot_fld(mySleuth, "sample3") + coord_cartesian(xlim = c(0,500), ylim = c(0, 0.035))
p4 <- plot_fld(mySleuth, "sample4") + coord_cartesian(xlim = c(0,500), ylim = c(0, 0.035))
p5 <- plot_fld(mySleuth, "sample5") + coord_cartesian(xlim = c(0,500), ylim = c(0, 0.035))
p6 <- plot_fld(mySleuth, "sample6") + coord_cartesian(xlim = c(0,500), ylim = c(0, 0.035))

cowplot::plot_grid(p1, p2, p3, p4, p5, p6, ncol = 2, labels = paste0(LETTERS[1:6], ") sample", 1:6))

plot_sample_heatmap(mySleuth)
plot_pca(mySleuth, text_labels = TRUE, units = "est_counts", use_filtered = F)
plot_pc_variance(mySleuth, text_labels = T, use_filtered = T)

plot_group_density(mySleuth, use_filtered = TRUE, units = "est_counts",
                   trans = "log", grouping = "condition", offset = 1)

dev.off()


# =============================================================================
# create heatmaps from your differentially expressed genes or transcripts

pdf("4.heatmaps_cluster_summaries.pdf")

# choose color pallette ----
# Some useful examples: colorpanel(40, "darkblue", "yellow", "white"); heat.colors(75);
# cm.colors(75); rainbow(75); redgreen(75); library(RColorBrewer); rev(brewer.pal(9,"Blues")[-1]).

# a color-blind friendly pallete
myheatcol <- colorRampPalette(colors=c("yellow","white","blue"))(100)

# cluster DEGs ----
#begin by clustering the genes (rows) in each set of differentially expressed genes by pearson correlation
hr <- hclust(as.dist(1-cor(t(diffGenes), method="pearson")), method="complete") 

#now cluster your samples (columns) by spearman correlation
#we may not acutally use this clustering result, but it's good to have just in case
hc <- hclust(as.dist(1-cor(diffGenes, method="spearman")), method="complete") 
#note: we use Spearman, instead of Pearson, for clustering samples because it gives equal weight to highly vs 
# lowly expressed transcripts or genes

# Cut the resulting tree and create color vector for clusters.  
#Vary the cut height to give more or fewer clusters, or you the 'k' argument to force n number of clusters
#we'll look at these clusters in more detail later
mycl <- cutree(hr, k=2)

#now assign a color to each cluster (makes it easy to identify and manipulate)
mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9) 
mycolhc <- mycolhc[as.vector(mycl)] 

# produce heatmap of DEGs ----
#plot the hclust results as a heatmap

heatmap.2(diffGenes, Rowv=as.dendrogram(hr), Colv=T, 
          col=myheatcol, scale='row', 
          density.info="none", trace="none", RowSideColors=mycolhc, 
          cexRow=1, cexCol=1, margins=c(8,20)) 
#what do the colors represent in this heatmap?
#what happens when you change scale=NULL

# you can annotate samples with any metadata available in your study design file
color.map <- function(groups) { 
  if (groups=="disease") "#FF0000" else if (groups=="normal") "#33A12B" else "#0000FF"}
color.map <- unlist(lapply(groups, color.map))

heatmap.2(diffGenes, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),
          col=myheatcol, scale="row", 
          density.info="none", trace="none", 
          RowSideColors=mycolhc, ColSideColors = color.map,
          cexRow=1, cexCol=1, margins=c(8,20)) 


dev.off()

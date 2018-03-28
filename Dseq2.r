library(DESeq2)
library(stringr)
library(ggplot2)
counts=read.csv("/home/alnour/h3bionet/RNA_Assesment/cleaned_counts_2.csv", row.names="Geneid")
samples = colnames(counts)
condtions = append(rep("control",3), rep("normal", 3))
all(rownames(samples) %in% colnames(counts))
coldata = data.frame(condtions)
row.names(coldata) = samples
dds=DESeqDataSetFromMatrix(countData = as.matrix(counts), colData = coldata,design = ~ condtions)
dds <- DESeq(dds)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
res = results(dds)
resultsNames(dds)
summary(res)
resultsNames(dds)
summary(res)
resLFC <- lfcShrink(dds, coef=2)
resLFC

library(DESeq2)
library(stringr)
library(ggplot2)
file_path = "test/counts.txt"
counts=read.csv(file_path, sep="", head=T, skip=1, row.names = "Geneid")
colnames(counts)[6:11]
samples_list = c()
for(counter in 1:6){
    samples_list[counter] =  paste("sample", counter)
}
 colnames(counts)[6:11] = samples_list
 counts = counts[,6:11]
 
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

#!/usr/bin/env Rscript

library(VennDiagram)
library(UpSetR)
library(readr)
library(purrr)
library(tidyverse)

dirs <- list.files()
myfiles <- as_tibble(matrix(nrow = 4, ncol = 3)) %>% setNames(dirs)

data <- c("ballgown", "limma", "sleuth_lrt", "slueth_wt")
dets <- vector()
for (dir in dirs){
  print(dir)
  myfiles[, dir] = paste0(dir, "/", list.files(path = dirs[2], pattern = "DEA*"))
  mydir <- paste0(dir, "_",data)
  dets[mydir[1]] <- read_csv(myfiles[[1,dir]]) %>%
    mutate(transcriptName = str_replace(transcriptName, "[.].*", "")) 
  dets[mydir[2]] <- read_csv(myfiles[[2,dir]]) 
  dets[mydir[3]] <- read_csv(myfiles[[3,dir]])
  dets[mydir[4]] <- read_csv(myfiles[[4,dir]]) 
}

venn.plot <- venn.diagram(
  x = list(
    sleuth_wt = sleuth_wt$target_id, 
    sleuth_lrt = sleuth_lrt$target_id, 
    limma =  limma$targetid, 
    ballgown =  ballgown$transcriptName),
  filename = "common_transcripts.tiff",
  col = "transparent",
  fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
  alpha = 0.50,
  cex = 1.5,
  fontfamily = "serif",
  fontface = "bold",
  # cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
  cat.cex = 1.5,
  cat.pos = 0,
  cat.dist = 0.07,
  cat.fontfamily = "serif"
  # rotation.degree = 270,
  # margin = 0.2
);

aggregate_gene_list = list(
  sleuth_wt = sleuth_wt$target_id, 
  sleuth_lrt = sleuth_lrt$target_id, 
  limma =  limma$targetid, 
  ballgown =  ballgown$transcriptName)


aggregate_gene_list <- dets
upset(nsets = 12, fromList(aggregate_gene_list), order.by = c("freq"), keep.order = F,  
      sets.bar.color = "#56B4E9")


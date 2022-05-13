args <- commandArgs(trailingOnly = TRUE)
rna_file <- args[1]

# load necessary libraries
suppressPackageStartupMessages({
	library(dplyr)
	library(Seurat)
	library(mclust)
})

# load initial objects
rna <- readRDS(rna_file)

# F1 score
expected <- rna@meta.data$major_cell_type
predicted  <- rna@meta.data$predicted.id    
expected <- ifelse(expected %in% c("DN","DP","HSC","iNKT","PLT"),"Other",expected)
predicted <- ifelse(predicted %in% c("DN","Basophil","intMono"),"Other",predicted)
cm = as.matrix(table(expected, predicted))
precision <- diag(cm) / colSums(cm)
recall <- diag(cm) / rowSums(cm)
f1 <-  ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
print("The F1 score for the overlapping populations is:\n")
print(f1)

# ARI
ari <- adjustedRandIndex(rna@meta.data$major_cell_type, rna@meta.data$predicted.id)  
print("The ARI score for this integration is:\n")
print(ari)

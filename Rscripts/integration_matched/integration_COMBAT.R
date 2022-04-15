args <- commandArgs(trailingOnly = FALSE)
rna_file <- args[1]
adt_file <- args[2]
cytof_file <- args[3]
outfile <- args[4]

# load necessary libraries
suppressPackageStartupMessages({
	library(SingleCellExperiment)
	library(Seurat)
	library(dplyr)
	library(tidyr)
})

# load initial objects
rna <- readRDS(rna_file)
cytof <- readRDS(cytof_file)
adt <- readRDS(adt_file)

# set B cell markers and common features
DefaultAssay(rna) <- "RNA"
common_features <- intersect(row.names(rna), row.names(cytof))

# filtering cells without expression in common genes
indx <- colSums(rna[common_features,])>5
rna <- rna[,indx]
adt <- adt[,colnames(rna)]

# find anchors
transfer.anchors_cyt <- FindTransferAnchors(reference = cytof, query = rna, 
                                                  features = common_features,
                                                  reference.assay = "Cytof", 
                                                  query.assay = "RNA", 
                                                  reduction = "cca",
                                                  dims=1:15)

# transfer labels and protein intensities         
imputation_cyt <- TransferData(anchorset = transfer.anchors_cyt, 
                                 refdata = GetAssayData(cytof, assay = "Cytof", slot = "data")[common_features, ],
                                 weight.reduction = rna[["pca"]],
                                 dims=1:15)
rna[["Cytof"]] <- imputation_cyt

celltype.predictions_cyt <- TransferData(anchorset = transfer.anchors_cyt, 
                                               refdata = cytof@meta.data[,"fine_populations"],
                                               weight.reduction =rna[["pca"]],
                                               dims=1:15)
rna <- AddMetaData(rna, metadata = celltype.predictions_cyt[,c("predicted.id","prediction.score.max")])

saveRDS(rna, file = outfile)
sessionInfo()

args <- commandArgs(trailingOnly = FALSE)
rna_file <- args[1]
adt_file <- args[2]
cytof_file <- args[3]
k_param <- args[4]
outfile <- args[5]

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
Bplasmamarkers <- c("CD19","HLA-DR","CTLA4","CD28","Ki-67","CCR7","CD127","CD99","CD103","PD1","CD161","CD27",     
	"CD45RO","CD57","CD39","BCL-2","CD45RA","GZB","CD69","CD11c","CD20","IgD","KLRG1","FOXP3","CD38","CD25","CLA","CX3CR1" )
DefaultAssay(rna) <- "RNA"
common_features <- intersect(row.names(rna),Bplasmamarkers)

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
                                 dims=1:15, k.weight=k_param)
rna[["Cytof"]] <- imputation_cyt

celltype.predictions_cyt <- TransferData(anchorset = transfer.anchors_cyt, 
                                               refdata = cytof@meta.data[,"fine_populations"],
                                               weight.reduction =rna[["pca"]],
                                               dims=1:15, k.weight=k_param)
rna <- AddMetaData(rna, metadata = celltype.predictions_cyt[,c("predicted.id","prediction.score.max")])

saveRDS(rna, file = outfile)
sessionInfo()

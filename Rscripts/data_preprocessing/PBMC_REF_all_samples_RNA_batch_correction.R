### load necessary libraries
library(Seurat)
library(plyr)
library(tidyverse)
library(SingleCellExperiment)
library(future)

# load 10x 3' GEX Data and metadata from GEO GSE164378 for PBMC REF
rna_mtx <- Read10X(data.dir = "./PBMC_RNA_3P/")
pbmc_metadata <- read.table(file="./GSE164378_sc.meta.data_3P.csv", sep = ",", row.names = 1, header=T)


# Create Seurat object and Add metadata
rna_seurat <- CreateSeuratObject(counts = rna_mtx, project = "SeuratProject")
rna_seurat <- AddMetaData(rna_seurat, metadata = list(pbmc_metadata$nCount_RNA, pbmc_metadata$nFeature_RNA, pbmc_metadata$orig.ident, pbmc_metadata$lane, pbmc_metadata$donor), 
                              col.name = c("nCount_RNA","nFeature_RNA","orig.ident","lane","donor"))

rna_seurat <- AddMetaData(rna_seurat, metadata = list(pbmc_metadata$time, pbmc_metadata$celltype.l1, pbmc_metadata$celltype.l2, pbmc_metadata$celltype.l3, pbmc_metadata$Phase, pbmc_metadata$Batch), 
                              col.name = c("time","celltype.l1","celltype.l2","celltype.l3","Phase","Batch"))
rna_seurat$donor_time <- paste(rna_seurat$donor,"_",rna_seurat$time,sep = "")
rm(rna_mtx, pbmc_metadata)

#Subset seurat object to remove cells annotated as Doublet in the original publication
Idents(rna_seurat) <- "celltype.l2"

rna_seurat <- subset(rna_seurat, idents = "Doublet", invert=T)

# Run Seurat v3 RNA log normalize and RPCA anchor based integration workflow

# split samples by each donor and  log normalize and scale for nCount_RNA and Run PCA
Idents(rna_seurat) <- "donor_time"
rna_list <- SplitObject(rna_seurat, split.by = "donor_time")
rna_list <- lapply(X=rna_list, FUN = function(x){
  x <- NormalizeData(x,verbose=F)
  x <- FindVariableFeatures(x, nfeatures = 3000, verbose=F)
})


features <- SelectIntegrationFeatures(object.list = rna_list, nfeatures = 3000)

rna_list <- lapply(X = rna_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE, vars.to.regress = ("nCount_RNA"))
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# Use all donors_time0 (unvaccinated) samples the reference samples
plan("multisession", workers=4, gc=T)
options(future.globals.maxSize= 16*1024*1024^2)
anchors <- FindIntegrationAnchors(object.list = pbmc_3p_list, reference = c(8,9,11,12,14,16,18,23), reduction = "rpca",
                                     dims = 1:50, anchor.features= features)
saveRDS(anchors, file="./anchors_pbmc_multimodal_all_rna_workdlow.rds")

plan("sequential", gc=T)
rna.integrated <- IntegrateData(anchorset = anchors, dims = 1:50, normalization.method = "LogNormalize")

DefaultAssay(rna.integrated) <- "integrated"
plan("multisession", workers=4, gc=T)
options(future.globals.maxSize= 16*1024*1024^2)
rna.integrated <- ScaleData(rna.integrated, verbose = FALSE)

plan("sequential", gc=T)
rna.integrated <- RunPCA(rna.integrated, verbose = FALSE)
rna.integrated<- RunUMAP(rna.integrated, dims = 1:50)

saveRDS(rna.integrated, file="pbmc_multimodal_all_rna_workflow_rpca.rds")

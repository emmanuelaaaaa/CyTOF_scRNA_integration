### load necessary libraries
library(Seurat)
library(plyr)
library(tidyverse)
library(SingleCellExperiment)
library(future)
library(SeuratDisk)


# ADD ADT Assay for all samples
adt <- LoadH5Seurat(file="./pbmc_multimodal.h5seurat")

DefaultAssay(pbmc_ref) <- "ADT"
pbmc_ref[["SCT"]] <- NULL
pbmc_ref@reductions$pca <- "NULL"
pbmc_ref@reductions$spca <- "NULL"
pbmc_ref@reductions$umap <- "NULL"
pbmc_ref@reductions$wnn.umap <- "NULL"


# collapse ADT celltypes to match COMBAT CyTOF celltypes
adt$celltypes <- as.factor(adt$celltype.l2)
adt$celltypes <- fct_collapse( adt$celltypes ,B = c("B intermediate" ,"B memory" ,"B naive"), 
                               CD4 = c("CD4 CTL","CD4 Naive","CD4 Proliferating","CD4 TCM","CD4 TEM","Treg") , 
                               CD8 = c("CD8 Naive" ,"CD8 Proliferating","CD8 TCM","CD8 TEM"), NK = c( "NK Proliferating" ,"NK_CD56bright" , "NK"), 
                               DC = c("ASDC","cDC1" ,"cDC2" ,"pDC"), PB = c("Plasmablast"),
                               DN = c("dnT"), CD16_Mono = c("CD16 Mono"), CD14_Mono = c("CD14 Mono"), GDT = c("gdT"))
adt$donor_time <- paste(adt$donor,"_",adt$time,sep = "")

# remove doublets
Idents(adt) <- "celltype.l2"

adt <- subset(adt, idents = "Doublet", invert=T)

saveRDS(adt, file = "PBMC_REF_ADT_all_samples.rds")


# Subset ADT to only un-vaccinated samples
Idents(adt) <- "donor_time"

adt_day0 <- subset(x = adt, idents = c("P1_0", "P2_0","P3_0","P4_0","P5_0","P6_0","P7_0","P8_0"))
saveRDS(adt_day0, file="./pbmc_ref_adt_day0.rds")


# READ in CyTOF Toy dataset 182K mixed to match all samples PBMC-REF
cytof <- readRDS("./CyTOF_mixed/CyTOFobj_toy_mixed182000.RDS")
cytof$celltypes <- as.factor(cytof$major_cell_type)
cytof$celltypes <- fct_collapse( cytof$celltypes, CD16_Mono = c("ncMono"), CD14_Mono = c("cMono"), B = c("B cells"), DC = c("cDC1","pDC", "cDC2"), MAIT = c("MAIT cells"), PB = c("Plasmablasts"), CD4 = c("CD4 T cells"), GDT = c("Vd2 gd T cells"),
                                 CD8 = c("CD8 T cells"),DN = c("DN T cells"), NK = c("NK cells"))
saveRDS(cytof, file = "./CyTOF_mixed/CyTOFobj_toy_mixed182000.RDS")

# READ in CyTOF Toy dataset 50K Healthy controls to match unvaccinated PBMC-REF samples
cytof <- readRDS("./CyTOFobj_toy_50000.RDS")
cytof$celltypes <- as.factor(cytof$major_cell_type)
cytof$celltypes <- fct_collapse( cytof$celltypes, CD16_Mono = c("ncMono"), CD14_Mono = c("cMono"))

saveRDS(cytof, file = "./CyTOFobj_toy_50000.RDS" )

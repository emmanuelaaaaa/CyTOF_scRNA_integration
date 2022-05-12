library(SingleCellExperiment)
library(Seurat)
library(ggplot2)
library(cowplot)
library(plyr)
library(tidyverse)
library(future)
library(Matrix)


# Load the extracted bcells saved as a list of datatsets
seu_pbmc <- readRDS(file="./hdf5_files/covid_sc_pbmc_bcells_only_v2.rds")
seu_wb <- readRDS(file="./hdf5_files/covid_sc_wb_bcells_only_v2.rds")

# remove "schulte-schrepping_2020" from seu_pbmc list as all the Bcells for pbmc & WB present in the seu_wb list
seu_pbmc[["schulte-schrepping_2020"]] <- NULL

# combine seu_pbmc & seu_wb lists for easier integration
seuList <- c(seu_pbmc, seu_wb)
rm(seu_pbmc,seu_wb)

# we also remove the kusanadi dataset as it only has less than 100 Bcells and this will not  be helpful for us and it is causing #integration errors based on these github posts: https://github.com/satijalab/seurat/issues/4803 and #https://github.com/satijalab/seurat/issues/3930
seuList[["kusnadi_2021"]] <- NULL

# We decided to use SU 2020 on its own as it had ADT data for independent validation of the COMBAT dataset
seuList[["su_2020"]] <- NULL

#Remove non-Bcells using paper criteria
#use the criteria from this paper: https://www.frontiersin.org/articles/10.3389/fimmu.2021.602539/full#h3

#criteria: Cells were selected for further analyses according to the following criteria: (i) express zero CD3E, GNLY, CD14, FCER1A, GCGR3A or FCGR3A , LYZ, PPBP and CD8A transcripts, to exclude any non-B cells

#As some genes were not present in the different data sets, we had to do this individually for each dataset in the list

DefaultAssay(seuList[[1]]) <- "RNA"
seuList_QC[[1]] <- subset(seuList[[1]],  slot= "counts", subset= CD3E == 0 & GNLY ==0 & CD14==0 & FCER1A==0 & FCGR3A ==0 & LYZ ==0 & PPBP ==0 & CD8A ==0)

DefaultAssay(seuList[[2]]) <- "RNA"
seuList_QC[[2]] <- subset(seuList[[2]],  slot= "counts", subset= CD3E == 0 & GNLY ==0 & CD14==0 & FCER1A==0 & FCGR3A ==0 & LYZ ==0 & PPBP ==0 & CD8A ==0)

DefaultAssay(seuList[[3]]) <- "RNA"
seuList_QC[[3]] <- subset(seuList[[3]],  slot= "counts", subset= CD3E == 0 & GNLY ==0 & CD14==0 & FCER1A==0 & FCGR3A ==0 & LYZ ==0 & PPBP ==0 & CD8A ==0)

DefaultAssay(seuList[[4]]) <- "RNA"
seuList_QC[[4]] <- subset(seuList[[4]],  slot= "counts", subset= CD3E == 0 & GNLY ==0 & CD14==0 & FCER1A==0 & FCGR3A ==0 & LYZ ==0 & PPBP ==0 & CD8A ==0)


DefaultAssay(seuList[[5]]) <- "RNA"
seuList_QC[[5]] <- subset(seuList[[5]],  slot= "counts", subset= CD3E == 0 & GNLY ==0 & CD14==0 & FCER1A==0 & FCGR3A ==0 & LYZ ==0 & PPBP ==0 & CD8A ==0)

DefaultAssay(seuList[[6]]) <- "RNA"
seuList_QC[[6]] <- subset(seuList[[6]],  slot= "counts", subset= CD3E == 0 & GNLY ==0 & CD14==0 & FCER1A==0 & FCGR3A ==0 & LYZ ==0 & PPBP ==0 & CD8A ==0)

DefaultAssay(seuList[[7]]) <- "RNA"
seuList_QC[[7]] <- subset(seuList[[7]],  slot= "counts", subset= CD3E == 0 & GNLY ==0 & CD14==0 & FCGR3A ==0 & LYZ ==0  & CD8A ==0)

DefaultAssay(seuList[[8]]) <- "RNA"
seuList_QC[[8]] <- subset(seuList[[8]],  slot= "counts", subset= CD3E == 0 & GNLY ==0 & CD14==0 & FCGR3A ==0 & LYZ ==0 & PPBP ==0 & CD8A ==0)

DefaultAssay(seuList[[9]]) <- "RNA"
seuList_QC[[9]] <- subset(seuList[[9]],  slot= "counts", subset= CD3E == 0 & GNLY ==0 & CD14==0 & FCER1A==0 & FCGR3A ==0 & LYZ ==0 & PPBP ==0 & CD8A ==0)

DefaultAssay(seuList[[10]]) <- "RNA"
seuList_QC[[10]] <- subset(seuList[[10]],  slot= "counts", subset= CD3E == 0 & GNLY ==0 & CD14==0 & FCER1A==0 & FCGR3A ==0 & LYZ ==0 & PPBP ==0 & CD8A ==0)

DefaultAssay(seuList[[11]]) <- "RNA"
seuList_QC[[11]] <- subset(seuList[[11]],  slot= "counts", subset= CD3E == 0 & GNLY ==0 & CD14==0 & FCER1A==0 & FCGR3A ==0 & LYZ ==0 & PPBP ==0 & CD8A ==0)

DefaultAssay(seuList[[12]]) <- "RNA"
seuList_QC[[12]] <- subset(seuList[[12]],  slot= "counts", subset= CD3E == 0 & GNLY ==0 & CD14==0 & FCER1A==0 & FCGR3A ==0 & LYZ ==0 & PPBP ==0 & CD8A ==0)


DefaultAssay(seuList[[13]]) <- "RNA"
seuList_QC[[13]] <- subset(seuList[[13]],  slot= "counts", subset= CD3E == 0 & GNLY ==0 & CD14==0 & FCER1A==0 & FCGR3A ==0 & LYZ ==0 & PPBP ==0 & CD8A ==0)

DefaultAssay(seuList[[14]]) <- "RNA"
seuList_QC[[14]] <- subset(seuList[[14]],  slot= "counts", subset= CD3E == 0 & GNLY ==0 & CD14==0 & FCER1A==0 & FCGR3A ==0 & LYZ ==0 & PPBP ==0 & CD8A ==0)


# remove cells  predicted  as plasmablasts


seuList_QC <- lapply(X=seuList_QC, FUN = function(x){
  Idents(x) <- "predicted.celltype.l2"
  x <- subset(x, idents = c("Plasmablast"), invert=T)
})

seuList_QC[["meckiff_2020"]] <- NULL
seuList_QC[["bacher_2020"]] <- NULL

saveRDS(seuList_QC, file="./hdf5_files/covid_sc_pbmc_bcells_only_v3.rds")




# Run Seurat RNA log normalise and RPCA integration workflow 
seuList <- readRDS(file="./hdf5_files/covid_sc_pbmc_bcells_only_v3.rds")

# Change Gene names to CyTOF Protein names in all the relevant RNA assay slots
gene_protein_table <- read.table(file="./cytoF_protein_gene_table_bcells.txt",header=T)

# wilk 2020
gene_protein_table_s <- filter(gene_protein_table, gene_name %in% row.names(seuList$wilk_2020@assays$RNA))
# count slot
m <- match(gene_protein_table_s$gene_name,row.names(seuList$wilk_2020@assays$RNA@counts))
row.names(seuList$wilk_2020@assays$RNA@counts)[m] <- gene_protein_table_s$cytof_protein
#data slot
m <- match(gene_protein_table_s$gene_name,row.names(seuList$wilk_2020@assays$RNA@data))
row.names(seuList$wilk_2020@assays$RNA@data)[m] <- gene_protein_table_s$cytof_protein

# wen 2020
#count slot
gene_protein_table_s <- filter(gene_protein_table, gene_name %in% row.names(seuList$wen_2020@assays$RNA))
m <- match(gene_protein_table_s$gene_name,row.names(seuList$wen_2020@assays$RNA@counts))
row.names(seuList$wen_2020@assays$RNA@counts)[m] <- gene_protein_table_s$cytof_protein
#data slot
m <- match(gene_protein_table_s$gene_name,row.names(seuList$wen_2020@assays$RNA@data))
row.names(seuList$wen_2020@assays$RNA@data)[m] <- gene_protein_table_s$cytof_protein

# lee 2020
#counts slot
gene_protein_table_s <- filter(gene_protein_table, gene_name %in% row.names(seuList$lee_2020@assays$RNA))
m <- match(gene_protein_table_s$gene_name,row.names(seuList$lee_2020@assays$RNA@counts))
row.names(seuList$lee_2020@assays$RNA@counts)[m] <- gene_protein_table_s$cytof_protein
# data slot
m <- match(gene_protein_table_s$gene_name,row.names(seuList$lee_2020@assays$RNA@data))
row.names(seuList$lee_2020@assays$RNA@data)[m] <- gene_protein_table_s$cytof_protein

# zhu 2020
#counts slot
gene_protein_table_s <- filter(gene_protein_table, gene_name %in% row.names(seuList$zhu_2020@assays$RNA))
m <- match(gene_protein_table_s$gene_name,row.names(seuList$zhu_2020@assays$RNA@counts))
row.names(seuList$zhu_2020@assays$RNA@counts)[m] <- gene_protein_table_s$cytof_protein
# data slot
m <- match(gene_protein_table_s$gene_name,row.names(seuList$zhu_2020@assays$RNA@data))
row.names(seuList$zhu_2020@assays$RNA@data)[m] <- gene_protein_table_s$cytof_protein
  

#arunachalam_2020
#counts slot
gene_protein_table_s <- filter(gene_protein_table, gene_name %in% row.names(seuList$arunachalam_2020@assays$RNA))
m <- match(gene_protein_table_s$gene_name,row.names(seuList$arunachalam_2020@assays$RNA@counts))
row.names(seuList$arunachalam_2020@assays$RNA@counts)[m] <- gene_protein_table_s$cytof_protein
# data slot
m <- match(gene_protein_table_s$gene_name,row.names(seuList$arunachalam_2020@assays$RNA@data))
row.names(seuList$arunachalam_2020@assays$RNA@data)[m] <- gene_protein_table_s$cytof_protein

#yu_2020
#counts slot
gene_protein_table_s <- filter(gene_protein_table, gene_name %in% row.names(seuList$yu_2020@assays$RNA))
m <- match(gene_protein_table_s$gene_name,row.names(seuList$yu_2020@assays$RNA@counts))
row.names(seuList$yu_2020@assays$RNA@counts)[m] <- gene_protein_table_s$cytof_protein
# data slot
m <- match(gene_protein_table_s$gene_name,row.names(seuList$yu_2020@assays$RNA@data))
row.names(seuList$yu_2020@assays$RNA@data)[m] <- gene_protein_table_s$cytof_protein

#stephenson_2021
#counts slot
gene_protein_table_s <- filter(gene_protein_table, gene_name %in% row.names(seuList$stephenson_2021@assays$RNA))
m <- match(gene_protein_table_s$gene_name,row.names(seuList$stephenson_2021@assays$RNA@counts))
row.names(seuList$stephenson_2021@assays$RNA@counts)[m] <- gene_protein_table_s$cytof_protein
# data slot
m <- match(gene_protein_table_s$gene_name,row.names(seuList$stephenson_2021@assays$RNA@data))
row.names(seuList$stephenson_2021@assays$RNA@data)[m] <- gene_protein_table_s$cytof_protein

#yao_2021
#counts slot
gene_protein_table_s <- filter(gene_protein_table, gene_name %in% row.names(seuList$yao_2021@assays$RNA))
m <- match(gene_protein_table_s$gene_name,row.names(seuList$yao_2021@assays$RNA@counts))
row.names(seuList$yao_2021@assays$RNA@counts)[m] <- gene_protein_table_s$cytof_protein
# data slot
m <- match(gene_protein_table_s$gene_name,row.names(seuList$yao_2021@assays$RNA@data))
row.names(seuList$yao_2021@assays$RNA@data)[m] <- gene_protein_table_s$cytof_protein

#silvin_2020
#counts slot
gene_protein_table_s <- filter(gene_protein_table, gene_name %in% row.names(seuList$silvin_2020@assays$RNA))
m <- match(gene_protein_table_s$gene_name,row.names(seuList$silvin_2020@assays$RNA@counts))
row.names(seuList$silvin_2020@assays$RNA@counts)[m] <- gene_protein_table_s$cytof_protein
# data slot
m <- match(gene_protein_table_s$gene_name,row.names(seuList$silvin_2020@assays$RNA@data))
row.names(seuList$silvin_2020@assays$RNA@data)[m] <- gene_protein_table_s$cytof_protein


#`schulte-schrepping_2020`
#counts slot
gene_protein_table_s <- filter(gene_protein_table, gene_name %in% row.names(seuList$`schulte-schrepping_2020`@assays$RNA))
m <- match(gene_protein_table_s$gene_name,row.names(seuList$`schulte-schrepping_2020`@assays$RNA@counts))
row.names(seuList$`schulte-schrepping_2020`@assays$RNA@counts)[m] <- gene_protein_table_s$cytof_protein
# data slot
m <- match(gene_protein_table_s$gene_name,row.names(seuList$`schulte-schrepping_2020`@assays$RNA@data))
row.names(seuList$`schulte-schrepping_2020`@assays$RNA@data)[m] <- gene_protein_table_s$cytof_protein

#combes_2021
#counts slot
gene_protein_table_s <- filter(gene_protein_table, gene_name %in% row.names(seuList$combes_2021@assays$RNA))
m <- match(gene_protein_table_s$gene_name,row.names(seuList$combes_2021@assays$RNA@counts))
row.names(seuList$combes_2021@assays$RNA@counts)[m] <- gene_protein_table_s$cytof_protein
# data slot
m <- match(gene_protein_table_s$gene_name,row.names(seuList$combes_2021@assays$RNA@data))
row.names(seuList$combes_2021@assays$RNA@data)[m] <- gene_protein_table_s$cytof_protein

#bost_2021
#counts slot
gene_protein_table_s <- filter(gene_protein_table, gene_name %in% row.names(seuList$bost_2021@assays$RNA))
m <- match(gene_protein_table_s$gene_name,row.names(seuList$bost_2021@assays$RNA@counts))
row.names(seuList$bost_2021@assays$RNA@counts)[m] <- gene_protein_table_s$cytof_protein
# data slot
m <- match(gene_protein_table_s$gene_name,row.names(seuList$bost_2021@assays$RNA@data))
row.names(seuList$bost_2021@assays$RNA@data)[m] <- gene_protein_table_s$cytof_protein




seuList <- lapply(X=seuList, FUN = function(x){
  x <- NormalizeData(x,verbose=F)
  x <- FindVariableFeatures(x, nfeatures = 2000, verbose=F)
})

features <- SelectIntegrationFeatures(object.list = seuList,nfeatures = 2000 )
seuList <- lapply(X = seuList, FUN = function(x) {
  x$sex <- as.factor(x$sex)
  if(length(levels(x$sex)) > 1){
    message(">> ScaleData while regressing nUMI,sample and sex....")
    x <- ScaleData(x, features = features, verbose = FALSE, vars.to.regress = c("nCount_RNA","sample","sex"))
    
  } else {
    message(">> ScaleData while regressing nUMI & sample")
    x <- ScaleData(x, features = features, verbose = FALSE, vars.to.regress = c("nCount_RNA","sample"))
    
  }
  
  # x <- RunPCA(x, features = features, verbose = FALSE, reduction.name="pca")
  #x <- RunUMAP(x ,dims = 1:15, reduction.name ="umap")
})

invisible(gc())

plan("multisession", workers=3, gc=TRUE)

options(future.globals.maxSize= 15*1024*1024^2)
features <- SelectIntegrationFeatures(object.list = seuList,nfeatures = 1000 )

anchors <- FindIntegrationAnchors(object.list = seuList, reduction = "cca",
                                  dims = 1:15, anchor.features= features, scale=FALSE)
saveRDS(anchors, file="./hdf5_files/anchors_covid_pbmc_ref_bcells_rna_workflow_cca_v6.rds")

rm(seuList)
plan("sequential")

pbmc_covid.integrated <- IntegrateData(anchorset = anchors, dims = 1:15, normalization.method = "LogNormalize")
invisible(gc())


saveRDS(pbmc_covid.integrated, file="./hdf5_files/pbmc_covid_bcells_NO_PB_rna_workflow_cca_v6.rds")

rm(anchors)
invisible(gc())

DefaultAssay(pbmc_covid.integrated) <- "integrated"
plan("multisession", workers=3, gc=TRUE)
options(future.globals.maxSize= 16*1024*1024^2)

pbmc_covid.integrated <- ScaleData(pbmc_covid.integrated,verbose = FALSE)

plan("sequential")
pbmc_covid.integrated <- RunPCA(pbmc_covid.integrated, verbose = FALSE,npcs = 20,reduction.name="pca")
pbmc_covid.integrated<- RunUMAP(pbmc_covid.integrated, dims = 1:15, reduction = "pca", reduction.name ="umap")
invisible(gc())
saveRDS(pbmc_covid.integrated, file="./hdf5_files/pbmc_covid_bcells_NO_PB_rna_workflow_cca_v6.rds")


## Subset to Covid only cells and re Scale, PCA and UMAP on the integrated assay

plan("sequential")
DefaultAssay(pbmc_covid.integrated) <- "integrated"
pbmc_covid.integrated$disease_combined <- paste(pbmc_covid.integrated$disease_status,"_",pbmc_covid.integrated$disease_severity,sep="")
Idents(pbmc_covid.integrated) <- "disease_combined"
pbmc_covid.integrated <- subset(pbmc_covid.integrated, idents =c("healthy_","other_"), invert=T)


options(future.globals.maxSize= 10*1024*1024^2)

pbmc_covid.integrated <- ScaleData(pbmc_covid.integrated,verbose = FALSE, vars.to.regress = c("sex"))

plan("sequential")
pbmc_covid.integrated <- RunPCA(pbmc_covid.integrated, verbose = FALSE,npcs = 20,reduction.name="pca")
pbmc_covid.integrated<- RunUMAP(pbmc_covid.integrated, dims = 1:15, reduction = "pca", reduction.name ="umap")
invisible(gc())


saveRDS(pbmc_covid.integrated, file="./hdf5_files/pbmc_covid_bcells_No_PB_rna_workflow_cca_v6_covid_only.rds")




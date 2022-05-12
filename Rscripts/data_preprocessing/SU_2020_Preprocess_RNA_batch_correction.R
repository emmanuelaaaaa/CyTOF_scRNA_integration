library(SingleCellExperiment)
library(Seurat)
library(plyr)
library(tidyverse)
library(future)
library(Matrix)

# Load Su_2020 B-cell subsetted object only 
su_2020 <- readRDS(file="./hdf5_files/su_2020_bcell_processed.rds")
su_2020[["predicted_ADT"]] <- NULL
su_2020[["combined_ADT"]] <- NULL
su_2020@reductions$wnn.umap <- NULL
su_2020@reductions$adt.umap <- NULL
su_2020@reductions$adt.pca <- NULL
su_2020@reductions$adt.umap <- NULL
su_2020@reductions$rna.umap <- NULL


# Remove non-Bcells using paper criteria
#use the criteria from this paper: https://www.frontiersin.org/articles/10.3389/fimmu.2021.602539/full#h3

#criteria: Cells were selected for further analyses according to the following criteria: (i) express zero CD3E, GNLY, CD14, FCER1A, GCGR3A or FCGR3A , LYZ, PPBP and CD8A transcripts, to exclude any non-B cells

#remove cells annotated as  Plasmablasts in "predicted.celltype.l2" annotations
DefaultAssay(su_2020) <- "RNA"
su_2020 <- subset(su_2020,  slot= "counts", subset= CD3E == 0 & GNLY ==0 & CD14==0 & FCER1A==0 & FCGR3A ==0 & LYZ ==0 & PPBP ==0 & CD8A ==0)
Idents(su_2020) <- "predicted.celltype.l2"
su_2020 <- subset(su_2020, idents = c("Plasmablast"), invert=T)



# Change Gene names to CyTOF Protein names in all the relevant RNA assay slots
gene_protein_table <- read.table(file="./cytoF_protein_gene_table_bcells.txt",header=T)


### renaming the genes into their protein counterparts (from the gene-ADT-CyTOF protein table)
gene_protein_table <- filter(gene_protein_table, gene_name %in% row.names(su_2020@assays$RNA))
#gene_protein_table$protein_name <- gsub("-",replacement="_",gene_protein_table$protein_name)
m <- match(gene_protein_table$gene_name,row.names(su_2020@assays$RNA@counts))
row.names(su_2020@assays$RNA@counts)[m] <- gene_protein_table$cytof_protein

m <- match(gene_protein_table$gene_name,row.names(su_2020@assays$RNA@data))
row.names(su_2020@assays$RNA@data)[m] <- gene_protein_table$cytof_protein


# Split the Su 2020 dataset by batch and then do RNA based seurat v3 batch correction
DefaultAssay(su_2020) <- "RNA"
seuList <- SplitObject(su_2020, split.by = "batch")

# Run Seurat RNA log normalise and CCA integration workflow 
seuList <- lapply(X=seuList, FUN = function(x){
  x <- NormalizeData(x,verbose=F)
  x <- FindVariableFeatures(x, nfeatures = 2000, verbose=F)
})

features <- SelectIntegrationFeatures(object.list = seuList,nfeatures = 1000 )
seuList <- lapply(X = seuList, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE, vars.to.regress = c("nCount_RNA","sex_standard"))
  #x <- RunPCA(x, features = features, verbose = FALSE)
  # x <- RunUMAP(x ,dims = 1:15)
})

invisible(gc())

plan("multisession", workers=3, gc=TRUE)

options(future.globals.maxSize= 15*1024*1024^2)

features <- SelectIntegrationFeatures(object.list = seuList,nfeatures = 1000 )

anchors <- FindIntegrationAnchors(object.list = seuList, reduction = "cca",
                                  dims = 1:15, anchor.features= features, scale = FALSE)
saveRDS(anchors, file="./anchors_su_covid_pbmc_ref_bcells_No_PB_rna_workflow_cca.rds")

rm(seuList)
plan("sequential")

pbmc_covid.integrated <- IntegrateData(anchorset = anchors, dims = 1:15, normalization.method = "LogNormalize")
invisible(gc())


saveRDS(pbmc_covid.integrated, file="su_int_pbmc_covid_bcells_No_PB_rna_workflow_cca.rds")

rm(anchors)

DefaultAssay(pbmc_covid.integrated) <- "integrated"
plan("multisession", workers=3, gc=TRUE)
options(future.globals.maxSize= 15*1024*1024^2)

pbmc_covid.integrated <- ScaleData(pbmc_covid.integrated, verbose = FALSE, vars.to.regress="sex_standard")

plan("sequential")
pbmc_covid.integrated <- RunPCA(pbmc_covid.integrated, verbose = FALSE,npcs = 20, reduction.name="pca")
pbmc_covid.integrated<- RunUMAP(pbmc_covid.integrated, dims = 1:15,reduction = "pca", reduction.name ="umap")
invisible(gc())

# Change ADT names to match CyTOF protein name
adt_cytof_table <- read.table(file="./Su_2020_adt_protein_combat_cytof_protein_table.txt",header=T)


### renaming the genes into their protein counterparts (from the gene-ADT-CyTOF protein table)
adt_cytof_table <- filter(adt_cytof_table, adt_protein %in% row.names(pbmc_covid.integrated@assays$ADT))
m <- match(adt_cytof_table$adt_protein,row.names(pbmc_covid.integrated@assays$ADT@counts))
row.names(pbmc_covid.integrated@assays$ADT@counts)[m] <- adt_cytof_table$cytof_name


# Normalize and SCALE ADT assay

DefaultAssay(pbmc_covid.integrated) <- 'ADT'

VariableFeatures(pbmc_covid.integrated) <- rownames(pbmc_covid.integrated[["ADT"]])
pbmc_covid.integrated <- NormalizeData(pbmc_covid.integrated, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')


saveRDS(pbmc_covid.integrated, file="./su_int_pbmc_covid_bcells_No_PB_rna_workflow_cca.rds")



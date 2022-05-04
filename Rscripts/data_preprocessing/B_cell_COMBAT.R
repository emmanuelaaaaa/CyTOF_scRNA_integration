args <- commandArgs(trailingOnly = FALSE)
rna_file <- args[1]
rna_md_file <- args[2]
adt_file <- args[3]
cytof_file <- args[4]
outdir <- args[5]

# load necessary libraries
suppressPackageStartupMessages({
	library(SingleCellExperiment)
	library(dplyr)
	library(ggplot2)
	library(Seurat)
	library(cowplot)
})

# set B cell markers
Bmarkers <- c("CD19","HLA-DR","CTLA4","CD28","Ki-67","CCR7","CD127","CD99","CD103","PD1","CD161","CD27",     
	 "CD45RO","CD57","CD39","BCL-2","CD45RA","GZB","CD69","CD11c","CD20","IgD","KLRG1","FOXP3","CD38","CD25","CLA","CX3CR1" )

# load initial objects
rna_mtx <- readRDS(rna_file)
rna_md <- readRDS(rna_md_file)
adt_obj <- readRDS(adt_file)
sce <- readRDS(cytof_file)

### CyTOF
## downsampling, converting sce object into Seurat and processing

# subsample cytof by sample_ID and extract overlapping samples with scRNA
sce_coldata_subs_smaller <- as_tibble(colData(sce)) %>% filter(major_cell_type=="B cells" & condition!="Sepsis") %>% 
                            group_by(patient_id)  %>% slice_sample(n=700) 
sce_coldata_subs_smaller <- as.data.frame(filter(sce_coldata_subs_smaller, COMBAT_ID_Time %in% rna_md$COMBAT_ID_Time))
row.names(sce_coldata_subs_smaller) <- sce_coldata_subs_smaller$cellID
sce_subs <- sce[,colData(sce)$cellID %in% sce_coldata_subs_smaller$cellID]

cytofdata <- assay(sce_subs)
colnames(cytofdata) <- colData(sce_subs)$cellID 

# renaming IgD_TCRgd into IgD so that it can be used for the B cell analysis
row.names(cytofdata)[row.names(cytofdata)=="IgD_TCRgd"] <- "IgD"

# Creating the Seurat obj
cytof_obj <- CreateSeuratObject(counts = cytofdata, assay = "Cytof", meta.data=sce_coldata_subs_smaller, project = "CyTOF_data")
Idents(cytof_obj) <- cytof_obj@meta.data$fine_populations
cytof_obj <- ScaleData(cytof_obj, features=Bmarkers)
cytof_obj <- RunPCA(cytof_obj, features=Bmarkers)
cytof_obj <- RunUMAP(cytof_obj, reduction = "pca", features=Bmarkers, dims = NULL)
cytof_obj <- RunTSNE(cytof_obj, reduction = "pca", features=Bmarkers, dims = NULL, do.fast = TRUE, check_duplicates = FALSE)

saveRDS(cytof_obj, file=paste0(outdir,"/COMBAT_CyTOFobj_Bcells.RDS"))

### RNA
## subsetting the original matrix for B cells only, creating the Seurat object, batch correcting and processing

cell_table_subs_smaller <- as.data.frame(rna_md %>% filter(major_cell_type=="B") )
cell_table_subs_smaller <- as.data.frame(filter(cell_table_subs_smaller, COMBAT_ID_Time %in% sce_coldata_subs_smaller$COMBAT_ID_Time))
row.names(cell_table_subs_smaller) <- cell_table_subs_smaller$barcode_id
rna_mtx_subs_smaller <- rna_mtx[,row.names(cell_table_subs_smaller)]

# Creating the Seurat obj
rna_obj <- CreateSeuratObject(rna_mtx_subs_smaller, assay = "RNA", meta.data=cell_table_subs_smaller, project = "10X_RNA")

# batch correction
Idents(rna_obj) <- "COMBAT_ID_Time"
rna_list <- SplitObject(rna_obj, split.by = "COMBAT_ID_Time")
rna_list <- lapply(X=rna_list, FUN = function(x){
    x <- NormalizeData(x,verbose=F)
    x <- FindVariableFeatures(x, nfeatures = 1000, verbose=F)
})

features <- SelectIntegrationFeatures(object.list = rna_list, features=1000)
rna_list_all <- rna_list
rna_list <- rna_list_all[unlist(lapply(rna_list_all, ncol)>80)]
rna_list <- lapply(X = rna_list, FUN = function(x) {
	 x <- ScaleData(x, features = features, verbose = FALSE, vars.to.regress = ("nCount_RNA"))
	 x <- RunPCA(x, features = features, verbose = FALSE, npcs=30)
})

# use the HC as reference
HC_samples <- colData(sce) %>% as_tibble () %>% filter(condition=="HV") %>% distinct(COMBAT_ID_Time)
rna_anchors <- FindIntegrationAnchors(object.list = rna_list, reference=which(names(rna_list) %in% HC_samples$COMBAT_ID_Time), 
   	    dims = 1:30, scale = F, reduction = "rpca")

all_features <- lapply(rna_list, row.names) %>% Reduce(intersect, .) 

# integrate
rna.integrated <- IntegrateData(anchorset = rna_anchors, dims = 1:50, k.weight = 80, normalization.method = "LogNormalize") 
rna.integrated <- ScaleData(rna.integrated, verbose = FALSE)
rna.integrated <- RunPCA(rna.integrated, verbose = FALSE)
rna.integrated <- RunUMAP(rna.integrated, dims = 1:15)

# also running scaling, PCA and UMAP on the uncorrected object
DefaultAssay(rna.integrated) <- "RNA"
rna.integrated <- ScaleData(rna.integrated, features = row.names(rna.integrated))
rna.integrated <- RunPCA(rna.integrated, features = row.names(rna.integrated),reduction.name = "rna.pca" )
rna.integrated <- RunUMAP(rna.integrated, dims = 1:15, reduction = "rna.pca",assay = "RNA", reduction.name = "rna.umap", reduction.key = "rnaUMAP_")

saveRDS(rna.integrated, file=paste0(outdir,"/COMBAT_RNAobj_Bcells.RDS"))

### ADT
## subsetting Seurat object for B cells only and processing

adt_obj <- adt_obj[,colnames(rna.integrated)]
adt_obj <- ScaleData(adt_obj, features = row.names(adt_obj))
adt_obj <- RunPCA(adt_obj, features = row.names(adt_obj))
adt_obj <- RunTSNE(object = adt_obj, dims.use = 1:30, do.fast = TRUE)
adt_obj <- RunUMAP(adt_obj, dims = 1:30)

saveRDS(adt_obj, file=paste0(outdir,"/COMBAT_ADTobj_Bcells.RDS"))

date()
sessionInfo()

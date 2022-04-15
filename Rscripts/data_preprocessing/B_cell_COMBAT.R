args <- commandArgs(trailingOnly = FALSE)
rna <- args[1]
adt <- args[2]
cytof <- args[3]

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
load("/t1-data/project/taylorlab/erepapi/Fellowship/COVID19_Integration/RData/CyTOF_ADT_initialobj_newannot.RData")
load("/t1-data/project/taylorlab/erepapi/Fellowship/COVID19_Integration/RData/CyTOF_ADT_subsets_newannot.RData") # this has the right tables where low quality samples have been removed for the ADT
rm(cell_table_subs_allgroups) # won't need that - to avoid confusion

# subsample cytof by sample_ID and extract overlapping samples with scRNA
sce_coldata_subs_smaller <- as_tibble(colData(sce)) %>% filter(major_cell_type=="B cells" & condition!="Sepsis") %>% 
                            group_by(patient_id)  %>% slice_sample(n=700) 
sce_coldata_subs_smaller <- as.data.frame(filter(sce_coldata_subs_smaller, COMBAT_ID_Time %in% cell_table_subs$COMBAT_ID_Time))
row.names(sce_coldata_subs_smaller) <- sce_coldata_subs_smaller$cellID
sce_subs <- sce[,colData(sce)$cellID %in% sce_coldata_subs_smaller$cellID]

# convert sce object into Seurat and process
cytofdata <- assay(sce_subs)
colnames(cytofdata) <- colData(sce_subs)$cellID 
row.names(cytofdata)[row.names(cytofdata)=="IgD_TCRgd"] <- "IgD"
cytof_obj <- CreateSeuratObject(counts = cytofdata, assay = "Cytof", meta.data=sce_coldata_subs_smaller, project = "CyTOF_data")
Idents(cytof_obj) <- cytof_obj@meta.data$fine_populations
cytof_obj <- ScaleData(cytof_obj, features=Bmarkers)
cytof_obj <- RunPCA(cytof_obj, features=Bmarkers)
cytof_obj <- RunUMAP(cytof_obj, reduction = "pca", features=Bmarkers, dims = NULL)
cytof_obj <- RunTSNE(cytof_obj, reduction = "pca", features=Bmarkers, dims = NULL, do.fast = TRUE, check_duplicates = FALSE)

saveRDS(cytof_obj, file="/t1-data/project/taylorlab/erepapi/Fellowship/COVID19_Integration/RData/CyTOFobj_toy_Bcells_55000_nosepsis.RDS")

####RNA
rna_obj <- readRDS("/t1-data/project/taylorlab/erepapi/Fellowship/COVID19_Integration/RData/RNAobjintegrated_toy_Bcells.RDS")
rna_obj <- subset(rna_obj, cells=colnames(rna_obj)[rna_obj@meta.data$COMBAT_ID_Time %in% cytof_obj@meta.data$COMBAT_ID_Time])
setdiff(cytof_obj@meta.data$COMBAT_ID_Time,rna_obj@meta.data$COMBAT_ID_Time)
[1] "S00056-Ja005" "S00099-Ja005"
# because in the rna object I had already done the filtering to remove samples with low numbers of cells for the integration
# 32728 cells remaining

Idents(rna_obj) <- "COMBAT_ID_Time"

rna_list <- SplitObject(rna_obj, split.by = "COMBAT_ID_Time")
rna_list <- lapply(X=rna_list, FUN = function(x){
    x <- NormalizeData(x,verbose=F)
    x <- FindVariableFeatures(x, nfeatures = 1000, verbose=F)
})

features <- SelectIntegrationFeatures(object.list = rna_list, features=1000)

table(unlist(lapply(rna_list, ncol)<80))
FALSE 
   89 
#rna_list_all <- rna_list
#rna_list <- rna_list_all[unlist(lapply(rna_list_all, ncol)>80)]
rna_list <- lapply(X = rna_list, FUN = function(x) {
	 x <- ScaleData(x, features = features, verbose = FALSE, vars.to.regress = ("nCount_RNA"))
	 x <- RunPCA(x, features = features, verbose = FALSE, npcs=30)
})

HC_samples <- colData(sce) %>% as_tibble () %>% filter(condition=="HV") %>% distinct(COMBAT_ID_Time)
rna_anchors <- FindIntegrationAnchors(object.list = rna_list, reference=which(names(rna_list) %in% HC_samples$COMBAT_ID_Time), 
   	    dims = 1:30, scale = F, reduction = "rpca")
saveRDS(rna_anchors, file="/t1-data/project/taylorlab/erepapi/Fellowship/COVID19_Integration/RData/RNAobjintegrated_toy_Bcells_nosepsis_rnaanchors.RDS")

#plan("sequential", gc=T,workers=4)
all_features <- lapply(rna_list, row.names) %>% Reduce(intersect, .) 

rna.integrated <- IntegrateData(anchorset = rna_anchors, dims = 1:50, k.weight = 80, normalization.method = "LogNormalize") #k.weight needs to be 80 because there are samples with less than 100 cells (which is the default)
rna.integrated <- ScaleData(rna.integrated, verbose = FALSE)
rna.integrated <- RunPCA(rna.integrated, verbose = FALSE)
rna.integrated <- RunUMAP(rna.integrated, dims = 1:15)#30

DefaultAssay(rna.integrated) <- "RNA"
rna.integrated <- ScaleData(rna.integrated, features = row.names(rna.integrated))
rna.integrated <- RunPCA(rna.integrated, features = row.names(rna.integrated),reduction.name = "rna.pca" )
rna.integrated <- RunUMAP(rna.integrated, dims = 1:15, reduction = "rna.pca",assay = "RNA", reduction.name = "rna.umap", reduction.key = "rnaUMAP_")

saveRDS(rna.integrated, file="/t1-data/project/taylorlab/erepapi/Fellowship/COVID19_Integration/RData/RNAobjintegrated_toy_Bcells_nosepsis.RDS")

# ADT
adt_obj <- readRDS("/t1-data/project/taylorlab/erepapi/Fellowship/COVID19_Integration/RData/ADTobj_toy_Bcells.RDS")
adt_obj <- adt_obj[,colnames(rna.integrated)]
adt_obj <- ScaleData(adt_obj, features = row.names(adt_obj))
adt_obj <- RunPCA(adt_obj, features = row.names(adt_obj))
adt_obj <- RunTSNE(object = adt_obj, dims.use = 1:30, do.fast = TRUE)
adt_obj <- RunUMAP(adt_obj, dims = 1:30)

saveRDS(adt_obj, file="/t1-data/project/taylorlab/erepapi/Fellowship/COVID19_Integration/RData/ADTobj_toy_Bcells_nosepsis.RDS")

date()
sessionInfo()

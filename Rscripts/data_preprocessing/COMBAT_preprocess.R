args <- commandArgs(trailingOnly = FALSE)
rna_file <- args[1]
adt_file <- args[2]
cytof_file <- args[3]
slice_rna <- args[4]
slice_cytof <- args[5]
outdir <- args[6]

# load necessary libraries
suppressPackageStartupMessages({
	library(SingleCellExperiment)
	library(dplyr)
	library(Seurat)
})

# load initial objects
rna_obj <- readRDS(rna_file)
rna_md <- rna_obj@meta.data
adt_obj <- readRDS(adt_file)
sce <- readRDS(cytof_file)

# extra necessary functions
labels_transform <- function(x) {
        x <- gsub(" cells", replacement="", x, fixed=T)
        x <- gsub(" T", replacement="", x, fixed=T)
        x <- gsub("Vd2 gd", replacement="GDT",  x)
        x <- gsub("Plasmablasts", replacement="PB", x)
        x <- gsub("[cp]DC[12]*", replacement="DC", x)
        return(x)
    }

### CyTOF
## downsampling, converting sce object into Seurat and processing

# subsample cytof by sample_ID and extract overlapping samples with scRNA
sce_coldata_subs_smaller <- as_tibble(colData(sce))%>% filter(major_cell_type!="Unclassified")%>% filter(condition!="Sepsis")  %>%  
                            group_by(patient_id)  %>% slice_sample(n=slice_cytof) 
sce_coldata_subs_smaller <- as.data.frame(filter(sce_coldata_subs_smaller, COMBAT_ID_Time %in% rna_md$COMBAT_ID_Time))
row.names(sce_coldata_subs_smaller) <- sce_coldata_subs_smaller$cellID
sce_subs <- sce[,colData(sce)$cellID %in% sce_coldata_subs_smaller$cellID]

# Creating the Seurat obj
cytofdata <- assay(sce_subs)
colnames(cytofdata) <- colData(sce_subs)$cellID 
cytof_obj <- CreateSeuratObject(counts = cytofdata, assay = "Cytof", meta.data=sce_coldata_subs_smaller, project = "CyTOF_data")
cytof_obj@meta.data$major_cell_type <- labels_transform(cytof_obj@meta.data$major_cell_type)
Idents(cytof_obj) <- cytof_obj@meta.data$major_cell_type
cytof_obj <- ScaleData(cytof_obj)
cytof_obj <- RunPCA(cytof_obj, features=row.names(cytof_obj))
cytof_obj <- RunUMAP(cytof_obj, reduction = "pca", dims = 1:30)
cytof_obj <- RunTSNE(cytof_obj, reduction = "pca", dims = 1:30, do.fast = TRUE, check_duplicates = FALSE) 

saveRDS(cytof_obj, file=paste0(outdir,"/COMBAT_CyTOFobj_", slice_cytof,".RDS"))

### RNA
## subsetting the integrated (batch corrected) Seurat object

rna_subs_smaller <- as_tibble(rna_md) %>% group_by(scRNASeq_sample_ID)  %>% slice_sample(n=slice_rna) 
rna.subsampled <- subset(rna_obj, cell=rna_subs_smaller$barcode_id)
saveRDS(rna.subsampled, file=paste0(outdir,"/COMBAT_RNAobj_", slice_rna,".RDS"))

### ADT
## subsetting the Seurat object 

adt.subsampled <- subset(adt_obj, cell=rna_subs_smaller$barcode_id)
saveRDS(adt.subsampled, file=paste0(outdir,"/COMBAT_ADTobj_", slice_rna,".RDS"))

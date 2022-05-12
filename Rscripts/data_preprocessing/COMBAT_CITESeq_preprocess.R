### load necessary libraries
library(Seurat)
library(tidyverse)
library(dplyr)
library(SingleCellExperiment)

### load data
sce <- readRDS("cytof_sce.RDS")
adt <- readRDS("adt_mtx.RDS")
rna_mtx <- readRDS("rna_mtx.RDS")
cell_table <- readRDS("citeseq_md.RDS")
gene_protein_table <- read.table("../input_tables/gene_protein_table_COMBAT.txt", sep="\t", header=T)

### renaming the genes into their protein counterparts (from the gene-ADT-CyTOF table)
gene_protein_table <- filter(gene_protein_table, Gene %in% row.names(rna_mtx))
gene_protein_table$CyTOF <- gsub("-",replacement="_",gene_protein_table$CyTOF)
m <- match(gene_protein_table$Gene,row.names(rna_mtx))
row.names(rna_mtx)[m] <- gene_protein_table$CyTOF

### Filtering CITE-Seq datasets

# filtering samples with low numbers of cells
cell_table_subs <- filter(cell_table, !(scRNASeq_sample_ID %in% c("G05092-Ja005E-PBCa","S00030-Ja003E-PBCa"))) 
cell_table_subs <- filter(cell_table_subs, COMBAT_ID_Time %in% colData(sce)$COMBAT_ID_Time) 

# filtering low quality cells and cells with uncertain annotation
cell_table_subs <- filter(cell_table_subs, 
                        nfeatures_ADT>quantile(nfeatures_ADT, 0.001) & 
			total_UMI_ADT>quantile(total_UMI_ADT, 0.001) & 
			ngenes>quantile(ngenes, 0.001) & 
			total_UMI>quantile(total_UMI, 0.001)) %>%
                    filter(trusted_phenotype=="yes") %>%
                    filter(!grepl("|",major_cell_type, fixed=T) & major_cell_type!="NA") 

# exctracting 2000 cells per sample from the metadata table
cell_table_subs_smaller <- as.data.frame(cell_table_subs %>% 
				filter(!major_cell_type %in% c("RET","Mast")) %>% 
				group_by(scRNASeq_sample_ID) %>% 
				slice_sample(n=2000, replace=F))
cell_table_subs_smaller <- as.data.frame(filter(cell_table_subs_smaller, COMBAT_ID_Time %in% as_tibble(colData(sce))$COMBAT_ID_Time))
row.names(cell_table_subs_smaller) <- cell_table_subs_smaller$barcode_id

### ADT
## subsetting the Seurat object 

adt_subs <- adt[,row.names(cell_table_subs_smaller)]
adt_obj <- CreateSeuratObject(adt_subs, assay = "ADT", meta.data=cell_table_subs_smaller, project = "10X_ADT")
# the data is already normalised so I just scale:
adt_obj <- ScaleData(adt_obj, features = row.names(adt_obj))
adt_obj <- RunPCA(adt_obj, features = row.names(adt_obj))
adt_obj <- RunTSNE(object = adt_obj, dims.use = 1:30, do.fast = TRUE)
adt_obj <- RunUMAP(adt_obj, dims = 1:30)

adt.subsampled <- subset(adt_obj, cell=rna_subs_smaller$barcode_id)
saveRDS(adt_obj, file="COMBAT_ADTobj_182000.RDS")

### RNA

# subsetting the Seurat object 
rna_mtx_subs_smaller <- rna_mtx[,row.names(cell_table_subs_smaller)]
rna_obj <- CreateSeuratObject(rna_mtx_subs_smaller, assay = "RNA", meta.data=cell_table_subs_smaller, project = "10X_RNA")

# batch correction
Idents(rna_obj) <- "COMBAT_ID_Time"
rna_list <- SplitObject(rna_obj, split.by = "COMBAT_ID_Time")
rna_list <- lapply(X=rna_list, FUN = function(x){
    x <- NormalizeData(x,verbose=F)
    x <- FindVariableFeatures(x, nfeatures = 3000, verbose=F)
})

features <- SelectIntegrationFeatures(object.list = rna_list)
rna_list <- lapply(X = rna_list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE, vars.to.regress = ("nCount_RNA"))
    x <- RunPCA(x, features = features, verbose = FALSE)
})

# use the HC as reference
HC_samples <- colData(sce) %>% as.tibble () %>% filter(condition=="HV") %>% distinct(COMBAT_ID_Time)
rna_anchors <- FindIntegrationAnchors(object.list = rna_list, reference=which(names(rna_list) %in% HC_samples$COMBAT_ID_Time), # using only HC as reference, otherwise it doesn't run: Cholmod error 'problem too large'
    dims = 1:50, scale = F, reduction = "rpca")

all_features <- lapply(rna_list, row.names) %>% Reduce(intersect, .) 

# integrate
rna.integrated <- IntegrateData(anchorset = rna_anchors, dims = 1:50, features.to.integrate = features, normalization.method = "LogNormalize")
rna.integrated <- ScaleData(rna.integrated, verbose = FALSE)
rna.integrated <- RunPCA(rna.integrated, verbose = FALSE)
rna.integrated <- RunUMAP(rna.integrated, dims = 1:30)

# also running scaling, PCA and UMAP on the uncorrected object
DefaultAssay(rna.integrated) <- "RNA"
rna.integrated <- ScaleData(rna.integrated, features = row.names(rna.integrated))
rna.integrated <- RunPCA(rna.integrated, features = row.names(rna.integrated),reduction.name = "rna.pca" )
rna.integrated <- RunUMAP(rna.integrated, dims = 1:30, reduction = "rna.pca",assay = "RNA", reduction.name = "rna.umap", reduction.key = "rnaUMAP_")

saveRDS(rna.integrated, file="COMBAT_RNAobjintegrated_182000.RDS")


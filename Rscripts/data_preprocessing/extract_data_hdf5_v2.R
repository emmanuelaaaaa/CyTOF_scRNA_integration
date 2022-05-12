
# script to extract only the Bcells from the downloaded HDF5 covid sc PBMC and WB  datatset HDF5 files from https://atlas.fredhutch.org/fredhutch/covid/

library(SingleCellExperiment)
library(Seurat)
library(SeuratDisk)
library(plyr)
library(tidyverse)
#library(randomcoloR)
library(DT)
library(data.table)
#setwd("/t1-data/project/sims-lab/dagarwal/Integration_paper/data/covid_sc_data")
setDTthreads(12)
data_dir <- "./hdf5_files"
all_meta <- fread("./all_meta.tsv")

# Deal with Whole Blood
pbmc_datasets <- unique(all_meta[sample_type == "Peripheral Blood Mononuclear Cells"]$dataset)
whole_blood_datasets <- unique(all_meta[sample_type == "Whole Blood"]$dataset)

# Takes vector of gene aliases and maps to HGNC gene symbol
hgncAlias2Symbol <- fread("hgncAlias2Symbol.tsv")
mapGeneSymbol <- function(gene_aliases) {
    dt <- data.table(ALIAS = gene_aliases)
    dt <- merge(dt, hgncAlias2Symbol, all.x = TRUE, all.y = FALSE, sort = FALSE)
    dt[is.na(SYMBOL), SYMBOL := ""]
    return(dt$SYMBOL)
}

extract_data <- function(seu, dataset_name, assay = "RNA") {
    
    # Map gene alias to symbol via HUGO db (as of 12/23/2020)
    # Advice is to NOT rename features: https://github.com/satijalab/seurat/issues/1049
    # We do it only for datasets to be merged.
    #seu[[assay]]@counts@Dimnames[[1]] <- mapGeneSymbol(seu[[assay]]@counts@Dimnames[[1]])
    #seu[[assay]]@data@Dimnames[[1]] <- mapGeneSymbol(seu[[assay]]@data@Dimnames[[1]])
    #rownames(seu[[assay]]@scale.data) <- mapGeneSymbol(rownames(seu[[assay]]@scale.data))
    #seu[[assay]]@var.features <- mapGeneSymbol(seu[[assay]]@var.features)
    
    # Subset to only genes which successfully mapped to aliases
    # seu[["RNA"]] <- subset(seu[["RNA"]], features = shared_genes)
    
    # Merge appropriate meta data
    merge_field <- unique(all_meta[dataset == dataset_name]$sample_source_field)
    seu_meta <- seu@meta.data[, c(merge_field,
                                  "predicted.celltype.l1",
                                  "predicted.celltype.l2",
                                  "predicted.celltype.l3",
                                  "predicted.celltype.l1.score",
                                  "predicted.celltype.l2.score",
                                  "predicted.celltype.l3.score",
                                  "mapping.score")]
    seu_meta$id <- rownames(seu_meta)
    seu_meta$dataset <- dataset_name
    seu_meta <- merge(seu_meta,
                      all_meta,
                      by.x = c(merge_field, "dataset"),
                      by.y = c("sample", "dataset"),
                      all.x = TRUE,
                      all.y = FALSE)
    seu_meta$sample <- seu_meta[[merge_field]]
    rownames(seu_meta) <- seu_meta$id
    seu_meta <- seu_meta[, c("dataset",
                             "tissue",
                             "sample_type",
                             "sample_type_note",
                             "patient",
                             "sample",
                             "race_reported",
                             "sex",
                             "age",
                             "disease_status",
                             "disease_severity",
                             "days_since_symptom_onset",
                             "predicted.celltype.l1",
                             "predicted.celltype.l2",
                             "predicted.celltype.l3",
                             "predicted.celltype.l1.score",
                             "predicted.celltype.l2.score",
                             "predicted.celltype.l3.score",
                             "mapping.score")]
    # Set it in the correct order
    seu_meta <- seu_meta[rownames(seu@meta.data),]
    seu@meta.data <- seu_meta
    DefaultAssay(seu) <- "RNA"
    seu[["ref.umap"]]@assay.used <- "SCT"
    
    # Subset only to Bcells
    Idents(seu) <- "predicted.celltype.l1"
    seu <- subset(x =seu,idents ="B")
    
    DietSeurat(seu,
               assays = c("RNA", "predicted_ADT","ADT"),
               dimreducs = c("ref.umap", "rna.umap","pca","ref.spca"))
}

#Just PBMCs
seuList_pbmc <- lapply(pbmc_datasets, function(dataset_name) {
    message("-----", dataset_name, "-----")
    #loadDir <- h5seuratDir
    
    seu <- LoadH5Seurat(file.path(data_dir, paste0(dataset_name, "_processed.HDF5")),
                        assays = c("RNA", "predicted_ADT"),
                        reductions = c("ref.umap", "rna.umap","pca","ref.spca"),
                        meta.data = TRUE,
                        verbose = FALSE)
    # Subset to only PBMCs
    pbmc_samples <- all_meta[dataset == dataset_name &
                                 sample_type == "Peripheral Blood Mononuclear Cells"]$sample
    sample_field <- unique(all_meta[dataset == dataset_name]$sample_source_field)
    seu <- seu[,seu[[sample_field]][,] %in% pbmc_samples]
    
    
    return(extract_data(seu, dataset_name, "RNA"))
})

names(seuList_pbmc) <- pbmc_datasets
saveRDS(seuList_pbmc,file=paste0(data_dir,"/","covid_sc_pbmc_bcells_only_v2.rds"))

# Just whole blood
seuList_wb <- lapply(whole_blood_datasets, function(dataset_name) {
    message("----- ", dataset_name, " whole blood -----")
    #path <- file.path(data_dir, list.files(wholeBloodDir, pattern = dataset_name))
    # if (length(path) != 1) stop("Found ", length(path), " whole blood files for ", dataset_name)
    seu <- LoadH5Seurat(file.path(data_dir, paste0(dataset_name, "_processed.HDF5")),
                        assays = c("ADT"),
                        reductions = c("ref.umap","rna.umap","pca","ref.spca"),
                        meta.data = TRUE,
                        verbose = FALSE)
#    seu <- seu[, !seu$monaco_main %in% c("Neutrophils", "Basophils")]
    
    # For 2 datasets where it is missing, add derived "sample" column
    # (normally would happen when writing metadata)
  #  if (dataset_name == "silvin_2020") {
     #   seu <- AddMetaData(seu,
     #                      paste0(trimws(seu$Characteristics.individual.),
      #                            "_",
       #                           seu$Factor.Value.sampling.time.point.),
        #                   "sample")
#    } else if (dataset_name == "bost_2021") {
 #       seu <- AddMetaData(seu,
  #                         paste0(trimws(seu$subject_id),
   #                               "_",
    #                              seu$tissue),
     #                      "sample")
    #} else if (dataset_name == "schulte-schrepping_2020") {
     #   seu <- AddMetaData(seu,
      #                     paste0(seu$sampleID, "_", seu$cells),
       #                    "sampleid_unique")
#    }
    
    return(extract_data(seu, dataset_name, "RNA"))
})

names(seuList_wb) <- whole_blood_datasets
saveRDS(seuList_wb,file=paste0(data_dir,"/","covid_sc_wb_bcells_only_v2.rds"))
#seuList <- c(seuList_pbmc, seuList_wb)

#saveRDS(seuList, file=paste0(data_dir,"/","covid_sc_pbmc_wb_bcells_only.rds"))


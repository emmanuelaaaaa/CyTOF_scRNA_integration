args <- commandArgs(trailingOnly = TRUE)
rna_file <- args[1]
adt_file <- args[2]
plotname <- args[3]

# load necessary libraries
suppressPackageStartupMessages({
	library(ggplot2)
	library(Seurat)
	library(cowplot)
	library(scuttle)
})

# load initial objects
rna <- readRDS(rna_file)
adt <- readRDS(adt_file)
adt <- adt[,colnames(rna)]

rna_sc <- as.SingleCellExperiment(rna, assay="Cytof")
adt_sc <- as.SingleCellExperiment(adt, assay="ADT")

rna_pseudobulk <- aggregateAcrossCells(rna_sc, id=colData(rna_sc)[,c("pseudobulk", "scRNASeq_sample_ID")], statistics="mean", use.assay.type="logcounts")
adt_pseudobulk <- aggregateAcrossCells(adt_sc, id=colData(adt_sc)[,c("pseudobulk", "scRNASeq_sample_ID")], statistics="mean", use.assay.type="logcounts")

if(all(colData(rna_pseudobulk)$pseudobulk==colData(adt_pseudobulk)$pseudobulk & colData(rna_pseudobulk)$scRNASeq_sample_ID==colData(adt_pseudobulk)$scRNASeq_sample_ID)) {
  pseudo_imputedvsreal <- data.frame(t(assay(rna_pseudobulk, "logcounts")[c("CD45RO","CD45RA"),]),t(assay(adt_pseudobulk, "logcounts")[c("CD45RO","CD45RA"),]), 
                                     colData(adt_pseudobulk)$pseudobulk, colData(adt_pseudobulk)$scRNASeq_sample_ID)
  names(pseudo_imputedvsreal) <- c("CyTOF_imputed_CD45RO","CyTOF_imputed_CD45RA","ADT_CD45RO","ADT_CD45RA","pseudobulk","scRNASeq_sample_ID")
} else {
  print("Records not matching between adt and rna pseudobulks!")
}
cor_RO <- round(cor(pseudo_imputedvsreal$CyTOF_imputed_CD45RO, pseudo_imputedvsreal$ADT_CD45RO), digits=3)
cor_RA <- round(cor(pseudo_imputedvsreal$CyTOF_imputed_CD45RA, pseudo_imputedvsreal$ADT_CD45RA), digits=3)

p1 <- ggplot(pseudo_imputedvsreal, aes(CyTOF_imputed_CD45RO, ADT_CD45RO, col=pseudobulk)) + geom_point(size=1) +
            annotate("text", x=2.25, y=0, label=paste0("R = ",cor_RO), color = "black")
p2 <- ggplot(pseudo_imputedvsreal, aes(CyTOF_imputed_CD45RA, ADT_CD45RA, col=pseudobulk)) + geom_point(size=1) +
              annotate("text", x=2, y=0, label=paste0("R = ",cor_RA), color = "black")

legend <- get_legend(p1 + theme(legend.position="bottom"))
p <- plot_grid(plot_grid(p1 + theme(legend.position="none"), p2 + theme(legend.position="none")), 
          legend, ncol = 1, rel_heights = c(1.5,1))
ggsave(p, filename=paste0(plotname,'.pdf'), width = 12, height = 7)
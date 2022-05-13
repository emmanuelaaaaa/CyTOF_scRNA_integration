args <- commandArgs(trailingOnly = TRUE)
rna_file <- args[1]
adt_file <- args[2]
cytof_file <- args[3]
plotname <- args[4]

# load necessary libraries
suppressPackageStartupMessages({
	library(ggplot2)
	library(Seurat)
	library(cowplot)
	library(ggrepel)
})

# load initial objects
rna <- readRDS(rna_file)
adt <- readRDS(adt_file)
adt <- adt[,colnames(rna)]
cytof <- readRDS(cytof_file)
common_features <- intersect(row.names(rna), row.names(cytof))
common_features_adt <- intersect(row.names(cytof), row.names(adt))

# run the correlations
cor_ADT_impCyt <- c()
for (i in common_features_adt) {
    cor_ADT_impCyt <- c(cor_ADT_impCyt,cor(as.vector(Assays(adt, slot="ADT")[i,]),as.vector(Assays(rna, slot="Cytof")[i,])))
}
names(cor_ADT_impCyt) <- common_features_adt

cor_ADT_RNA <- c()
for (i in intersect(common_features, common_features_adt)) {
    cor_ADT_RNA <- c(cor_ADT_RNA,cor(as.vector(Assays(adt, slot="ADT")[i,]),as.vector(Assays(rna, slot="RNA")[i,])))
}
names(cor_ADT_RNA) <- intersect(common_features, common_features_adt)

cor_impCyt_RNA <- c()
for (i in common_features) {
    cor_impCyt_RNA <- c(cor_impCyt_RNA,cor(as.vector(Assays(rna, slot="RNA")[i,]), as.vector(Assays(rna, slot="Cytof")[i,])))
}
names(cor_impCyt_RNA) <- common_features

# create the tables of the correlations
cor_RNA <- merge(cor_ADT_RNA,cor_impCyt_RNA, by.x=0, by.y=0)
names(cor_RNA) <- c("feature", "Corr_RNA_ADT", "Corr_ImputedCyTOF_RNA")
cor_ADT <- merge(cor_ADT_RNA,cor_ADT_impCyt, by.x=0, by.y=0)
names(cor_ADT) <- c("feature", "Corr_RNA_ADT", "Corr_ImputedCyTOF_ADT")

# plot
p1 <- ggplot(cor_RNA,aes(Corr_RNA_ADT, Corr_ImputedCyTOF_RNA, label=feature)) + xlim(-0.3, 1) + ylim(-0.15,1) +
                xlab(paste0("correlations of RNA and ADT values")) +
                ylab(paste0("correlations of RNA \nand Seurat imputed CyTOF values")) + geom_point(colour="darkgreen") +
                geom_text_repel(size=5) + geom_abline(slope=1, linetype="dotted") + theme_minimal()+ theme(text = element_text(size=15))

p2 <- ggplot(cor_ADT,aes(Corr_RNA_ADT, Corr_ImputedCyTOF_ADT, label=feature)) + xlim(-0.3, 1) + ylim(-0.15,1) +
                xlab(paste0("correlations of RNA and ADT values")) +
                ylab(paste0("correlations of Seurat imputed CyTOF markers \n and ADT values")) + geom_point(colour="darkgreen") +
                geom_text_repel(size=5) + geom_abline(slope=1, linetype="dotted") + theme_minimal()+ theme(text = element_text(size=15))
p <- plot_grid(p1, p2)
ggsave(p, filename=paste0(plotname,'.pdf'), width = 12, height = 6)

args <- commandArgs(trailingOnly = TRUE)
rna_file <- args[1]
plotname <- args[2]

# load necessary libraries
suppressPackageStartupMessages({
	library(Seurat)
	library(flipPlots)
	library(htmlwidgets)
})

# load initial objects
rna <- readRDS(rna_file)

library(flipPlots)
edges <- as.data.frame(table(rna@meta.data$major_cell_type, rna@meta.data$predicted.id), stringAsFactors=F)
colnames(edges)[1] <- "Original_Celltype"
colnames(edges)[2] <- "Predicted_Celltype"
colnames(edges)[3] <- "Value"

edges$Original_Celltype <- as.factor(edges$Original_Celltype)
edges$Predicted_Celltype <- as.factor(edges$Predicted_Celltype)
 
p1 <- SankeyDiagram(edges[,-3],
              link.color = "First variable",
              weights = edges$Value, max.categories = 23, hovertext.show.percentages = T, label.show.counts = T, font.size = 10, variables.share.values =T, label.show.percentages = T) 
saveWidget(p1 , file=paste0(plotname,'.html'))

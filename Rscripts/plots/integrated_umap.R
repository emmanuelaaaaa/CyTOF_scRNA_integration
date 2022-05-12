args <- commandArgs(trailingOnly = TRUE)
rna_file <- args[1]
cytof_file <- args[2]
plotname <- args[3]

# load necessary libraries
suppressPackageStartupMessages({
	library(ggplot2)
	library(Seurat)
	library(cowplot)
})

# colours used
c38 <- c("dodgerblue2", "#E31A1C", # red
                "green4", 
                "#6A3D9A", # purple
                "#FF7F00", # orange
                "black", "gold1",
                "skyblue2", "#FB9A99", # lt pink
                "palegreen2",
                "#CAB2D6", # lt purple
                "#FDBF6F", # lt orange
                "darkgrey", "khaki2",
                "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
                "darkturquoise", "green1", "yellow4", "yellow3",
                "darkorange4", "grey","darkorchid4","brown",
                "azure3","burlywood1","aquamarine","cadetblue","coral2",
                "darkolivegreen1","darkgoldenrod2","brown1","darkgreen","blueviolet","darkred"
              )

# load initial objects
rna <- readRDS(rna_file)
cytof <- readRDS(cytof_file)

# Merging the two datasets
seurat_integrated <- merge(x = rna, y = cytof)
DefaultAssay(seurat_integrated) <- "Cytof"
# Centering the data matrix
common_features <- intersect(row.names(rna), row.names(cytof))
seurat_integrated <- ScaleData(seurat_integrated, features = common_features, do.scale = FALSE)
seurat_integrated <- RunUMAP(seurat_integrated, features = common_features, reduction.name="umap.cytof", verbose = F)

# plotting the embeddings
p1 <- DimPlot(seurat_integrated, reduction = "umap.cytof", pt.size=0.1, group.by="orig.ident") + 
	theme(legend.title = element_text( size=5), legend.text=element_text(size=5)) + facet_wrap(~orig.ident) + theme(legend.position = "none")
p2 <- DimPlot(seurat_integrated, reduction = "umap.cytof", pt.size=0.1, group.by="major_cell_type") + scale_color_manual(values=c38) + 
	theme(legend.title = element_text( size=5), legend.text=element_text(size=5))
p <- plot_grid(p1, p2, rel_widths=c(1,0.6))
ggsave(p, filename=paste0(plotname,'.png'), width = 20, height = 10, dpi = 300, type="cairo-png")


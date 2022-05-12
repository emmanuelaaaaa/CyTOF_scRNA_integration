args <- commandArgs(trailingOnly = TRUE)
rna_file <- args[1]
plotname <- args[2]

# load necessary libraries
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(Seurat)
	library(cowplot)
})

# load initial objects
rna <- readRDS(rna_file)

freq_comp <- rna@meta.data[,c("COMBAT_ID_Time","major_cell_type","predicted.id","prediction.score.max")]
freqannot1 <- freq_comp %>%
            group_by(COMBAT_ID_Time, major_cell_type) %>%
            summarise(n = n()) %>%
            mutate(freq = n / sum(n)) %>%
            mutate(var1 = paste(COMBAT_ID_Time,major_cell_type))

freqannot2 <- freq_comp %>%
            group_by(COMBAT_ID_Time, predicted.id) %>%
            summarise(n = n()) %>%
            mutate(freq = n / sum(n))  %>%
            mutate(var1 = paste(COMBAT_ID_Time,predicted.id))

freqannot <- merge(freqannot1, freqannot2, by.x="var1", by.y="var1")

p1 <- ggplot(freqannot, aes(freq.x, freq.y)) + geom_point() + 
    labs(y="Frequencies of predicted clusters", x="Frequencies of RNA clusters")
p <- p1 + facet_wrap(~major_cell_type, scales = "free") + geom_abline(slope=1,linetype = 2, col="grey") 
ggsave(p, filename=paste0(plotname,'.pdf'), width = 12, height = 8)
library(Seurat)
library(dplyr)
ags<- commandArgs(trailingOnly = TRUE)
if(length(ags)<3)
{
  cat("Usage: TopFeatures RDSin PNGout #Features\n")
  q(save = "no")
}
rds <- readRDS(file = ags[1])
features=as.integer(ags[3])
png(filename = ags[2],width = 1000,height = 200+10*features)
top <- rds.markers %>% group_by(cluster) %>% top_n(features, avg_logFC)
# setting slim.col.label to TRUE will print just the cluster IDS instead of
# every cell name
DoHeatmap(object = rds, genes.use = top$gene, slim.col.label = TRUE, remove.key = TRUE)
dev.off()
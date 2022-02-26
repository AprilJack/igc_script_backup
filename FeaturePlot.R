library(Seurat)
ags<- commandArgs(trailingOnly = TRUE)
if(length(ags)<3)
{
  cat("Usage: FeaturePlot RDSin PNGout Feature1 [Feature2] ... [FeatureN]\n")
  q(save = "no")
}
rds <- readRDS(file = ags[1])
png(filename = ags[2],width = min(3000,200*(sqrt(length(ags)-2)+1)),min(3000,height = 200*(sqrt(length(ags)-2)+1)))
FeaturePlot(object = rds, features.plot = ags[3:length(ags)],cols.use = c("lightgrey", "blue"))
dev.off()
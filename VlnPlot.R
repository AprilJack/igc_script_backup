library(Seurat)
ags<- commandArgs(trailingOnly = TRUE)
if(length(ags)<3)
{
  cat("Usage: VlnPlot RDSin PNGout Feature1 [Feature2] ... [FeatureN]\n")
  q(save = "no")
}
rds <- readRDS(file = ags[1])
png(filename = ags[2],width = min(3000,100+200*(length(ags)-2)),height = 800)
VlnPlot(object = rds, features.plot = ags[3:length(ags)], x.lab.rot = TRUE)
dev.off()
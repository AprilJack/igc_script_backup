library(Seurat)
ags<- commandArgs(trailingOnly = TRUE)
if(length(ags)<3)
{
  cat("Usage: DimHeatmap RDSin PNGout #DimsToPlot\n")
  q(save = "no")
}
rds <- readRDS(file = ags[1])
dims=as.integer(ags[3])
png(filename = ags[2],width = 200+200*sqrt(dims),height = 200+200*sqrt(dims))
DimHeatmap(rds, dim.use = 1:dims, cells.use = 500, do.balanced = TRUE)
dev.off()
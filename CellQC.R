library(Seurat)
ags<- commandArgs(trailingOnly = TRUE)
if(length(ags)<3)
{
  cat("Usage: CellQC RDSin RDSout PNGout\n")
  cat("Assumes you have genes that start with MT-\n")
  q(save = "no")
}
rds <- readRDS(file = ags[1])
mito.genes <- grep(pattern = "^MT-", x = rownames(x = rds@data), value = TRUE)
percent.mito <- Matrix::colSums(rds@raw.data[mito.genes, ])/Matrix::colSums(rds@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
rds <- AddMetaData(object = rds, metadata = percent.mito, col.name = "percentMito")

saveRDS(rds,file=ags[2])
png(filename = ags[3],width = 800,height = 600)
VlnPlot(rds, features = c("nUMI", "nGene", "percentMito"),nCol = 3)
dev.off()
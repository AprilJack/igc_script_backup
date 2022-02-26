

ags<- commandArgs(trailingOnly = TRUE)
if (length(ags)==0){
  cat("Usage: rlog raw.txt\n",stderr())
  ags=c("raw.txt")
}

data <- as.matrix(read.delim(file = ags[1],header = TRUE, sep = "\t"))
if (!require("DESeq2")) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("DESeq2")
  library(DESeq2)
}
if(dim(data)[2]>8)
{
  if(!is.numeric(data[1,8]))
  {
    write("Assuming HOMER formated raw count table\n",stderr())
    data2 = data[,9:dim(data)[2]]
    for( i in 1:dim(data)[1]) 
    {
      data[i,8]= strsplit(data[i,8],"[|]")[[1]][1]
    }
    data2= matrix(as.integer(data2),nrow=dim(data2)[1],ncol=dim(data2)[2])
    rownames(data2) <- data[,8]
    colnames(data2) <- colnames(data[,9:dim(data)[2]])
    colData=data.frame(factor(rep("X",dim(data2)[2])))
    colnames(colData)<-"Condition"
    dds <- DESeqDataSetFromMatrix(data2, colData, ~ 1)
    dds <- estimateSizeFactors(dds)
    dds <- estimateDispersions(dds)
    write("Running rlog normalization. This could take a moment",stderr())
    norm=rlog(dds,blind = TRUE)
    norm_matrix=assay(norm)
    rownames(norm_matrix)<-rownames(data2)
    colnames(norm_matrix)<-colnames(data2)
    write.table(cbind(data[,1:8],norm_matrix),quote = FALSE,sep="\t",row.names = FALSE,col.names = NA, stdout())
  }  
} else {
  write("Assuming just one column for ids and one row for sample names\n",stderr())
  data2 = data[,2:dim(data)[2]]
  data2= matrix(as.integer(data2),nrow=dim(data2)[1],ncol=dim(data2)[2])
  rownames(data2) <- data[,1]
  colnames(data2) <- colnames(data[,2:dim(data)[2]])
  colData=data.frame(factor(rep("X",dim(data2)[2])))
  colnames(colData)<-"Condition"
  dds <- DESeqDataSetFromMatrix(data2, colData, ~ 1)
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  write("Running rlog normalization. This could take a moment",stderr())
  norm=rlog(dds,blind = TRUE)
  norm_matrix=assay(norm)
  rownames(norm_matrix)<-rownames(data2)
  colnames(norm_matrix)<-colnames(data2)
  write.table(norm_matrix,quote = FALSE,sep="\t",row.names = TRUE,col.names = NA, stdout())
}

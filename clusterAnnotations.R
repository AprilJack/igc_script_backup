if (!require("gplots")) {
  install.packages("gplots", repos="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", repos ="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(RColorBrewer)
}

ags<- commandArgs(trailingOnly = TRUE)
if(length(ags)==0){
  cat("clusterAnnotations will cluster the normalized expression values of TranscriptAnn\noutput and produce a data-table of the log-transformed values as well as PNG image of the clustering\n\nUsage:clusterAnnotations annFile1.txt [[annFile2.txt]...]\n")
}else
{
  for( f in ags)
  {
    if(file.info(f)$size > 0)
    {
      data <- as.matrix(read.delim(file = f,header = TRUE, sep = "\t",))
      data[is.na(data)] <- 0
      data2 = data[,10:(dim(data)[2]-1)]
      data2 = data2[apply(data[,10:(dim(data)[2]-1)],1,sd)>0,apply(data[,10:(dim(data)[2]-1)],2,sd)>0]
      if(min(as.numeric(data2)) >=0 ) {
        data2 = log(matrix(as.numeric(data2)+0.000001,nrow=dim(data2)[1],ncol=dim(data2)[2]))
      }else {
        data2 = matrix(as.numeric(data2),nrow=dim(data2)[1],ncol=dim(data2)[2])
      }
      rownames(data2) <- data[apply(data[,10:(dim(data)[2]-1)],1,sd)>0,1]
      trimmed = data[,10:(dim(data)[2]-1)]
      colnames(data2) <- colnames(trimmed[,apply(trimmed,2,sd)>0])
      write.table(data2,file = sprintf("%s_clustering.txt",f),sep = "\t", eol = "\r\n",col.names=NA)
      my_palette <- colorRampPalette(c("white","white","white","lightpink", "red","black"))(n=299)
      png(sprintf("%s.png",f), width=max(300*5,20*dim(data2)[2]), max(300*5,height=20*dim(data2)[1]),res=300,pointsize=8)
      heatmap.2(data2, main = "Ln(NormExp) Clustering", trace="none", margins=c(10,10),col=my_palette,cexRow = 0.3,cexCol=0.3,distfun=function(c) as.dist(1-cor(t(c))), ,key.title = "",hclustfun = function(x) hclust(x,method="ward.D2"))
      dev.off()
    }else cat(sprintf("The file %s was empty\n",f))
  }
  cat("Finished\n")
}

if (!require("gplots")) {
  install.packages("gplots", repos="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", repos ="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(RColorBrewer)
}

if (!require("DESeq2")) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("DESeq2")
}
ags<- commandArgs(trailingOnly = TRUE)
if (length(ags)==0){
  cat("Usage: normalizeRaw raw_file [clusters to output as separate text files] [Simple format?TRUE/FALSE]\n")
  ags=c("raw.txt",3,FALSE)
}
if(length(ags)==1) {
  ags=c(ags[1],"3",FALSE)
}
if(length(ags)==2) {
      ags=c(ags[1],ags[2],FALSE)
}
data <- as.matrix(read.delim(file = ags[1],header = TRUE, sep = "\t"))
if(ags[3]==TRUE)
{
  data2 = data[,2:dim(data)[2]]
  data2b = data[,2:dim(data)[2]]
  data2= matrix(as.numeric(data2),nrow=dim(data2)[1],ncol=dim(data2)[2])
  data2b= matrix(as.numeric(data2b),nrow=dim(data2b)[1],ncol=dim(data2b)[2])  
  data3=
  rownames(data2) <- data[,1]
  colnames(data2) <- colnames(data[,2:dim(data)[2]])
  rownames(data2b) <- data[,1]
  colnames(data2b) <- colnames(data[,2:dim(data)[2]])
  cat("Assuming just one column for ids and one row for sample names")
} else {
  data2 = data[,9:dim(data)[2]]
  data2b = data[,9:dim(data)[2]]
  for( i in 1:dim(data)[1]) 
  {
    data[i,8]= strsplit(data[i,8],"[|]")[[1]][1]
    if( i %% 1000==0)
      cat(".")
  }
  data2= matrix(as.numeric(data2),nrow=dim(data2)[1],ncol=dim(data2)[2])
  data2b= matrix(as.numeric(data2b),nrow=dim(data2b)[1],ncol=dim(data2b)[2])
  rownames(data2) <- data[,8]
  colnames(data2) <- colnames(data[,9:dim(data)[2]])
  rownames(data2b) <- data[,1]
  colnames(data2b) <- colnames(data[,9:dim(data)[2]])
  cat("Assuming HOMER formated fpkm table (otherwise set third arg to TRUE")
}
write.table(data2,file = sprintf("%s_fpkm_plus5.txt",ags[1]),sep = "\t", eol = "\r\n",col.names=NA)
if (length(ags) > 1)
{
  #user wants to get the top x clusters
  hc<-hclust(as.dist(1-cor(t(data3b))),method="ward.D2")
  ct = cutree(hc,k=ags[2])
  clusters=cbind(ct,data3b)
  cColors = rep(c("red","green","blue","yellow","cyan","gray","black","orange","purple","brown"),10)
  for ( i in 1:ags[2])
  {
    cluster=as.data.frame(clusters[clusters[,1]==i,])
    rownames(cluster)=rownames(clusters)[clusters[,1]==i]
    colnames(cluster)[1]="Cluster ID"
    cluster=rbind(colnames(cluster),cluster)
    rownames(cluster)[1]="Genes"
    write.table(cluster,file=sprintf("%s_cluster_%s.txt",ags[1],cColors[i]),sep = "\t",col.names = FALSE,quote=FALSE)
  }
  cat(sprintf("\nClassified the genes into %s clusters.\n",ags[2]))
  for(i in 1:ags[2])
  {
    ct[ct==i]<-cColors[i]
  }
}

my_palette <- colorRampPalette(c("green","black","red"))(n=299)
png(sprintf("%s_bigCluster.png",ags[1]), width=3000, height=3000,res=300,pointsize=8)
heatmap.2(data3, main = "Ln(FPKM+5) Clustering",trace="none", col=my_palette,margins = c(12,12),cexRow = 0.25,cexCol=0.5,distfun=function(c) as.dist(1-cor(t(c))),scale="row", hclustfun = function(x) hclust(x,method="ward.D2"),labRow=c(""),key.title = "")
dev.off()
png(sprintf("%s_cluster.png",ags[1]), width=1200, height=1200,res=300,pointsize=8)
heatmap.2(data3, main = "Ln(FPKM+5) Clustering", RowSideColors = ct,  margins = c(10,10) ,trace="none", col=my_palette,cexRow = 0.1,cexCol=0.5,distfun=function(c) as.dist(1-cor(t(c))),scale="row",labRow=c(""),key.title = "",hclustfun = function(x) hclust(x,method="ward.D2"))
dev.off()
cat("Finished\n")

if (!require("gplots")) {
  install.packages("gplots", repos="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", repos ="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(RColorBrewer)
}
ags<- commandArgs(trailingOnly = TRUE)
if (length(ags)==0){
  cat("Usage: fpkmKmeans fpkm_file [clusters to output as separate text files] [Simple format?TRUE/FALSE]\n")
  ags=c("fpkm.txt",3,FALSE)
}
if(length(ags)==1) {
  ags=c(ags[1],"3",FALSE)
}
if(length(ags)==2) {
      ags=c(ags[1],ags[2],FALSE)
}
cat(sprintf("Reading in %s\n",ags[1]))
data <- as.matrix(read.delim(file = ags[1],header = TRUE, sep = "\t"))
if(dim(data)[2]>8 & !is.numeric(as.numeric(data[1,8])))
{
  data2 = data[,2:dim(data)[2]]
  data2b = data[,2:dim(data)[2]]
  data2= matrix(as.numeric(data2),nrow=dim(data2)[1],ncol=dim(data2)[2])
  data2b= matrix(as.numeric(data2b),nrow=dim(data2b)[1],ncol=dim(data2b)[2])  
  rownames(data2) <- data[,1]
  colnames(data2) <- colnames(data[,2:dim(data)[2]])
  rownames(data2b) <- data[,1]
  colnames(data2b) <- colnames(data[,2:dim(data)[2]])
  data3= log(data2[apply(data2, 1, function(x) !all(x<16)&var(x)>0),]+1,base=2)
  dim(data3)
  data3b= log(data2b[apply(data2b, 1, function(x) !all(x<16)&var(x)>0),]+1,base=2)
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
  data3= log(data2[apply(data2, 1, function(x) !all(x<16)&var(x)>0),]+1,base=2)
  dim(data3)
  data3b= log(data2b[apply(data2b, 1, function(x) !all(x<16)&var(x)>0),]+1,base=2)
  cat("Assuming HOMER formated fpkm table (otherwise set third arg to TRUE")
}
write.table(data2,file = sprintf("%s_LogPlus1.txt",ags[1]),sep = "\t", eol = "\r\n",col.names=NA)
if (length(ags) > 1)
{
  #user wants to get the top x clusters
  hc<-kmeans(data3,centers=ags[2],trace=FALSE)
  clusters=cbind(hc$cluster,rownames(data3),data3b)
  cColors = rep(c("red","green","blue","yellow","cyan","gray","black","orange","purple","brown"),10)
  for ( i in 1:ags[2])
  {
    cluster=as.data.frame(clusters[clusters[,1]==i,])
    rownames(cluster)=rownames(clusters)[clusters[,1]==i]
    colnames(cluster)[1]="Cluster ID"
    colnames(cluster)[2]="Symbol"
    write.table(cluster,file=sprintf("%s_cluster_%s.txt",ags[1],cColors[i]),sep = "\t",col.names = NA,quote=FALSE,row.names = TRUE)
  }
  cat(sprintf("\nClassified the genes into %s clusters.\n",ags[2]))
  for(i in 1:ags[2])
  {
    hc$color[hc$cluster==i]<-cColors[i]
  }
}

my_palette <- colorRampPalette(c("green","black","red"))(n=299)
if(dim(data3)[2] > 30)
{
  my_palette <- colorRampPalette(c("green","green","black","red","red"))(n=299)
  png(sprintf("%s_bigCluster.png",ags[1]), width=max(50*dim(data3)[2],3000), height=3000,res=300,pointsize=8)
  heatmap.2(data3[order(hc$cluster),], main = sprintf("Log2(FPKM) %g",dim(data3)[1]),trace="none", RowSideColors = hc$colors,scale="row", col=my_palette,margins = c(12,12),cexRow = 0.25,cexCol=0.5,Rowv = FALSE,Colv = FALSE,labRow=c(""),key.title = "")
  dev.off()
} else {
  png(sprintf("%s_cluster.png",ags[1]), width=1500, height=1500,res=300,pointsize=8)
  if(dim(data3)[1] < 100)
  {
    heatmap.2(data3, main = sprintf("Log2(FPKM) %g",dim(data3)[1]), RowSideColors = ct,  margins = c(10,10) ,scale="row",trace="none", col=my_palette,cexRow = 0.5,cexCol=0.5,distfun=function(c) as.dist(1-cor(t(c))),key.title = "",hclustfun = function(x) hclust(x,method="ward.D2"))
  }else {
    heatmap.2(data3, main = sprintf("Log2(FPKM) %g",dim(data3)[1]), RowSideColors = ct,  margins = c(10,10) ,scale="row",trace="none", col=my_palette,cexRow = 0.1,cexCol=0.5,distfun=function(c) as.dist(1-cor(t(c))),key.title = "",labRow = c(""), hclustfun = function(x) hclust(x,method="ward.D2"))
  }
  dev.off()
}
cat("Finished\n")

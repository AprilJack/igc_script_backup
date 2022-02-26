library(gplots)
library(RColorBrewer)
ags<- commandArgs(trailingOnly = TRUE)
if (length(ags)==0){
  cat("Usage: clusterR homer_file [cluster samples?]\n")
  ags=c("tpm.txt",TRUE)
}
if(length(ags)==1) {
  ags=c(ags[1],TRUE)
}
if(length(ags)==2) {
      ags=c(ags[1],ags[2])
}
clusterSamples=TRUE
if(length(ags)>1){
  clusterSamples=as.logical(ags[2])
} 
cat(sprintf("Generating Clustering for %s\n",ags[1]))
data <- as.matrix(read.delim(file = ags[1],header = TRUE, sep = "\t",quote = "",blank.lines.skip = TRUE,stringsAsFactors = FALSE))
cat(sprintf("Read in a table containing %d rows and %d columns\n",dim(data)[1],dim(data)[2]))
if(startsWith(x=colnames(data)[1], prefix="Transcript.RepeatID"))
{
  cat("Assuming HOMER formated table\n")
  data2 = data[,9:dim(data)[2]]
  data2b = data[,9:dim(data)[2]]
  for( i in 1:dim(data)[1]) 
  {
    data[i,8]= strsplit(data[i,8],"[|]")[[1]][1]
    if( i %% 1000==0)
      cat(".")
  }
  data2= matrix(as.double(data2),nrow=dim(data2)[1],ncol=dim(data2)[2])
  rownames(data2) <- data[,8]
  colnames(data2) <- colnames(data[,9:dim(data)[2]])
  data3=data2
  if(min(data2)>=0)
  {
	  cat("Taking the log\n")
  	data3= log(data3[apply(data3, 1, function(x) var(x)>0),]+1,base=2)
  } else {
	  cat("Negative numbers detected so assuming already log transformed\n")
	  data3= data3[apply(data3, 1, function(x) var(x)>0),]
  }
  o=order(rowMeans(data3),decreasing = T)
  data4= data3[o[1:min(dim(data3)[1],2000)],]
  rownames(data4)<-rownames(data3)[o[1:min(dim(data3)[1],2000)]]
} else {
  cat("Assuming just one column for ids and one row for sample names\n")
  data2 = data[,2:dim(data)[2]]
  data2= matrix(as.numeric(data2),nrow=dim(data2)[1],ncol=dim(data2)[2])
  rownames(data2) <- data[,1]
  colnames(data2) <- colnames(data[,2:dim(data)[2]])
  data3=data2
  if(min(data2)>=0)
  {
    cat("Taking the log\n")
    data3= log(data3[apply(data3, 1, function(x) var(x)>0),]+1,base=2)
  } else {
    cat("Negative numbers detected so assuming already log transformed\n")
    data3= data3[apply(data3, 1, function(x) var(x)>0),]
  }
  o=order(rowMeans(data3),decreasing = T)
  data4= data3[o[1:min(dim(data3)[1],2000)],]
  rownames(data4)<-rownames(data3)[o[1:min(dim(data3)[1],2000)]]
}
rownames(data4)=make.unique(rownames(data4), sep = ".")
data3=data4

my_palette <- colorRampPalette(c("blue","white","red"))(n=299)
if(dim(data3)[2] > 20)
{
  my_palette <- colorRampPalette(c("blue","blue","white","red","red"))(n=299)
  pdf(sprintf("%s_cluster.pdf",ags[1]), width=8, height=12,title =sprintf("%s_cluster.pdf",ags[1]),pointsize = 8)

  if (dim(data3)[1] > 100) {
  	heatmap.2(data3, Colv = clusterSamples, main = sprintf("Top %d genes",dim(data3)[1]),trace="none", scale="row", col=my_palette,margins = c(12,12),cexRow = 0.2,cexCol=0.5,distfun=function(c) as.dist(1-cor(t(c))), hclustfun = function(x) hclust(x,method="ward.D2"),key.title = "") 
  } else {
  	heatmap.2(data3, Colv = clusterSamples, main = sprintf("Top %d genes",dim(data3)[1]),trace="none", scale="row", col=my_palette,margins = c(12,12),cexRow = 0.5,cexCol=0.5,distfun=function(c) as.dist(1-cor(t(c))), hclustfun = function(x) hclust(x,method="ward.D2"),key.title = "") 
  }
  dev.off()
} else {
  pdf(sprintf("%s_cluster.pdf",ags[1]), width=8, height=12,title =sprintf("%s_cluster.pdf",ags[1]),pointsize = 8)
  if(dim(data3)[1] < 100)
  {
    heatmap.2(data3, Colv = clusterSamples, main = sprintf("Top %d genes",dim(data3)[1]),  margins = c(10,10) ,scale="row",trace="none", col=my_palette,cexRow = 0.5,cexCol=0.5,distfun=function(c) as.dist(1-cor(t(c))),key.title = "",hclustfun = function(x) hclust(x,method="ward.D2"))
  }else {
    heatmap.2(data3, Colv = clusterSamples, main = sprintf("Top %d genes",dim(data3)[1]),  margins = c(10,10) ,scale="row",trace="none", col=my_palette,cexRow = 0.1,cexCol=0.5,distfun=function(c) as.dist(1-cor(t(c))),key.title = "", hclustfun = function(x) hclust(x,method="ward.D2"))
  }
  dev.off()
}
cat("Finished\n")

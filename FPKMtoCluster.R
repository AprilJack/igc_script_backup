library(gplots)
library(RColorBrewer)
ags<- commandArgs(trailingOnly = TRUE)
if (length(ags)==0){
  cat("Usage: fpkmCluster fpkm_file [clusters to output as separate text files] \n")
  ags=c("fpkm.txt",3)
}
if(length(ags)==1) {
  ags=c(ags[1],"3",2000)
}
if(length(ags)==2) {
      ags=c(ags[1],ags[2],2000)
}
if(length(ags)==3) {
  ags=c(ags[1],ags[2],as.numeric(ags[3]))
}
cat(sprintf("Generating Clustering for %s\n",ags[1]))
data <- as.matrix(read.delim(file = ags[1],header = TRUE, sep = "\t",quote = "",blank.lines.skip = TRUE,stringsAsFactors = FALSE))
cat(sprintf("Read in a table containing %d rows and %d columns\n",dim(data)[1],dim(data)[2]))
if(startsWith(x=colnames(data)[1], prefix="Transcript.RepeatID"))
{
  cat("Assuming HOMER formated fpkm table\n")
  data2 = data[,9:dim(data)[2]]
  data2b = data[,9:dim(data)[2]]
  for( i in 1:dim(data)[1]) 
  {
    data[i,8]= strsplit(data[i,8],"[|]")[[1]][1]
    if( i %% 1000==0)
      cat(".")
  }
  data2= matrix(as.double(data2),nrow=dim(data2)[1],ncol=dim(data2)[2])
  data2b= matrix(as.numeric(data2b),nrow=dim(data2b)[1],ncol=dim(data2b)[2])
  rownames(data2) <- data[,8]
  colnames(data2) <- colnames(data[,9:dim(data)[2]])
  rownames(data2b) <- data[,1]
  colnames(data2b) <- colnames(data[,9:dim(data)[2]])
  data3=data2
  data3b=data2b
  if(min(data2)>=0)
  {
	  cat("Taking the log\n")
  	data3= log(data3[apply(data3, 1, function(x) var(x)>0),]+1,base=2)
  	data3b= log(data3b[apply(data3b, 1, function(x) var(x)>0),]+1,base=2)
  } else {
	  cat("Negative numbers detected so assuming already log transformed\n")
	  data3= data3[apply(data3, 1, function(x) var(x)>0),]
	  data3b= data3b[apply(data3b, 1, function(x) var(x)>0),]
  }
  o=order(rowMeans(data3),decreasing = T)
  data4= data3[o[1:min(dim(data3)[1],as.numeric(ags[3]))],]
  data4b= data3b[o[1:min(dim(data3)[1],as.numeric(ags[3]))],]
  rownames(data4)<-rownames(data3)[o[1:min(dim(data3)[1],as.numeric(ags[3]))]]
  rownames(data4b)<-rownames(data3b)[o[1:min(dim(data3b)[1],as.numeric(ags[3]))]]
} else {
  cat("Assuming just one column for ids and one row for sample names\n")
  data2 = data[,2:dim(data)[2]]
  data2b = data[,2:dim(data)[2]]
  data2= matrix(as.numeric(data2),nrow=dim(data2)[1],ncol=dim(data2)[2])
  data2b= matrix(as.numeric(data2b),nrow=dim(data2b)[1],ncol=dim(data2b)[2])  
  rownames(data2) <- data[,1]
  colnames(data2) <- colnames(data[,2:dim(data)[2]])
  rownames(data2b) <- data[,1]
  colnames(data2b) <- colnames(data[,2:dim(data)[2]])
  data3=data2
  data3b=data2b
  if(min(data2)>=0)
  {
    cat("Taking the log\n")
    data3= log(data3[apply(data3, 1, function(x) var(x)>0),]+1,base=2)
    data3b= log(data3b[apply(data3b, 1, function(x) var(x)>0),]+1,base=2)
  } else {
    cat("Negative numbers detected so assuming already log transformed\n")
    data3= data3[apply(data3, 1, function(x) var(x)>0),]
    data3b= data3b[apply(data3b, 1, function(x) var(x)>0),]
  }
  o=order(rowMeans(data3),decreasing = T)
  data4= data3[o[1:min(dim(data3)[1],as.numeric(ags[3]))],]
  data4b= data3b[o[1:min(dim(data3)[1],as.numeric(ags[3]))],]
  rownames(data4)<-rownames(data3)[o[1:min(dim(data3)[1],as.numeric(ags[3]))]]
  rownames(data4b)<-rownames(data3b)[o[1:min(dim(data3b)[1],as.numeric(ags[3]))]]
}
rownames(data4)=make.unique(rownames(data4), sep = ".")
rownames(data4b)=make.unique(rownames(data4b), sep = ".")
data3=data4
data3b=data4b
if (length(ags) > 1)
{
  #user wants to get the top x clusters
  hc<-hclust(as.dist(1-cor(t(data3))),method="ward.D2")
  ct = cutree(hc,k=ags[2])
  clusters=cbind(ct,rownames(data3),data3)
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
    ct[ct==i]<-cColors[i]
  }
}

my_palette <- colorRampPalette(c("green","black","red"))(n=299)
if(dim(data3)[2] > 20)
{
  my_palette <- colorRampPalette(c("green","green","black","red","red"))(n=299)
  png(sprintf("%s_bigCluster.png",ags[1]), width=max(50*dim(data3)[2],3000), height=3000,res=300,pointsize=8)
  if (dim(data3)[1] > 100) {
  	heatmap.2(data3, main = sprintf("Top %d genes",dim(data3)[1]),trace="none", RowSideColors = ct,scale="row", col=my_palette,margins = c(12,12),cexRow = 0.25,cexCol=0.5,distfun=function(c) as.dist(1-cor(t(c))), hclustfun = function(x) hclust(x,method="ward.D2"),labRow=c(""),key.title = "") 
  } else {
  	heatmap.2(data3, main = sprintf("Top %d genes",dim(data3)[1]),trace="none", RowSideColors = ct,scale="row", col=my_palette,margins = c(12,12),cexRow = 0.5,cexCol=0.5,distfun=function(c) as.dist(1-cor(t(c))), hclustfun = function(x) hclust(x,method="ward.D2"),key.title = "") 
  }
  dev.off()
} else {
  png(sprintf("%s_cluster.png",ags[1]), width=1500, height=1500,res=300,pointsize=8)
  if(dim(data3)[1] < 100)
  {
    heatmap.2(data3, main = sprintf("Top %d genes",dim(data3)[1]), RowSideColors = ct,  margins = c(10,10) ,scale="row",trace="none", col=my_palette,cexRow = 0.5,cexCol=0.5,distfun=function(c) as.dist(1-cor(t(c))),key.title = "",hclustfun = function(x) hclust(x,method="ward.D2"))
  }else {
    heatmap.2(data3, main = sprintf("Top %d genes",dim(data3)[1]), RowSideColors = ct,  margins = c(10,10) ,scale="row",trace="none", col=my_palette,cexRow = 0.1,cexCol=0.5,distfun=function(c) as.dist(1-cor(t(c))),key.title = "",labRow = c(""), hclustfun = function(x) hclust(x,method="ward.D2"))
  }
  dev.off()
}
cat("Finished\n")

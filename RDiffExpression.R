if (!require("gplots")) {
  install.packages("gplots", repos="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(gplots)
}

if (!require("edgeR")) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("edgeR")
  library(edgeR)
}

if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", repos ="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(RColorBrewer)
}



#helper functions
topTwo <- function(x){
  x2= x %*% diag(1000000/colSums(x))
  x2= log(x2[apply(x2, 1, function(x3) !all(x3<100)&var(x3)>0),]+1,base=2)
  #rm <- rowMeans(x2)
  #x2 <- sweep(x2, 1, rm)
  #sx <- apply(x2, 1, sd)
  #x2 <- sweep(x2, 1, sx, "/")
  hc<-hclust(as.dist(1-cor(x2)),method="ward.D2")
  ct = cutree(hc,k=2)
  result= list()
  result[[1]]=as.matrix(x[,ct==1])
  result[[2]]=as.matrix(x[,ct==2])
  if(sum(ct==1)==1 | sum(ct==2)==1)
  {
    colnames(result[[1]])=colnames(x)[ct==1]
    colnames(result[[2]])=colnames(x)[ct==2]
  }
  return(result)
}


rDiff <- function(x,name){
  if(typeof(x) =="list")
  {
    samples1=combineTree(x[[1]])
    samples2=combineTree(x[[2]])
    cat(sprintf("Comparing clusters with %d and %d samples\n",dim(samples1)[2],dim(samples2)[2]))
    combined=cbind(samples1,samples2)
    dim1=1
    dim2=1
    if(!is.null(dim(samples1)))
    {
      dim1=dim(samples1)[2]
    }
    if(!is.null(dim(samples2)))
    {
      dim2=dim(samples2)[2]
    }
    if(dim1>1 & dim2>1)
    {
      group<-c(rep("A",dim1),rep("B",dim2))
      y <- DGEList(counts=combined, group=group)
      y <- calcNormFactors(y)
      y <- estimateDisp(y,design = model.matrix(~group))
      et <- exactTest(y)
      upgenes = et$table[et$table[,3]<0.001 & et$table[,1]< -1,]
      downgenes = et$table[et$table[,3]<0.001 & et$table[,1]> 1,]
      cat(sprintf("Cluster with %d vs %d there were %d genes up, and %d genes down (p-value < 0.001 and log2fold >1)\n",dim(samples1)[2],dim(samples2)[2],dim(upgenes)[1],dim(downgenes)[1]))
      combined = combined %*% diag(1000000/colSums(combined))
      uptable = combined[rownames(combined) %in% rownames(upgenes),]
      downtable = combined[rownames(combined) %in% rownames(downgenes),]
      colnames(uptable)=paste(group,c(colnames(samples1),colnames(samples2)),sep="_")
      colnames(downtable)=paste(group,c(colnames(samples1),colnames(samples2)),sep="_")
      write.table(uptable,file=sprintf("%s_%d_vs_%d_Up.txt",name,dim1,dim2),col.names = NA, row.names = TRUE,sep="\t",quote = FALSE)
      write.table(downtable,file=sprintf("%s_%d_vs_%d_Down.txt",name,dim1,dim2),col.names = NA, row.names = TRUE,sep="\t", quote=FALSE)
      forcluster = log(combined[apply(combined, 1, function(x3) !all(x3<100)&var(x3)>0),]+1,base=2)
      #rm <- rowMeans(forcluster)
      #forcluster <- sweep(forcluster, 1, rm)
      #forcluster <- apply(forcluster, 1, sd)
      #forcluster <- sweep(forcluster, 1, sx, "/")
      colnames(forcluster)=paste(group,c(colnames(samples1),colnames(samples2)),sep="_")
      my_palette <- colorRampPalette(c("green","black","red"))(n=299)
      if(dim(forcluster)[2] > 30)
      {
        my_palette <- colorRampPalette(c("green","green","black","red","red"))(n=299)
      }
      cat("Generating a heatmap from the log2(CPM+1)\n")
      png(sprintf("%s_cluster_%d_vs_%d.png",name,dim1,dim2), width=max(25*dim(forcluster)[2],1500), height=2000,res=300,pointsize=8)
      heatmap.2(forcluster, main = "Log2(CPM+1) Clustering",trace="none", col=my_palette,margins = c(12,12),cexRow = 0.25,cexCol=0.5,distfun=function(c) as.dist(1-cor(t(c))),scale="row", hclustfun = function(xx) hclust(xx,method="ward.D2"),labRow=c(""),key.title = "")
      dev.off();
      rDiff(x[[1]],paste(name,"A",sep=""))
      rDiff(x[[2]],paste(name,"B",sep=""))
    }
    else
    {
      cat(sprintf("%s had too few samples in one or both clusters: %d %d\n",name,dim1,dim2))
      if(dim1 > 1)
      {
        rDiff(x[[1]],paste(name,"A",sep=""))
      }
      if(dim2 > 1)
      {
        rDiff(x[[2]],paste(name,"B",sep=""))
      }
    }
  }
}

combineTree <-function(node)
{
  if(typeof(node)=="double")
  {
    return(as.matrix(node))
  } else {
    if(typeof(node)=="list")
    {
        return(cbind(combineTree(node[[1]]),combineTree(node[[2]])))
    } else {
      message("Weird type detected!\n")
    }
  }
}

rClust <- function(dataset,level=1)
{
  if(level == 0)
    return(dataset)
  #ok now lets run diff expression on these two clusters to see what the main genes are
  if(dim(dataset)[2]<4)
  {
    return(dataset)
  } else {
    result = list()
    two=topTwo(dataset)
    result[[1]]=rClust(two[[1]],level-1)
    result[[2]]=rClust(two[[2]],level-1)
    return(result)
  }
}

#ags<- commandArgs(trailingOnly = TRUE)
ags=c("raw_5-9-17.txt",6,FALSE)
if (length(ags)==0){
  cat("Usage: RDiffExpression raw_file clusterDepth [Simple format?TRUE/FALSE]\n")
  ags=c("raw.txt",4,FALSE)
}
cat(sprintf("Loading %s ...\n",ags[1]))
data <- as.matrix(read.delim(file = ags[1],header = TRUE, sep = "\t"))
if(length(ags)>2 & ags[3]==TRUE) #simple
{
  cat("Assuming just one column for ids and one row for sample names\n")
  data2 = data[,2:dim(data)[2]]
  data2= matrix(as.numeric(data2),nrow=dim(data2)[1],ncol=dim(data2)[2])
  rownames(data2) <- data[,1]
  colnames(data2) <- colnames(data[,2:dim(data)[2]])

} else {
  cat("Assuming HOMER formated raw table (otherwise set third arg to TRUE\n")
  data2 = data[,9:dim(data)[2]]
  for( i in 1:dim(data)[1]) 
  {
    data[i,8]= strsplit(data[i,8],"[|]")[[1]][1]
    if( i %% 1000==0)
      cat(".")
  }
  data2= matrix(as.numeric(data2),nrow=dim(data2)[1],ncol=dim(data2)[2])
  rownames(data2) <- data[,8]
  colnames(data2) <- colnames(data[,9:dim(data)[2]])

}
#ok now lets figure out the clustering for each cluster:)
cat("Determining clustering structure\n")
clustering=rClust(data2,as.double(ags[2]))
cat("Running Differential Expression and Outputing to Files...\n")
rDiff(clustering,paste(ags[1],"_",sep=""))
cat("Finished\n")


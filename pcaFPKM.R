#library(tm)
ags<- commandArgs(trailingOnly = TRUE)
#ags= c("fpkm_combined_normed.txt","")
if (length(ags)==0){
  cat("Usage: fpkmPCA fpkm_file [term1] [term2] ... [termN]\n")
  cat("Will color dots with prefixes you specify, otherwise will try to cluster")
  ags=c("fpkm.txt") 
  exit(1)
  
} 
library(Rtsne)
library(cluster)
library(gplots)
library(wordcloud)
terms=ags[2:length(ags)]
cat(sprintf("Generating PCA plot for %s",ags[1]))
data <- as.matrix(read.delim(file = ags[1],header = TRUE, sep = "\t",quote = "",blank.lines.skip = T,stringsAsFactors = F))
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
  data3= log(data2[apply(data2, 1, function(x) var(x)>0),]+1,base=2)
  data3b= log(data2b[apply(data2b, 1, function(x) var(x)>0),]+1,base=2)
  o=order(rowMeans(data3),decreasing = T)
  data3= data3[o[1:min(dim(data3)[1],1000)],]
  data3b= data3b[o[1:min(dim(data3)[1],1000)],]
  
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
  data3= log(data2[apply(data2, 1, function(x) var(x)>0),]+1,base=2)
  data3b= log(data2b[apply(data2b, 1, function(x) var(x)>0),]+1,base=2)
}
pr=prcomp(t(data3),scale = TRUE,center=TRUE)
variances<-pr$sdev^2
best=0
bestClusters=1
ncol(pr$x)
nrow(pr$x)
if(length(ags)==1)
{
	for(i in 2:min(12,nrow(pr$x)-1))
	{
	  km=kmeans(pr$x[,1:2],i,nstart = 100)
	  score=mean(silhouette(km$cluster,dist(pr$x[,1:2]))[,3])
	  if(score > best)
	  {
	    best=score
	    bestClusters=i
	  }
	}
	cat(sprintf("PCA best clusters is %d with score %3.2f\n",bestClusters,best))
	clust = kmeans(pr$x[,1:2],bestClusters,nstart = 100)
	cs = c("black","red","blue","green","cyan","orange","gray","brown") 
	colors=cs[clust$cluster]
} else {
 cat("Coloring by terms specified in the arguments\n")
 colors =rep("black",length(colnames(data3)))
 cs = c("red","blue","green","cyan","orange","gray","brown") 
 for( i in 1:length(terms))
 {
   colors[grep(terms[i],colnames(data3))]<-cs[i]	
 }
}
png(sprintf("%s_pca.png",ags[1]), width=1200, height=1200,res=300,pointsize=8)
#apply colored labels:
wordcloud::textplot(x=pr$x[,1],y=pr$x[,2],strtrim(colnames(data3),20),col=colors,xlim=c(min(pr$x[,1])-1,max(pr$x[,1])+1),ylim=c(min(pr$x[,2])-1,max(pr$x[,2])+1),cex=0.5,xlab=sprintf("PC1 (%3.2f%%)",100*(variances[1]/sum(variances))),ylab=sprintf("PC2 (%3.2f%%)",100*(variances[2]/sum(variances))), main="Plotting along first two principal components")
points(x=pr$x[,1],y=pr$x[,2],col=colors,xlim=c(min(pr$x[,1])-1,max(pr$x[,1])+1),ylim=c(min(pr$x[,2])-1,max(pr$x[,2])+1),cex=0.5,xlab=sprintf("PC1 (%3.2f%%)",100*(variances[1]/sum(variances))),ylab=sprintf("PC2 (%3.2f%%)",100*(variances[2]/sum(variances))))
dev.off()

png(sprintf("%s_pca_dots.png",ags[1]), width=1200, height=1200,res=300,pointsize=8)
#apply colored labels:
plot(x=pr$x[,1],y=pr$x[,2],pch=18,col=colors,xlim=c(min(pr$x[,1])-1,max(pr$x[,1])+1),ylim=c(min(pr$x[,2])-1,max(pr$x[,2])+1),cex=0.5,xlab=sprintf("PC1 (%3.2f%%)",100*(variances[1]/sum(variances))),ylab=sprintf("PC2 (%3.2f%%)",100*(variances[2]/sum(variances))),main="Plotting along first two principal components")
dev.off()
out=cbind(pr$x[,1],pr$x[,2])
rownames(out)=rownames(pr$x)
colnames(out)=c("PC1","PC2")
write.table(out,sprintf("%s_pca_dots.txt",ags[1]),sep="\t",col.names=NA,quote=FALSE)



set.seed(42)
#dim(data3)
tsne<-Rtsne(t(data3), perplexity = dim(data3)[2]/4,verbose = TRUE)
best=0
bestClusters=1
for(i in 2:min(12,nrow(tsne$Y)-1))
{
  km=kmeans(tsne$Y,i)
  score=mean(silhouette(km$cluster,dist(tsne$Y))[,3])
  if(score > best)
  {
    best=score
    bestClusters=i
  }
}
cat(sprintf("tSNE best clusters is %d with score %3.2f\n",bestClusters,best))
clust = kmeans(tsne$Y,bestClusters)
png(sprintf("%s_tsne.png",ags[1]), width=1200, height=1200,res=300,pointsize=8)
#plot(tsne$Y,pch=20,col="white",)
wordcloud::textplot(x=tsne$Y[,1],y=tsne$Y[,2],colnames(data3),xlim=c(min(tsne$Y[,1]-1),max(tsne$Y[,1])+1),ylim=c(min(tsne$Y[,2])-1,max(tsne$Y[,2])+1),col = colors,cex=0.5,main="tSNE" ,xlab="tSNE 1",ylab="tSNE 2")
dev.off()

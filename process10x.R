if (!require("cellrangerRkit")) {
  source("https://s3-us-west-2.amazonaws.com/10x.files/supp/cell-exp/rkit-install-1.1.0.R")
  library(cellrangerRkit)
  packageVersion("cellrangerRkit")
}

if (!require("edgeR")) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("edgeR")
  library(edgeR)
}
if (!require("preprocessCore")) {
  install.packages("preprocessCore", repos ="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(preprocessCore)
}
if (!require("Rtsne")) {
  install.packages("Rtsne", repos ="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(Rtsne)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", repos ="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(RColorBrewer)
}
if (!require("mygene")) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("mygene")
  library(mygene)
}

# Function to plot color bar
color.bar <- function(data,col=c("white","black"), xmin, xmax, y, nticks=11) {
  par(xpd=TRUE)
  min=min(data)
  max=max(data)
  pallete=colorRampPalette(col)(n=nticks)
  for (i in 1:nticks) {
    x = xmin+((i-1)*(xmax-xmin))/(nticks)
    rect(x,y,x+(xmax-xmin)/(nticks),y+2,col=pallete[i], border=NA)
    text(x+3,y-1,labels = sprintf("%3.1f",min+i*(max-min)/nticks),cex=0.8)
  }
}


ags<- commandArgs(trailingOnly = TRUE)
if (length(ags)==0){
  cat("Usage: process10x 10xDir [clusters=4] \n")
  ags=c(".",4)
}
cellranger_pipestance_path= ags[1]
cat(sprintf("Loading in %s\n",ags[1]))
gbm <- load_cellranger_matrix(cellranger_pipestance_path)
gbm2=as.matrix(exprs(gbm))
startingCells=dim(gbm2)[2]
startingGenes=dim(gbm2)[1]
colnames(gbm2)=colnames(gbm)
rownames(gbm2)=rownames(gbm)
quantileThreshold=quantile(colSums(gbm2),probs=0.3)
gbm2 <- gbm2[,(colSums(gbm2)>quantileThreshold)&gbm2[which(rownames(gbm)=="ENSMUSG00000070570"),]>0]
goodCells=dim(gbm2)[2]
cat(sprintf("Filtering cells in %s. Went from %g to %g cells which had counts greater than %g\n",ags[1],startingCells,goodCells,quantileThreshold))
gbm3 = gbm2 %*% diag(1000000/colSums(gbm2))
gbm3 = normalize.quantiles(gbm3)
colnames(gbm3)=colnames(gbm2)
rownames(gbm3)=rownames(gbm2)
#vars = apply(gbm3,1,var)
#varThreshold=quantile(vars,probs=0.99)
data <- log(gbm3[apply(gbm3, 1, function(x) (sum(x>5)>dim(gbm3)[2]/10)&var(x)>0),]+1,base=2)

goodGenes=dim(data)[1]
cat(sprintf("Filtering genes in %s. Went from %g to %g genes which had more than 5 normalized counts in %g cells\n",ags[1],startingGenes,goodGenes,dim(gbm3)[2]/10))
set.seed(42)
tsne<-Rtsne(t(data), dims = 2, perplexity = 10,
            theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 1000,
            verbose = FALSE, is_distance = FALSE)
km=kmeans(tsne$Y,centers = as.integer(ags[2]))
colors=c("red","blue","green","yellow","cyan","orange","purple","brown","gray","black")
for(i in 1:dim(tsne$Y)[1])
{
  tsne$Color[i]=colors[km$cluster[i]]
}
png(sprintf("TSNE_with_%s_clusters.png",ags[2]), width=2000, height=2000,res=300,pointsize=8)
plot(tsne$Y,pch=19,col=tsne$Color,xlab="tSNE 1",ylab="tSNE 2")
legend("topright",legend=as.character(1:as.integer(ags[2])),col=colors,pch = 19)
dev.off()
#get diff expression for all groups

#lets do edgeR between clusters
for(i in 1:as.integer(ags[2]))
{
  group = factor(km$cluster==i)
  design= model.matrix(~group)
  cat(sprintf("Finding genes specific to cluster %s\n",colors[i]))
  y <- DGEList(counts=gbm2, group=group)
  y <- calcNormFactors(y)
  y <- estimateDisp(y,design)
  norm.counts <- cpm(y)
  et <- exactTest(y)
  upgenes = et$table[et$table[,3]<0.05 & et$table[,1]> 1,]
  downgenes = et$table[et$table[,3]<0.05 & et$table[,1]< -1,]
  cat(sprintf("Annotating differential genes in cluster %s\n",colors[i]))
  if(startsWith(rownames(upgenes)[1],"ENSMUSG"))
  {
    upannotation=queryMany(rownames(upgenes),scopes=("ensemblgene"),species="mouse",returnall=TRUE)$response
    downannotation=queryMany(rownames(downgenes),scopes=("ensemblgene"),species="mouse",returnall=TRUE)$response
  } else {
    upannotation=queryMany(rownames(upgenes),scopes=("ensemblgene"),species="human",returnall=TRUE)$response
    downannotation=queryMany(rownames(downgenes),scopes=("ensemblgene"),species="human",returnall=TRUE)$response
  }
  
  upmatrix=norm.counts[rownames(norm.counts)%in%rownames(upgenes),]
  colnames(upmatrix)=paste(group,colnames(upmatrix),sep="_")
  downmatrix=norm.counts[rownames(norm.counts)%in%rownames(downgenes),]
  colnames(downmatrix)=paste(group,colnames(downmatrix),sep="_")
  upmatrix= upmatrix[, order(colnames(upmatrix))]
  downmatrix= downmatrix[, order(colnames(downmatrix))]
  write.table(file=sprintf("Cluster_%s_Up.txt",colors[i]),cbind(upgenes,upannotation,upmatrix),row.names=TRUE,col.names =NA, sep="\t")
  write.table(file=sprintf("Cluster_%s_Down.txt",colors[i]),cbind(downgenes,downannotation,downmatrix),row.names=TRUE,col.names =NA, sep="\t")
  cat(sprintf("Cluster %s had %d genes up, and %d genes down (p-value < 0.05 and log2fold >1)\n",colors[i],dim(upgenes)[1],dim(downgenes)[1]))
  for(g in 1:dim(upmatrix)[1])
  {
    cat(sprintf("Generating expression plot for gene %s\n",rownames(upgenes)[g]))
    maxCount=max(upmatrix[g,])
    minCount=min(upmatrix[g,])
    my_palette <- colorRampPalette(c("white",colors[i]))(n=maxCount-minCount+1)
    for(c in 1:dim(tsne$Y)[1])
    {
      tsne$Exp[c]=my_palette[upmatrix[g,c]-minCount+1]
    }
    png(sprintf("%s_%s.png",colors[i],rownames(upmatrix)[g]), width=1200, height=1200,res=300,pointsize=8)
    plot(tsne$Y,pch=21,cex=0.7,bg=tsne$Exp,col="gray",xlab="tSNE 1",ylab="tSNE 2",main=sprintf("%s expression",rownames(upmatrix)[g]))
    color.bar(data=upmatrix[g,],col=c("white",colors[i]),xmin = min(tsne$Y[,1]),xmax=max(tsne$Y[,1]),y=max(tsne$Y[,2])+2)
    dev.off()
  }
}

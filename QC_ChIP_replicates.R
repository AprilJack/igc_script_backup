#!/usr/bin/env Rscript
# plot the PCA and heatmap for the downstream of Quantify_ChIP_Replicates.sh
# special points about this script: use norm2total for ChIP data.

# assume I will automatically load CalcPCA function
options(stringsAsFactors=F)
#source("/gpfs/analyses/ling/Rcode/PCA.R")
library(ggplot2)
library(ggrepel)
library(cowplot)
args<-commandArgs(trailingOnly=T)

if(length(args)!=1) {
  print("please provide the results of Quantify_ChIP_Replicates.sh as the only input for this script")
  stop()
}

if(!file.exists(args[1])) {
  print("input file does not exist")
  stop()
}


# put the two functions here so that the script is no longer dependent on the PCA.R

CalcPCA<-function(df,top_ratio=0.1,logTransform=F) {
    # if top_ratio > 1, I think it refers to the absolute number of genes to be used
  if(sum(grepl("ANNO",names(df)))==1) {
    rownames(df)<-df$ANNO
    df<-df[,-(which(names(df)=="ANNO"))]
  }
    if(logTransform) {
        dflog2<-log2(df+5)
        df<-dflog2
    }
    if(top_ratio > 1) {
        ntop<-top_ratio
    } else {
        ntop=round(nrow(df)*top_ratio,0)
    }
    # choose the candidates based on hard expression cutoff, not used anymore
    # exprcutoff<-quantile(rowSums(dflog2),prob=1-top_ratio)
    # dflog2expressed<-dflog2[rowSums(dflog2)>=exprcutoff,]
    # numexpr=nrow(dflog2expressed)
  # choose the top candidates based on variation...DESeq2 setting
  rv <- genefilter::rowVars(df)
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  topdf<-df[select,]
  PCA=prcomp(t(topdf),center=T,scale.=F)
  return(PCA)
}

PlotPCA<-function(PCA,group,top_ratio) {
    # if top_ratio > 1, I think it refers to the absolute number of genes to be used
    if(top_ratio > 1) {
        main=paste("PCA on top ",top_ratio," most variable genes",sep='')
    } else {
        main=paste("PCA on top ",top_ratio*100,"% most variable genes",sep='')
    }
  PC1=round(PCA$sdev[1]^2/sum(PCA$sdev^2)*100,1)
  PC2=round(PCA$sdev[2]^2/sum(PCA$sdev^2)*100,1)
  xlab=paste("PC1   ",PC1,"% of total variance",sep='')
  ylab=paste("PC2   ",PC2,"% of total variance",sep='')
  dfPCA<-data.frame(PC1=PCA$x[,1],PC2=PCA$x[,2],Group=group,label=rownames(PCA$x))
  p<-ggplot(dfPCA,aes(x=PC1,y=PC2,col=Group,text=label,label=label)) + geom_point() + theme_bw() +
    xlab(xlab) + ylab(ylab) + ggtitle(main)
  return(p)
}



df<-read.table(args[1],sep="\t",header=T,quote="",comment.char="")
# the first 19 columns will be from annotatedPeaks.pl

# calculate the normDepth factor
TotalReads<-colnames(df)[grepl("Tag.Count.in.given.bp",colnames(df))]
nSample<-length(TotalReads)
sampleName<-sapply(TotalReads,function(x) strsplit(x,".Tag.Count.in.given.bp",fixed=T)[[1]][1])
TotalReads<-sapply(TotalReads,function(x) strsplit(x,".Tag.Count.in.given.bp..",fixed=T)[[1]][2])
TotalReads<-sapply(TotalReads,function(x) strsplit(x,".Total..normalization.factor",fixed=T)[[1]][1])
TotalReads<-as.numeric(TotalReads)
# normDepth<-TotalReads/mean(TotalReads)

dfmatrix<-df[,20:ncol(df)]
colnames(dfmatrix)<-sampleName
rownames(dfmatrix)<-df[,1]

# calculate simpleNorm counts per million total mapped reads
dfmatrixNorm<-t(apply(dfmatrix,1,function(x) x/TotalReads*1000000))
dfmatrixNormLog<-log2(dfmatrixNorm+5)

dfpca<-CalcPCA(dfmatrixNormLog,top_ratio = 500,logTransform = F)
dfpca.p<-PlotPCA(dfpca,group=sampleName,top_ratio=500) + ggrepel::geom_text_repel() + theme_cowplot()

outName<-gsub(".txt",".jpeg",args[1])
ggsave(outName,dfpca.p,width=30,height=10,units="cm",dpi=600)
print("done")


options(width=160)
ags<- commandArgs(trailingOnly = TRUE)
#ags= c("raw.txt","design.txt")
if (length(ags)==0){
  cat("Usage: getDiffR raw.txt sampleSheet.txt [design formula] [fdr=0.05] [log2fold=1] [contrasts_file.txt]\n
      Checks for column labels in sampleSheet. If it sees \"Batch\" it will get all pairwise diffs with additive batch correction using formula ~0+Batch+Condition
      Otherwise expects a column labeled \"Condition\" and will run the formula ~0+Condition.
      Custom design formuals are possible by specifying on arg3. Then it will ask for manual user contrast entry.\n
      If you specify a contrasts_file as the 6th arg, then getDiffR will perform those comparisons automatically.\n
      The comparison list should have one line per comparison with the format \"A-B\" or \"1 -1 0 0\"")
  quit()
}
if (!require("BiocParallel")) {
  install.packages("BiocParallel", dependencies= TRUE,quiet = TRUE)
  library("BiocParallel",quietly = TRUE,warn.conflicts = FALSE)
}
if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies= TRUE,quiet=TRUE)
  library(gplots,quietly = TRUE,warn.conflicts = FALSE)
}
if (!require("ggfortify")) {
  install.packages("ggfortify", dependencies= TRUE,quiet=TRUE)
  library(ggfortify,quietly = TRUE,warn.conflicts = FALSE)
}

if (!require("DESeq2")) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("DESeq2")
  library(DESeq2,quietly = TRUE, warn.conflicts = FALSE)
}
if(!require("limma")){
  BiocManager::install("limma")
  library("limma",quietly = TRUE, warn.conflicts = FALSE)  
}
data <- as.matrix(read.delim(file = ags[1],header = TRUE, sep = "\t",quote = "",blank.lines.skip = T,stringsAsFactors = F))
annot = NULL
if(startsWith(x=colnames(data)[1], prefix="Transcript.RepeatID"))
{
  cat("Assuming HOMER formated raw count table\n")
  data2 = data[,9:dim(data)[2]]
  for( i in 1:dim(data)[1]) 
  {
    data[i,8]= strsplit(data[i,8],"[|]")[[1]][1]
    if( i %% 1000==0)
      cat(".")
  }
  data2= matrix(as.integer(data2),nrow=dim(data2)[1],ncol=dim(data2)[2])
  rownames(data2) <- data[,8]
  colnames(data2) <- colnames(data[,9:dim(data)[2]])
  dim(data2)
  annot=data[,1:7]

} else {
  cat("Assuming just one column for ids and one row for sample names\n")
  data2 = data[,2:dim(data)[2]]
  data2= matrix(as.integer(as.matrix(data2)),nrow=dim(data2)[1],ncol=dim(data2)[2])
  rownames(data2) <- data[,1]
  colnames(data2) <- colnames(data[,2:dim(data)[2]])
  dim(data2)
}
register(MulticoreParam(16))
sample=as.matrix(read.delim(file = ags[2],header = TRUE, sep="\t",quote="",blank.lines.skip = TRUE,stringsAsFactors = FALSE))
print(cbind(colnames(data2),sample))
customFormula=TRUE
if(length(ags)<3)
{
  #We are doing default design (e.g. look out for Batch )
  if("Batch" %in% colnames(sample))
  {
    ags[3]="~0+Condition+Batch"
  } else {
    ags[3]="~0+Condition"
  }
  customFormula=FALSE
} 
if(length(ags)<4) {
  fdr=0.05
} else {
  fdr <- as.numeric(ags[4])
}
if(length(ags)<5) {
  log2fold=1
} else {
  log2fold <- as.numeric(ags[5])
}
ds = as.formula(ags[3])
dds <- DESeqDataSetFromMatrix(data2, sample, design = ds)
cat("Estimating size factors...\n")
dds <- estimateSizeFactors(dds)
idx <- rowSums( counts(dds, normalized=TRUE) >= 10 ) >= min(1,dim(data2)[2]/4)
if(!is.null(annot)){
  annot=annot[idx,]
}
cat(sprintf("Filter all but the top %d genes from the original %d\n",sum(idx),dim(data2)[1]))
dds <- dds[idx,]
cat("Running DESeq2\n")
dds <- DESeq(dds,parallel=TRUE)
normed <- counts(dds, normalized=TRUE)
rnames=substring(resultsNames(dds),first = nchar(colnames(sample)[2])+1,last=99)
big=NULL
if(!customFormula)
{

  cat(sprintf("Running DESeq2 on %s using samplesheet %s and %s design with fdr %3.3f and log2fold %3.3f\n",ags[1],ags[2],ags[3],fdr,log2fold))
  for(i in 1:(length(rnames)-1))
  {
    for(j in (i+1):length(rnames))
    {
      cat(sprintf("%s vs %s: ",rnames[i],rnames[j]))
      x=vector(mode="numeric",length=length(rnames))
      x[i]=1
      x[j]=-1
      cat("Testing...")
      res <- results(dds,contrast =x, parallel=TRUE, alpha=fdr)
      if(dim(res)[1]>0)
      {
        big=cbind(big,res$baseMean, res$log2FoldChange,res$padj)
        if(length(colnames(big))==0) colnames(big)=c("","","")
        colnames(big)[(length(colnames(big))-2):length(colnames(big))]=c(sprintf("%s_vs_%s",rnames[i],rnames[j]),"log2fold","adj.p")
        rownames(big)=rownames(res)
        up <- res[!is.na(res$padj)&res$padj<fdr&res$log2FoldChange>log2fold,]
        up <- up[order(up$padj,decreasing = FALSE),]
        down <- res[!is.na(res$padj)&res$padj<fdr&res$log2FoldChange < -log2fold,]
        down <- down[order(down$padj,decreasing = FALSE),]
        png(sprintf("%s_vs_%s_Volcano.png",rnames[i],rnames[j]), width=2000, height=2000,res=300,pointsize=8)
        plot(res$log2FoldChange,-log10(res$padj),pch=20,main=sprintf("%s_vs_%s",rnames[i],rnames[j]),cex=0.6,xlab="log2FC",ylab="-log10(p.adj)",col=alpha("black",0.5),bg=alpha("black",0.5))
        points(up$log2FoldChange,-log10(up$padj),pch=20,col=alpha("red",0.5),bg=alpha("red",0.5),cex=0.7)
        points(down$log2FoldChange,-log10(down$padj),pch=20,col=alpha("green",0.5),bg=alpha("green",0.5),cex=0.7)
        legend("bottomright",sprintf("Up:%d Down:%d",dim(up)[1],dim(down)[1]))
        dev.off()
        cat(sprintf("UP:%d DOWN:%d\n",dim(up)[1],dim(down)[1]))
        if(dim(up)[1]>0) {
          colnames(up)[6]=sprintf("UP %s vs %s p.adj",rnames[i],rnames[j])
          write.table(x = up,file = sprintf("UP_%s_vs_%s.txt",rnames[i],rnames[j]),quote = FALSE,sep = "\t",col.names = NA)
        }
        if(dim(down)[1]>0) {
          colnames(down)[6]=sprintf("DOWN%s vs %s p.adj",rnames[i],rnames[j])
          write.table(x = down,file = sprintf("DOWN_%s_vs_%s.txt",rnames[i],rnames[j]),quote = FALSE,sep = "\t",col.names=NA)
        }
      }
    }
  }
} else {
  if (length(ags)>5)
  {
    cat(sprintf("Running DESeq2 on %s using samplesheet %s and manual design %s with fdr %3.3f and log2fold %3.3f with comparisons taken from %s\n",ags[1],ags[2],ags[3],fdr,log2fold,ags[6]))
    comparisons=read.table(ags[6],header = FALSE,sep = " ",blank.lines.skip = TRUE,stringsAsFactors = FALSE)
    for(i in 1:dim(comparisons)[1])
    {
      x=trimws(comparisons[i,])
      x = as.numeric(x)
      cat("Testing...")
      cat(x)
      #if(!grep(x,pattern ="[A-Za-z]"))
      #{
      #  x = as.numeric(unlist(strsplit(userString," ")))
      #  res <- results(dds,contrast =x, parallel=TRUE,alpha=fdr)
      #  x=paste(x[x != 0],rnames[x!=0],collapse = "_",sep = "")
      #} else {
        res <- results(dds,contrast =x, parallel=TRUE,alpha=fdr)
      #}
      cat("Done...")
      if(dim(res)[1]>0)
      {
        big=cbind(big,res$baseMean, res$log2FoldChange,res$padj)
        if(length(colnames(big))==0) colnames(big)=c("","","")
        colnames(big)[(dim(big)[2]-2):dim(big)[2]]=c(sprintf("%s",x),"log2fold","adj.p")
        rownames(big)=rownames(res)
        up <- res[!is.na(res$padj)&res$padj<fdr&res$log2FoldChange>log2fold,]
        if(length(up)>0) {
          up <- up[order(up$padj,decreasing = FALSE),]
        }
        down <- res[!is.na(res$padj)&res$padj<fdr&res$log2FoldChange < -log2fold,]
        if(length(down)>0) {
          down <- down[order(down$padj,decreasing = FALSE),]
        }
        cat("Plotting...")
        png(sprintf("Volcano_%s.png",x), width=2000, height=2000,res=300,pointsize=8)
        plot(res$log2FoldChange,-log10(res$padj),pch=20,main=sprintf("%s",x),cex=0.6,xlab="log2FC",ylab="-log10(p.adj)",col=alpha("black",0.5),bg=alpha("black",0.5))
        points(up$log2FoldChange,-log10(up$padj),pch=20,col=alpha("red",0.5),bg=alpha("red",0.5),cex=0.7)
        points(down$log2FoldChange,-log10(down$padj),pch=20,col=alpha("green",0.5),bg=alpha("green",0.5),cex=0.7)
        legend("bottomright",sprintf("Up:%d Down:%d",dim(up)[1],dim(down)[1]))
        dev.off();
        cat(sprintf("UP:%d DOWN:%d %s\n",dim(up)[1],dim(down)[1],x))
        if(dim(up)[1]>0) {
          colnames(up)[6]=sprintf("UP_%s p.adj",x)
          write.table(x = up,file = sprintf("UP_%s.txt",x),quote = FALSE,sep = "\t",col.names = NA) 
        }
        if(dim(down)[1]>0) {
          colnames(down)[6]=sprintf("DOWN_%s p.adj",x)
          write.table(x = down,file = sprintf("DOWN_%s.txt",x),quote = FALSE,sep = "\t", col.names=NA)
        }
      }
    }
  } else {
    cat(sprintf("Running DESeq2 on %s using samplesheet %s and manual design %s with fdr %3.3f and log2fold %3.3f\n",ags[1],ags[2],ags[3],fdr,log2fold))
    while(TRUE)
    {
      print(rnames)
      x=c()
      while(TRUE)
      {
        cat("Please type in a contrast vector (e.g. 1 -1 0 0) +Ctrl-D or empty+Ctrl-D to exit ")
        userString <- readLines("stdin",n=1)
        if(nchar(userString)==0)
        {
          break
        }
        x = as.numeric(unlist(strsplit(userString," ")))
        if(length(x)==length(rnames)) {
          break
        } else cat(sprintf("Please provide a contrast vector with %d elements!\n",length(rnames)))
      }
      if(length(x) != length(rnames))
      {
        break
      }
      cat("Testing...")
      res <- results(dds,contrast =x, parallel=TRUE,alpha=fdr)
      cat("Done...")
      x=paste(x[x != 0],rnames[x!=0],collapse = "_",sep = "")
      if(dim(res)[1]>0)
      {
        big=cbind(big,res$baseMean, res$log2FoldChange,res$padj)
        if(length(colnames(big))==0) colnames(big)=c("","","")
        colnames(big)[(dim(big)[2]-2):dim(big)[2]]=c(sprintf("%s",x),"log2fold","adj.p")
        rownames(big)=rownames(res)
        up <- res[!is.na(res$padj)&res$padj<fdr&res$log2FoldChange>log2fold,]
        if(length(up)>0) {
          up <- up[order(up$padj,decreasing = FALSE),]

        }
        down <- res[!is.na(res$padj)&res$padj<fdr&res$log2FoldChange < -log2fold,]
        if(length(down)>0) {
          down <- down[order(down$padj,decreasing = FALSE),]

        }
        cat("Plotting...")
        png(sprintf("Volcano_%s.png",x), width=2000, height=2000,res=300,pointsize=8)
        plot(res$log2FoldChange,-log10(res$padj),pch=20,main=sprintf("%s",x),cex=0.6,xlab="log2FC",ylab="-log10(p.adj)",col=alpha("black",0.5),bg=alpha("black",0.5))
        points(up$log2FoldChange,-log10(up$padj),pch=20,col=alpha("red",0.5),bg=alpha("red",0.5),cex=0.7)
        points(down$log2FoldChange,-log10(down$padj),pch=20,col=alpha("green",0.5),bg=alpha("green",0.5),cex=0.7)
        legend("bottomright",sprintf("Up:%d Down:%d",dim(up)[1],dim(down)[1]))
        dev.off();
        cat(sprintf("UP:%d DOWN:%d %s\n",dim(up)[1],dim(down)[1],x))
        if(dim(up)[1]>0) {
          colnames(up)[6]=sprintf("UP_%s p.adj",x)
          write.table(x = up,file = sprintf("UP_%s.txt",x),quote = FALSE,sep = "\t",col.names = NA) 
        }
        if(dim(down)[1]>0) {
          colnames(down)[6]=sprintf("DOWN_%s p.adj",x)
          write.table(x = down,file = sprintf("DOWN_%s.txt",x),quote = FALSE,sep = "\t", col.names=NA)
        }
      }
    }
  }
}
cat("Done! Saving big tables and version info.\n")
v=sessionInfo()
write.table(unlist(v),file=sprintf("Versions_%s_%s_%s",ags[1],ags[2],ags[3]))
if(is.null(annot)){
  write.table(x = cbind(big,normed),file = sprintf("Diff_%s_%s_%s.txt",ags[1],ags[2],ags[3]),quote = FALSE,sep = "\t",col.names = NA) 
} else {
  write.table(x = cbind(annot,big,normed),file = sprintf("Diff_%s_%s_%s.txt",ags[1],ags[2],ags[3]),quote = FALSE,sep = "\t",col.names = NA) 
}
write.table(x = normed,file = sprintf("Normed_%s_%s_%s.txt",ags[1],ags[2],ags[3]),quote = FALSE,sep = "\t",col.names = NA)



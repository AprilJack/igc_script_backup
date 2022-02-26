#!/usr/bin/env Rscript
# use DESeq2 to call the DE peaks
# NOTE: it does not consider the input background ChIP. it is for spike-in DE analysis only.
# try to add the batch effect correction

options(stringsAsFactors=F)
args <- commandArgs(trailingOnly = TRUE)


if(length(args)!=6) {
	print("command line arguments:")
	print("1. homer quantified peak raw counts table <annotatePeaks_noadj_peaks.txt>")
	print("2. results from CalcNormFactor.R <spikeIn_raw.tsv>")
	print("3. design table <design_table.txt>")
	print("4. name for the output <A_vs_B>")
	print("5. DE cutoff of fdr: <0.05>")
	print("6. DE cutoff of logFC: <1>")
	print("Usage example: DEpeaks_spikeIn.R annotatePeaks_noadj_peaks.txt spikeIn_raw.tsv design_table.txt A_vs_B 0.05 1")
	print("The ratio of spike-in reads to total mapped reads will be used as norm.factors for all samples. The library size will also be adjusted based on the total mapped reads per sample.")
	print("Because the normFactor is a relative scaling for library size. It changes with the input samples. So the norm factor is specific to each design_table.")
	stop()
}

infile<-args[1]
if(!file.exists(infile)) {
	stop("input file does not exists")
}
df<-read.table(infile,sep="\t",header=T,quote="",comment.char="")
# infile="tags/merged_pooled_H12ctrl_and_H12_93KD_test_regions.txt"
# remove the last line that contains "done", a bug has been corrected
# df<-df[-82127,]

spikeInFile<-args[2]
if(!file.exists(spikeInFile)) {
	stop("spike-in results not found")
}
dfspikein<-read.table(spikeInFile,sep="\t",header=T,quote="",comment.char="")

designTable<-args[3]
if(!file.exists(designTable)) {
	stop("design table not found")
}
dfdesign<-read.table(designTable,sep="\t",header=F,quote="",comment.char="")

outName<-args[4]

fdr<-as.numeric(args[5])
if(fdr>0.2 | fdr < 0) {
	stop("Please make sure the fdr is at the range (0,0.2)")
}

logFC<-as.numeric(args[6])
if(logFC<0.3 | logFC > 5) {
	stop("Please make sure the logFC is at the range (0.3,5)")
}


######################################## functions #######################
library(VennDiagram)
library(ggrepel)
GetSampleNameDesignTable<-function(oneDir) {
	if(grepl("/",oneDir,fixed=T)) {
		sampleName<-rev(strsplit(oneDir,split="\\/")[[1]])[1]
	} else {
		sampleName<-oneDir
	}
	sampleName<-gsub(".fastq.gz","",sampleName)
	sampleName<-gsub(".fastq","",sampleName)
	return(sampleName)
}

GetSampleNamePeak<-function(oneName,samples) {
	# we do not need to grab total mapped reads now
	sampleStrings<-strsplit(oneName,"\\.")[[1]]
	sampleName<-sampleStrings[sampleStrings %in% samples]
	# it is possible that this sample is not inlcuded in the design table
	if(length(sampleName)<1) {
		sampleName=paste0("NotUsed_",runif(1,min=1,max=1000))
	}
	totalRead<-sampleStrings[which(sampleStrings=="Total")-2]
	return(c(sampleName,totalRead))
}

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
    xlab(xlab) + ylab(ylab) + ggtitle(main) + geom_text_repel(size=2)
  return(p)
}

Plot_Venn<-function(alist,filename,pal="default",CategoryNames=F) {
	# # area1 must be greater than area2, which must be greater than area3
	if(pal=="default") {
		pal=c("cornflowerblue","grey50")
	}
	if(CategoryNames==F) {
		CN = rep("",length(alist))
	} else if (CategoryNames==T) {
		CN = names(alist)
	}
  venn.diagram(alist,filename,imagetype="png",width=10,height=10,units="cm",resolution=600,fontfamily="Arial",
	             main.fontfamily="Arial",sub.fontfamily = "Arial",cat.fontfamily="Arial",alpha=0.5,col="transparent",category.names=CN,
               fill=pal,margin=0.1)
}

######################################### end of functions ####################################

dfdesign$sample<-sapply(dfdesign$V1,GetSampleNameDesignTable)
colnames(dfspikein)<-sapply(colnames(dfspikein),function(x) rev(strsplit(x,".",fixed=T)[[1]])[1])



df.part1<-df[,c(1:19)]
df.part2<-df[,c(20:ncol(df))]
rownames(df.part2)<-df.part1[,1]
rownames(df.part1)<-df.part1[,1]
dfInfo<-as.data.frame(t(sapply(colnames(df.part2),function(x) GetSampleNamePeak(x,dfdesign$sample))))
rownames(dfInfo)<-dfInfo[,1]
colnames(dfInfo)<-c("Sample","TotalReads")
dfInfo$TotalReads<-as.numeric(dfInfo$TotalReads)
colnames(df.part2)<-dfInfo$Sample
dfInfo$PeakReads<-colSums(df.part2)


# 1. reorganize the column order to match the one in the design table
if(sum(colnames(df.part2)%in%dfdesign$sample)!=length(dfdesign$sample)) {
	stop("Please make sure the samples included in the design table are provided in the input count table")
}
df.part2<-df.part2[,dfdesign$sample]
dfInfo<-dfInfo[dfdesign$sample,]
dfspikein<-dfspikein[,dfdesign$sample]
group<-factor(dfdesign$V2)
group<-relevel(group,ref=dfdesign$V2[1]) # always use the first factor as reference level
if(ncol(dfdesign)>3) {
	Batch=T
} else {
	Batch=F
}

# 1.5 calculate the norm factor 
dfInfo$SpikeInReads<-colSums(dfspikein)
dfInfo$normDepth<-dfInfo$TotalReads/mean(dfInfo$TotalReads) # norm factor for sequencing depth
dfInfo$SpikeInRatio<-dfInfo$SpikeInReads/dfInfo$TotalReads # ratio of spike-in to total mapped reads, used for IP efficiency norm.
dfInfo$normSpikeIn<-dfInfo$SpikeInRatio/mean(dfInfo$SpikeInRatio) # norm factor for spike-in (IP efficiency)

# 2. DE analysis
################ DESeq2 #########################
library(DESeq2)
library(grid)
library(gridExtra)
design_edgeR_to_deseq2<-function(rawCount,group_factor,Batch=FALSE) {
	# Batch is either FALSE or a vector of the same length of group_factor, indicating the batches
	if(identical(Batch,FALSE)) {
		colData<-data.frame(condition=group_factor)
	} else {
		colData<-data.frame(condition=group_factor,batch=Batch)
	}
  rownames(colData)<-colnames(rawCount)
  return(colData)
}

if(Batch) {
	batchFactor<-factor(dfdesign$V3)
	batchFactor<-relevel(batchFactor,ref=dfdesign$V3[1])
	colData<-design_edgeR_to_deseq2(round(df.part2,0),group,Batch=batchFactor)
	dds<-DESeqDataSetFromMatrix(countData=round(df.part2,0),colData=colData,design=~batch+condition)
} else {
	colData<-design_edgeR_to_deseq2(round(df.part2,0),group,Batch=F)
	dds<-DESeqDataSetFromMatrix(countData=round(df.part2,0),colData=colData,design=~condition)
}

# NOTE: DESeq(dds) actually has three steps: estimateSizeFactors(dds); estimateDispersions(dds); nbinomWaldTest(dds)
# sizeFactors(dds)<- dfInfo$normDepth # this option produces identical result as homer getDifferentialPeaksReplicates.pl
sizeFactors(dds)<- dfInfo$normDepth * dfInfo$normSpikeIn
# they have extra 470 DE down peaks, include all 504 DE down peaks of simple normDepth, make sense.
dds<-estimateDispersions(dds)
dds<-nbinomWaldTest(dds)
res<-results(dds)
deseq2<-as.data.frame(apply(res,2,as.vector))
rownames(deseq2)<-rownames(res)

deseq2$sig<-ifelse(!is.na(deseq2$padj)&deseq2$padj<fdr&abs(deseq2$log2FoldChange)>logFC,1,0)
# identical(rownames(df.part1),rownames(deseq2)) ## TRUE
deseq2<-cbind(df.part1,deseq2)
# also add the normalized reads
identical(rownames(deseq2),rownames(counts(dds,normalized=T,replaced=F))) # TRUE
deseq2<-cbind(deseq2,counts(dds,normalized=T,replaced=F))
write.table(deseq2,paste0(outName,"_DiffPeaks_result.tsv"),sep="\t",row.names=F,quote=F)
DEdeseq2<-subset(deseq2,sig==1)
if(nrow(DEdeseq2)>0) {
	write.table(DEdeseq2,paste0(outName,"_DiffPeaks_DEresult.tsv"),sep="\t",row.names=F,quote=F)
} else {
	print("No DE peaks found at current cutoff")
}



# what about without spikein?
# sizeFactors(dds)<- dfInfo$normDepth # this option produces identical result as homer getDifferentialPeaksReplicates.pl
ddsNS<-dds
sizeFactors(ddsNS)<-dfInfo$normDepth
ddsNS<-estimateDispersions(ddsNS)
ddsNS<-nbinomWaldTest(ddsNS)
resNS<-results(ddsNS)
deseq2NS<-as.data.frame(apply(resNS,2,as.vector))
rownames(deseq2NS)<-rownames(resNS)
deseq2NS$sig<-ifelse(!is.na(deseq2NS$padj)&deseq2NS$padj<fdr&abs(deseq2NS$log2FoldChange)>logFC,1,0)
DEdeseq2NS<-subset(deseq2NS,sig==1)
# output the result without SpikeIn
deseq2NS<-cbind(df.part1,deseq2NS)
# add the normalized counts
deseq2NS<-cbind(deseq2NS,counts(ddsNS,normalized=T,replaced=F))
write.table(deseq2NS,paste0(outName,"_DiffPeaks_withoutSpikeIn_result.tsv"),sep="\t",row.names=F,quote=F)



# 3. QC statistics
print(paste0(nrow(DEdeseq2)," spike-in DE peaks"))
print(paste0(nrow(DEdeseq2NS)," without spike-in DE peaks"))
DEoverlap<-length(intersect(rownames(DEdeseq2),rownames(DEdeseq2NS)))
print(paste0(DEoverlap," shared DE results between the two methods"))
if(DEoverlap>0) {
	alist=list(spikeInNormalized=rownames(DEdeseq2),noSpikeInNormalized=rownames(DEdeseq2NS))
	alist=alist[order(sapply(alist,length),decreasing=T)]
	Plot_Venn(alist=alist,
	filename="venn.jpeg",pal="default",CategoryNames=T)
}

pdata<-CalcPCA(counts(dds,normalized=T,replaced=F),top_ratio=0.1,logTransform=T)
pdata<-PlotPCA(pdata,group,top_ratio=0.1) + ggtitle("spike-in normalized")
pdataNS<-CalcPCA(counts(ddsNS,normalized=T,replaced=F),top_ratio=0.1,logTransform=T)
pdataNS<-PlotPCA(pdataNS,group,top_ratio=0.1) + ggtitle("without spike-in normalized")
pdataTotal<-arrangeGrob(pdata,pdataNS,nrow=1,ncol=2)
ggsave("PCA.jpeg",pdataTotal,width=24,height=8,units="cm",dpi=600)

jpeg("hist_pval.jpeg",width=16,height=8,units="cm",res=600)
par(mfrow=c(1,2))
hist(deseq2$pvalue,xlab="raw p-value",main="spike-in normalized")
hist(deseq2NS$pvalue,xlab="raw p-value",main="without spike-in normalized")
dev.off()




save.image(paste0(outName,".RData"))
print("done")

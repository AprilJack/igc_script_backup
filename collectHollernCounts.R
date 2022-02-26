#!/gpfs/tools/R-4.0.3/bin/Rscript 
## R-4.0.3 this version has all packages installed
# Ling Huang Salk IGC 07-21-2021
# to collect the reads from hollern lab STAR-salmon pipeline
# it will iteratively go through all sub-folders and collect viral counts and salmon results, it should also convert salmon transcripts counts to gene level counts by tximport package suggested by salmon tutorial.
# regarding the viral reads, we would need to collect them from the custom java scripts instead of salmon. This is the new updates after checking the nextflow script from Dan.

args<- commandArgs(trailingOnly = TRUE)
if (length(args)!=3){
  cat("Usage: collectHollernCounts.R [working dir] [star dir name] [virus dir name]\n")
  cat("Example: collectHollernCounts.R ./ mm10.hollern.star virus.hollern.star \n")
  stop()
}

wdir=args[1]
starName=args[2]
virusName=args[3]

# 1. checking all the files needed.
cat("reading in summary.csv in the working directory \n")
summary.file<-paste0(wdir,"/summary.csv")
if(file.exists(summary.file)) {
  dfsummary<-read.table(summary.file,sep=",",header=T)
} else {
  stop("there is no summary.csv in the working directory!")
}
nSamples<-ncol(dfsummary)-1
# remove the Log in colnames
colnames(dfsummary)<-gsub("Log","",colnames(dfsummary))



cat("checking virus.gtf in the virus reference folder \n")
gtf.file<-paste0("/gpfs/genomes/",virusName,"/virus.gtf")
if(!file.exists(gtf.file)) {
  stop(paste0(gtf.file," does not exist!\n"))
}
virusgtf<-read.table(gtf.file,sep="\t",header=F,quote="",comment.char="")


cat("checking gene id and symbol info in the STAR reference folder \n")

tx.file<-paste0("/gpfs/genomes/",starName,"/tx2gene.tsv")
if(!file.exists(tx.file)) {
  stop(paste0(tx.file," does not exist!\n"))
}
symbol.file<-paste0("/gpfs/genomes/",starName,"/gene2symbol.tsv")
if(!file.exists(symbol.file)) {
  stop(paste0(symbol.file," does not exist!\n"))
}

library(tximport)
library(readr)
tx2gene <- read_delim(tx.file,delim="\t")
gene2symbol <- read.table(symbol.file,sep="\t",header=T)



# 2. collect viral and salmon output files in the wdir
cat("searching for salmon output quant.sf \n")
allsalmon.file<-list.files(wdir,pattern="quant.sf",recursive = T,full.names = T)
salmon.file<-allsalmon.file[grepl("_salmon",allsalmon.file)]

cat("searching for viral_read_counts.txt \n")
virus.file<-list.files(wdir,pattern="viral_read_counts.txt",recursive = T,full.names = T)

if(length(virus.file)<1) {
  stop("no viral_read_counts.txt found in the working directory!\n")
}
if(length(salmon.file)<1) {
  stop("no _salmon/quant.sf found in the working directory!\n")
}

# check if the salmon.file and virus.file are matched...ignore for now.
#if(length(virus.file)!=length(salmon.file)) {
#  stop("viral_read_counts.txt number and _salmon/quant.sf number is not the same!\n")
#}
# test if the files names are the same
#viralSamples<-sapply(virus.file,function(x) strsplit(x,"_virus",fixed=T)[[1]][1])
#test.salmon.file<-paste0(viralSamples,"_salmon/quant.sf")
#if(!identical(test.salmon.file,salmon.file)) {
#  stop("viralSamples and _salmon/quant.sf are not matched!\n")
#}

# 3. summarizing viral counts
cat("summarizing viral counts \n")
cat("matching viral sample names with star sample names \n")
dfvirus.meta<-data.frame(directory=virus.file)
dfvirus.meta$fileName<-sapply(dfvirus.meta$directory,function(x) rev(strsplit(x,"/",fixed=T)[[1]])[1])
dfvirus.meta$sampleName<-gsub("_viral_read_counts.txt","",dfvirus.meta$fileName)
dfvirus.meta$compatibleName<-make.names(dfvirus.meta$sampleName)
# make sure that the compatible names are found in summary file
if(!all(dfvirus.meta$compatibleName %in% colnames(dfsummary))){
  cat("some virus sample names are not found in the STAR output!\n")
  cat("not found samples: \n")
  cat(paste0(setdiff(dfvirus.meta$compatibleName,colnames(dfsummary)),collapse = "\n"))
  cat("\n")
  stop()
}

if(!all(colnames(dfsummary)[-1] %in%  dfvirus.meta$compatibleName)){
  cat("some star sample names are not found in the viral samples!\n")
  cat("not found samples: \n")
  cat(paste0(setdiff(colnames(dfsummary),dfvirus.meta$compatibleName),collapse = "\n"))
  cat("\n")
  stop()
}

virus.list<-lapply(virus.file,function(x) read.table(x,sep="\t",header=F))
virus.counts<-do.call(rbind,virus.list)
rownames(virus.counts)<-dfvirus.meta$compatibleName
virus.counts<-virus.counts[,-1]
virus.counts<-as.data.frame(t(virus.counts))
virus.counts$VirusName<-virusgtf$V1
virus.counts<-virus.counts[,c("VirusName",dfvirus.meta$compatibleName)]
cat("writing viral count output to viral_counts.tsv \n")
write.table(virus.counts,"viral_counts.tsv",sep="\t",row.names=F)


# 3. summarizing salmon counts
cat("summarizing salmon counts \n")
cat("matching salmon sample names with star sample names \n")
salmon.list<-lapply(salmon.file,function(x) read.table(x,sep="\t",header=T,quote="",comment.char=""))
dfsalmon.meta<-data.frame(directory=salmon.file)
dfsalmon.meta$sampleName<-gsub("/quant.sf","",dfsalmon.meta$directory,fixed=T)
dfsalmon.meta$sampleName<-sapply(dfsalmon.meta$sampleName,function(x) rev(strsplit(x,"/",fixed=T)[[1]])[1])
dfsalmon.meta$sampleName<-gsub("_salmon","",dfsalmon.meta$sampleName)
dfsalmon.meta$compatibleName<-make.names(dfsalmon.meta$sampleName)
if(!all(dfsalmon.meta$compatibleName %in% colnames(dfsummary))){
  cat("some salmon sample names are not found in the STAR output!\n")
  cat("not found samples: \n")
  cat(paste0(setdiff(dfsalmon.meta$compatibleName,colnames(dfsummary)),collapse = "\n"))
  cat("\n")
  stop()
}

if(!all(colnames(dfsummary)[-1] %in% dfsalmon.meta$compatibleName )){
  cat("some star sample names are not found in the salmon samples!\n")
  cat("not found samples: \n")
  cat(paste0(setdiff(colnames(dfsummary),dfsalmon.meta$compatibleName),collapse = "\n"))
  cat("\n")
  stop()
}

cat("summarizing the salmon count table by tximport package \n")
names(salmon.file)<-dfsalmon.meta$compatibleName
txi <- tximport(salmon.file, type = "salmon", tx2gene = tx2gene)
# names(txi)
# head(txi$counts) 
## countsFromAbundance are just flags indicating whether the counts are scaled (default: no)
dfsalmon.counts<-data.frame(GeneID=rownames(txi$counts))
dfsalmon.counts$Symbol<-gene2symbol$gene_name[match(dfsalmon.counts$GeneID,gene2symbol$gene_id)]
dfsalmon.tpm<-dfsalmon.counts
dfsalmon.counts<-cbind(dfsalmon.counts,txi$counts)
dfsalmon.tpm<-cbind(dfsalmon.tpm,txi$abundance)

cat("salmon counts output to salmon_counts.tsv \n")
cat("salmon tpm output to salmon_tpm.tsv \n")
write.table(dfsalmon.counts,"salmon_counts.tsv",sep="\t",row.names = F)
write.table(dfsalmon.tpm,"salmon_tpm.tsv",sep="\t",row.names=F)



# 4. collect viral count and salmon count info into the summary table
# NOTE: the sample order in dfsummary might not be the same as salmon.file
dfvirus.meta$viralReads<-colSums(virus.counts[,-1])
dfsalmon.meta$salmonReads<-colSums(dfsalmon.counts[,c(-1,-2)])

dfsummary<-as.data.frame(t(dfsummary))
dfsummary<-dfsummary[-1,]
colnames(dfsummary)<-c("UniqueRatio","MultipleRatio","UnmappedRatio","TotalReads")
dfsummary$salmonReads<-dfsalmon.meta$salmonReads[match(rownames(dfsummary),dfsalmon.meta$compatibleName)]
dfsummary$salmonRatio<-round(dfsummary$salmonReads/as.numeric(dfsummary$TotalReads)*100,2)
# NOTE: salmon deals with multi-mappers, so salmonRatio should be higher than uniqueRatio
dfsummary$viralReads<-dfvirus.meta$viralReads[match(rownames(dfsummary),dfvirus.meta$compatibleName)]

cat("save summary file to summary_salmon.tsv \n")
write.table(dfsummary,"summary_salmon.tsv",sep="\t",col.names=NA)

save.image("collectHollernCounts.RData")

cat("done \n")


############# 
# below is the code used to generate the reference files needed by tximport. Those files will need to be transferred to the star reference folder to run this script. 
# They are the tx2gene and gene2symbol files

# we first need to build the tx2gene table to contain mapping between transcript and gene
# they only need to be build once
#mm10db<-GenomicFeatures::makeTxDbFromGFF(file="/gpfs/genomes/mm10.hollern/gencode.vM27.chr_patch_hapl_scaff.annotation.gtf",format = "gtf",dataSource = "gencode.vM27.all",organism="Mus musculus")
#k <- AnnotationDbi::keys(mm10db, keytype = "TXNAME")
#tx2gene.mm10 <- AnnotationDbi::select(mm10db, k, "GENEID", "TXNAME")
#write.table(tx2gene.mm10,"/gpfs/genomes/mm10.hollern.star/tx2gene.tsv",sep="\t",row.names=F)
# extract gene symbol info as well?
#mm10.gtf <- rtracklayer::import("/gpfs/genomes/mm10.hollern/gencode.vM27.chr_patch_hapl_scaff.annotation.gtf")
#mm10.genes <- mm10.gtf[mm10.gtf$type == "gene"] # 55416, correct
#write.table(mm10.genes[,c("gene_id","gene_type","gene_name")],"/gpfs/genomes/mm10.hollern.star/gene2symbol.tsv",sep="\t",row.names=F)

#hg38db<-GenomicFeatures::makeTxDbFromGFF(file="/gpfs/genomes/hg38.hollern/gencode.v38.chr_patch_hapl_scaff.annotation.gtf",format = "gtf",dataSource = "gencode.v38.all",organism="Homo sapiens")
#khs <- AnnotationDbi::keys(hg38db, keytype = "TXNAME")
#tx2gene.hg38 <- AnnotationDbi::select(hg38db, khs, "GENEID", "TXNAME")
#write.table(tx2gene.hg38,"/gpfs/genomes/hg38.hollern.star/tx2gene.tsv",sep="\t",row.names=F)
#hg38.gtf <- rtracklayer::import("/gpfs/genomes/hg38.hollern/gencode.v38.chr_patch_hapl_scaff.annotation.gtf")
#hg38.genes <- hg38.gtf[hg38.gtf$type == "gene"] # 67049, correct
#write.table(hg38.genes[,c("gene_id","gene_type","gene_name")],"/gpfs/genomes/hg38.hollern.star/gene2symbol.tsv",sep="\t",row.names=F)





#!/usr/bin/env Rscript 
# calculate the norm_factor.txt for each tag directory based on the unique spike-in mapped reads.
# update: 09-01-2020 add support for flag --force to force the match based on input orders. 
# for yfarsakoglu we can use --spikeInTags tags/*_spike --targetTags tags/*[0-9]
options(stringsAsFactors=F)
args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 4) {
	print("command line arguments:")
	print("1. spike-in tag directories")
	print("2. tag directories that need to be normalized")
	print("Usage example: CalcNormFactor.R --spikeInTags spike-in_tags/* --targetTags tagDir1 tagDir2 ...")
	print("The normalization factor will be calculated by based on the ratio of reads mapped to spike-in and target reference. For bigWig visualization scaling.")
	print("use flag --spikeInTags and --targetTags to input the tag directories, also support multiple tag directories and * in bash style.")
	print("NOTE: in order to match the spike-in directories with the target tag directories, the program assume the two directories have the same name but in different path. So it doesn't work if they are named differently. To bypass this behavior, use flag --force so the samples will be matched based on the order of input.")
	stop()
}
# args<-c("--spikeInTags","/home/styagi/cutrundata/rawdata/yalignments/tags/*","--targetTags","/home/styagi/cutrundata/rawdata/tags/*")

spikeInIndex<-which(args=="--spikeInTags")
targetIndex<-which(args=="--targetTags")
forceIndex<-which(args=="--force")
matchOrder<-FALSE
if(length(forceIndex)==1) {
	# turn on the matchOrder mode and remove it from args
	matchOrder<-TRUE
	args<-args[(0-forceIndex)]
}


GetDirectory<-function(oneDir){
	# first check if it is a single directory or contains wildcard
	if(!grepl("*",oneDir,fixed=T)) {
		if(file.exists(oneDir)&file.exists(paste0(oneDir,"/tagInfo.txt"))) {
			return(oneDir)
		} else {
			return(NA)
		}
	} 
	CMD=paste0("ls ",oneDir) # it takes care of the regular expression part 
	dirs=system(CMD,intern=T)
	# remove the column from the returned output
	dirs=sapply(dirs,function(x) gsub(":","",x))
	dirs<-dirs[sapply(dirs,dir.exists)]
	dirsTagInfo<-sapply(dirs,function(x) file.exists(paste0(x,"/tagInfo.txt")))
	dirs<-dirs[dirsTagInfo]
	names(dirs)<-NULL
	return(dirs)	
}


SpikeIn<-vector()
Target<-vector()
if(spikeInIndex < targetIndex) {
	for(i in seq(spikeInIndex+1,targetIndex-1)) {
		SpikeIn<-c(GetDirectory(args[i]),SpikeIn)
	}
	for(i in seq(targetIndex+1,length(args))) {
		Target<-c(GetDirectory(args[i]),Target)
	}
} else {
	for(i in seq(targetIndex+1,spikeInIndex-1)) {
		Target<-c(GetDirectory(args[i]),Target)
	}
	for(i in seq(spikeInIndex+1,length(args))) {
		SpikeIn<-c(GetDirectory(args[i]),SpikeIn)
	}
}

SpikeIn<-SpikeIn[!is.na(SpikeIn)]
Target<-Target[!is.na(Target)]

print("Spike-In tag directories:")
print(SpikeIn)
print("Target tag directories:")
print(Target)

if(length(SpikeIn)<1|length(Target)<1) {
	stop("no tag directories available")
}
if(length(SpikeIn)!=length(Target)) {
	stop("number of tag directories differs in spike-in and target samples")
}

#library(edgeR) it is not used in this script

# first, get the spike-in expression matrix
SpikeInData<-lapply(SpikeIn,function(x) read.table(paste0(x,"/tagInfo.txt"),sep="\t",header=T,quote="",comment.char=""))
# the problem is that we cannot assume that every tagDirectory mapped to the spike-in reference has expression in all chromosomes.
spikeInChrom<-unlist(lapply(SpikeInData,nrow))
if(!all(spikeInChrom==max(spikeInChrom))) {
	print("Warning! some spike-in tag directories do not contain reads mapped to all chromosomes")	
}
spikeInChromMax<-which(spikeInChrom==max(spikeInChrom))[1] # get the first index
# NOTE: new version of Homer produces 1 more lines here...
test1<-SpikeInData[[1]]
GCindex<-which(grepl("averageFragmentGCcontent",test1$name))
dfspikein<-data.frame(matrix(ncol=length(SpikeIn),nrow=(max(spikeInChrom)-GCindex-1)))
colnames(dfspikein)<-make.names(SpikeIn)
rownames(dfspikein)<-SpikeInData[[spikeInChromMax]][(GCindex+1):(nrow(SpikeInData[[spikeInChromMax]])-1),1]
for(i in seq(1,length(SpikeIn))) {
	tempData<-SpikeInData[[i]][9:(nrow(SpikeInData[[i]])-1),c(1,3)]
	dfspikein[,i]<-tempData$Total.Tags[match(rownames(dfspikein),tempData$name)]
}
# replace the NA with zeros
dfspikein[is.na(dfspikein)]<-0

# prepare the final output data frame
dfnorm<-data.frame(SpikeInSample=SpikeIn,LibSizeSpikeIn=colSums(dfspikein))
# we need to pair the Target and SpikeIn
dfnorm$sample<-sapply(dfnorm$SpikeInSample,function(x) rev(strsplit(x,"/",fixed=T)[[1]])[1])
dftarget<-data.frame(targetDir=Target)
dftarget$sample<-sapply(dftarget$targetDir,function(x) rev(strsplit(x,"/",fixed=T)[[1]])[1])
# add support to match just based on order
if(matchOrder) {
	dfnorm$TargetSample<-dftarget$targetDir
} else {
	dfnorm$TargetSample<-dftarget$targetDir[match(dfnorm$sample,dftarget$sample)]
}

if((nrow(dfnorm)==1&sum(is.na(dfnorm$TargetSample))>0)|sum(is.na(dfnorm$TargetSample))>1) {
	print("Error: The name of the tag directories of spike-in and target do not match with each other")
	print("spike-in directories:")
	print(paste(dfnorm$sample[is.na(dfnorm$TargetSample)],sep="\n"))
	print("target directories:")
	print(paste(dftarget$sample[is.na(dfnorm$TargetSample)],sep="\n"))
	stop("please rename the spike-in tag directories using the same name as their corresponding target directories")
}

# we also need the total mapped human reads for lib.size
TargetData<-lapply(dfnorm$TargetSample,function(x) read.table(paste0(x,"/tagInfo.txt"),sep="\t",header=T,quote="",comment.char=""))
dfnorm$LibSizeTotal<-unlist(lapply(TargetData, function(x) x[1,3]))

# we need to check if there are too few reads mapped to the spike-in directories, we need to throw out an Error
if(all(dfnorm$LibSizeSpikeIn<1000)) { # 1000 is an arbitaray number...we might want even more reads.
	print("Warning! you may have very low coverage to the spike-in reference. It might be dangerous to use the low reads as spike-in.")
}

dfnorm$ratio<-dfnorm$LibSizeSpikeIn/dfnorm$LibSizeTotal # use the total number of mapped reads to spike-in to calculate the ratio
# also test the ratio, if the difference between ratios > 10, caution!
if(max(dfnorm$ratio) > 10 * min(dfnorm$ratio)) {
	print("Warning! The spike-in ratios are very different between your samples. It might be caused by unexpected high/low spike-in fragments in some of the samples. It might be dangerous to proceed using the calculated ratios. Needs manual examination!")
}

ratioExpected<-0.001 # this is an arbitary number, to get the scaling working even for multiple batches. It does not matter as long as the relatvie difference between samples are kept.
dfnorm$BigWigFactor<-1000000 * ratioExpected/dfnorm$ratio

# print a warning message if overwritting the output 
if(file.exists("spikeIn_norm.tsv") | file.exists("spikeIn_raw.tsv")) {
	stop("Warning! found spikeIn_norm.tsv or spikeIn_raw.tsv already exists in the current working directory. Please delete them if they are no longer needed and re-run this command.")
}

write.table(dfnorm,"spikeIn_norm.tsv",sep="\t",row.names=F,quote=F)
print("Output is written in spikeIn_norm.tsv")

# we need the expression table for each different comparison
# change the colnames of spikeIn_raw from spike-in to target
if(identical(rownames(dfnorm),colnames(dfspikein))) {
	colnames(dfspikein)<-make.names(dfnorm$TargetSample)
}



write.table(dfspikein,"spikeIn_raw.tsv",sep="\t",quote=F)
print("spikeIn matrix is written in spikeIn_raw.tsv")

#jpeg("scatter_ratio_NormFactor.jpeg",width=12,height=12,units="cm",res=600)
#plot(dfnorm$ratio,dfnorm$NormFactor)
#dev.off()


for(i in seq(1,nrow(dfnorm))) {
	outfile<-paste0(dfnorm$TargetSample[i],"/norm_factor.txt")
	normFactor<-dfnorm$BigWigFactor[i]
	write.table(normFactor,outfile,sep="\t",quote=F,row.names=F,col.names=F)
}
print("BigWig norm factor is written in norm_factor.txt under each tag directory")







print("done")

if (!require("gplots")) {
  install.packages("gplots", repos="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", repos ="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(RColorBrewer)
}

if (!require("DESeq2")) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("DESeq2")
}
ags<- commandArgs(trailingOnly = TRUE)
ags<- c("raw_no145.txt","design_matrix.txt")
if (length(ags)==0){
  cat("Takes a HOMER raw data matrix and a design matrix to ask specific questions of the data using GLM fits and contrasts.\n")
  cat("Usage: GLM HomerRAWTable design.txt\n")
  ags=c("raw.txt","design.txt")
}
data <- as.matrix(read.delim(file = ags[1],header = TRUE, sep = "\t"))
datafirstColumn = dim(data)[2]
for(i in dim(data)[2]:1)
{
  if(is.na(suppressWarnings(as.integer(data[2,i]))))
  {
    firstColumn = i+1
    break
  }
}
data2 = data[,firstColumn:dim(data)[2]]
rownames(data2) <- data[,firstColumn-1]
colnames(data2) <- colnames(data[,firstColumn:dim(data)[2]])
#Now let's go through and ask labels for each column
design=as.matrix(read.delim(file = ags[2],header = TRUE, sep = "\t"))
rownames(design)=design[,1]
design2=design[,2:length(colnames(design))]
design=design2
summary(design,maxsum=10)
model= readline("Please input a GLM model: e.g. ~Treatment+Batch:")
cat(sprintf("Model: %s\n",model))
dds<-DESeqDataSetFromMatrix(countData = data2, colData = design, design=as.formula(model))
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds = DESeq(dds)
while(TRUE) {
  query=readline(sprintf("Type in a(%s).Enter \"exit\" when done:\n",paste(resultsNames(dds),collapse=" ")))
  if(query == "exit")
    break
}
cat("Finished GLM. If you like the results send beer to ...")

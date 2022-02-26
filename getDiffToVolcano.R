if (!require("ggplot2")) {
  install.packages("ggplot2", repos="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(gplots)
}
if (!require("ggfortify")) {
  install.packages("ggfortify", repos="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(ggfortify)
}
#ags<- commandArgs(trailingOnly = TRUE)
ags= c("raw.txt","sampleSheet.txt")
if (length(ags)==0){
  cat("Usage: getDiffR raw.txt sampleSheet.txt design\n
      e.g. getDiffR raw.txt sampleSheet.txt ~Condition+Batch\n")
  ags=c("raw.txt","sampleSheet.txt")
}
cat(sprintf("Running DESeq2 on %s using samplesheet %s and design %s\n",ags[1],ags[2],ags[3]))
data <- as.matrix(read.delim(file = ags[1],header = TRUE, sep = "\t"))
if(dim(data)[2]>8 & !is.numeric(data[1,8]))
{
  cat("Assuming HOMER formated raw count table\n")
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
  dim(data2)

} else {
  cat("Assuming just one column for ids and one row for sample names\n")
  data2 = data[,2:dim(data)[2]]
  data2= matrix(as.numeric(data2),nrow=dim(data2)[1],ncol=dim(data2)[2])
  rownames(data2) <- data[,1]
  colnames(data2) <- colnames(data[,2:dim(data)[2]])
  dim(data2)
}
sample=read.delim(file = ags[2],header = TRUE, sep="\t")

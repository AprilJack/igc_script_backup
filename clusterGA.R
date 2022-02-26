# Clustering stats files
if (!require("gplots")) {
  install.packages("gplots", repos="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", repos ="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(RColorBrewer)
}
cat("Will cluster homer annotatePeaks -annStats files by genomic locatin enrichment\n")
ags<- commandArgs(trailingOnly = TRUE)
if (length(ags)==0){
  ags = list.files(path = ".", pattern="*_stats.txt")
  cat(sprintf("Loading %d files that end with _stats.txt by default.\n",length(ags)))  
}
data= array(dim = c(33,length(ags)))
counts = array(dim = c(33,length(ags)))
colnames(data)<-ags
for(i in 1:length(ags))
{
  x<-(read.delim(ags[i],sep="\t", header=TRUE,quote=""))
  data[,i]=as.numeric(data.matrix(x[14:46,4]))
  counts[,i]=as.numeric(data.matrix(x[14:46,2]))
  rownames(data)<-x[14:46,1]
}
data[is.na(data)]<-0
counts[is.na(counts)]<-0
data2<-matrix(data,nrow = 33,ncol = length(ags))
data2<-data2[rowSums(counts)>1,colSums(counts)>50]
colnames(data2)<-colnames(data)[colSums(counts)>50]
rownames(data2)<-rownames(data)[rowSums(counts)>1]

write.table(data2,file = "genomic_clustering.txt",sep = "\t", eol = "\r\n",col.names=NA)
my_palette <- colorRampPalette(c("white","white","white", "red", "black"))#colorRampPalette(c("white","red","red","red","red","red","red","red","red","red"))(n=599)
png("genomic_clusters.png", height=8*500, width=8*600,res=600,pointsize=8)
heatmap.2(data2, main = "Genomic Enrichment", trace="none", vline=, col=my_palette,margins=c(12,7),hclustfun = function(x) hclust(x,method="ward.D2"),keysize = 1,key.title=NA, colsep=seq(1:dim(data2)[2]),sepcolor = "gray")
dev.off()

# Clustering processFindGO terms
if (!require("gplots")) {
  install.packages("gplots", repos="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", repos ="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(RColorBrewer)
}

if (!require("wordcloud")) {
  install.packages("wordcloud", repos ="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(wordcloud)
}
ags<- commandArgs(trailingOnly = TRUE)
if (length(ags)==0){
  cat("Using enrichment.txt as default input file\n")
  ags=c("enrichment.txt")
}

data <- read.delim(file=ags[1],header = TRUE, sep = "\t")
mat=array(dim = c(dim(data)[1],dim(data)[2]/5))
colnames(mat)<-colnames(data[,seq(2,dim(data)[2],5)])
rownames(mat)<-data[order(data[,1]),1]
count=1
for( i in seq(1,dim(data)[2],5))
{
  mat[,count]=-log(data[order(data[,i]),i+1])
  count=count+1
}
rowMaxes = array(c(dim(data)[1]))
for(i in 1:dim(mat)[1])
{
  rowMaxes[i] = max(mat[i,])
}
data2<-t(mat[order(rowSums(mat),decreasing = TRUE),apply(mat, 2, sd)>0])
data2<-data2[apply(data2,1,sd)>0,]
write.table(t(data2),file = "enrichment_clustering.txt",sep = "\t", eol = "\r\n",col.names=NA)
my_palette <- colorRampPalette(c("white","red","red","red","red","red","red","red","red","red"))(n=599)
png("enrichment_cluster.png", height=max(1000,30*dim(data2)[1]), width=30*dim(data2)[2],res=300,pointsize=8)
heatmap.2(data2, main = "Enrichment Clustering",dendrogram = 'column', trace="none",Rowv = FALSE, col=my_palette,margins=c(10,10),sepwidth =c(0.01,0.01), colsep = seq(1,dim(data2)[2]), sepcolor = "gray",cexRow = 0.4,cexCol=0.4,distfun=function(c) as.dist(1-cor(t(c))), hclustfun = function(x) hclust(x,method="ward.D2"))

dev.off()

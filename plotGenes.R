if (!require("gplots")) {
  install.packages("gplots", repos="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", repos ="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(RColorBrewer)
}
ags<- commandArgs(trailingOnly = TRUE)
if (length(ags)<2){
  stop("Usage: plotGenes fpkm_file Genes.txt\n")
}

data <- as.matrix(read.delim(file = ags[1],header = TRUE, sep = "\t"))
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
genes <- as.matrix(read.delim(file=ags[2],header = FALSE))
data2<-data2[rownames(data2) %in% genes,]
expression= log(data2[apply(data2, 1, function(x) var(x)>0),]+1,base=2)
my_palette <- colorRampPalette(c("green","black","red"))(n=299)
my_palette2 <- colorRampPalette(c("white","red"))(n=299)
if(dim(data2)[2] > 30)
{
  my_palette <- colorRampPalette(c("green","green","black","red","red"))(n=299)
}
write.table(data2,file = sprintf("%s_LNFPKM_%s.txt",ags[1],ags[2]),sep = "\t", eol = "\r\n",col.names=NA)
png(sprintf("%s_%s.png",ags[1],ags[2]), width=max(50*dim(data2)[2],1200), height=max(12*dim(data2)[1],1200),res=300,pointsize=8)
heatmap.2(expression, main = "Ln(FPKM+1)",trace="none", col=my_palette,margins = c(12,12),cexRow = 0.25,cexCol=0.5,distfun=function(c) as.dist(1-cor(t(c))),scale="row", hclustfun = function(x) hclust(x,method="ward.D2"),key.title = "")
dev.off()
png(sprintf("%s_%s2.png",ags[1],ags[2]), width=max(50*dim(data2)[2],1200), height=max(12*dim(data2)[1],1200),res=300,pointsize=8)
heatmap.2(expression, main = "Ln(FPKM+1) not scaled",trace="none", col=my_palette2,margins = c(12,12),cexRow = 0.25,cexCol=0.5,distfun=function(c) as.dist(1-cor(t(c))),scale="none", hclustfun = function(x) hclust(x,method="ward.D2"),key.title = "")
dev.off()




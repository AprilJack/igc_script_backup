####
#  Avg. Log2TPM plot R Script. 
#  Set the input file, Gene to Highlight, Highlight Color, and image resolution etc. below.
#  Assumes your input has one column for gene names and the rest are Avg. Log2TPM counts for each condition. 
#
#  Written by Maxim N. Shokhirev, (C) Salk IGC 2017
####
setwd("/Volumes/viper/analyses/maxsh/Hsu/171204_rnaseq")
data<-read.delim(file="avgLog2TPM.txt",header=TRUE)
data2=data[,2:dim(data)[2]]
rownames(data2)<-data[,1]
data3=data2[rowMeans(data2)>1.5,]
cat(sprintf("Filtered genes with average Avg. Log2TPM less than 1.5 leaving %d of %d",dim(data3)[1],dim(data2)[1]))

######## GENE HIGHLIGHTING ############
GeneToHighlight="B4GALNT1"
HighlightColor="green2"

############ Generating a bit matric of all pairwise comparisons ###########
png(file = "AllvsAll.png",width = 4000,height=4000,res = 500)
par(mfrow=c(dim(data2)[2],dim(data2)[2]))
par(mar = c(2, 2, 0.5, 0.5))
par(oma = c(3, 3.5, 0, 0))
for(i in 1:dim(data2)[2])
{
  for(j in 1:dim(data2)[2])
  {
    data3$Color=rgb(red = 0,green=0,blue=0,alpha=0.1)
    plot(data3[,i],data3[,j], pch=16,col=data3$Color,cex=0.5)
    dataUp= data3[data3[,i]-data3[,j]>1.5,]
    dataDown= data3[data3[,i]-data3[,j]< -1.5,]
    geneData= data3[rownames(data3) == GeneToHighlight,]
    write.csv(file=sprintf("%svs%s.csv",colnames(data3)[i],colnames(data3)[j]),x=rbind(dataUp,dataDown,geneData),row.names = TRUE,quote = FALSE)
    points(dataUp[,i],dataUp[,j], pch=16,col="red",cex=0.6)
    points(dataDown[,i],dataDown[,j], pch=16,col="blue",cex=0.6)
    points(geneData[,i],geneData[,j], pch=16, col=HighlightColor,cex=0.7)
    text(geneData[,i],geneData[,j]+1,sprintf("%s",GeneToHighlight),col=HighlightColor,cex=0.3)
    text(2,10,sprintf("%d",sum(data3[,i]-data3[,j]< -1.5)),col="blue")
    text(10,2,sprintf("%d",sum(data3[,i]-data3[,j]> 1.5)),col="red")
    cat(sprintf("%d %d\n",i,j))
  }
}
################ LABELS #####################
mtext(paste(colnames(data2),collapse=","), side = 1, outer = TRUE, line = 2)
mtext(paste(rev(colnames(data2)),collapse=","), side = 2, outer = TRUE, line = 2)
dev.off()

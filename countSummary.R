ags<- commandArgs(trailingOnly = TRUE)
if (length(ags)==0){
  cat("Usage: countSummary [summary.csv file]\n")
  ags=c("summary.csv")
}
data=read.csv(ags[1], header = TRUE)
data2=data[,2:dim(data)[2]]
data2[5,]=data2[5,]/1000000.0
data2[6,]=data2[5,]*data2[1,]/100.0
data2[7,]=data2[4,]/1000000.0
png(sprintf("%s.png",substr(ags[1],1,nchar(ags[1])-4)), width=1200, height=1200,res=300,pointsize=8)
par(mar=c(10,7,2,1))
barplot(as.matrix(data2[5:7,]),names.arg=colnames(data2),beside=TRUE, main="Read Count Summary",ylab="Read # (millions)",cex.names=0.5,las=2)
legend("topright", c("Total","Uniquely Mapping","In exons"), cex=0.7, pt.cex = 1,fill = c("black","gray","white"))
dev.off();
cat("Finished\n")

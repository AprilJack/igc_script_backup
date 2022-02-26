if (!require("gplots")) {
  install.packages("gplots", repos="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", repos ="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(RColorBrewer)
}
if (!require("Rtsne")) {
  install.packages("Rtsne", repos ="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(Rtsne)
}
if (!require("minimist")) {
  install.packages("minimist", repos ="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(minimist)
}

if (!require("SLICER")) {
  library("devtools")
  install_github("jw156605/SLICER")
  library(SLICER)
}

#ags=minimist(c("raw_5-9-17_justAge.txt","-c","decadesC.txt","--keepSamples"))
ags=minimist(commandArgs(trailingOnly = TRUE))
if (length(ags$'_')==0){
  cat("Usage: plotTSNE raw.txt [options] \n
      -c\tclusters.txt:\tspecify that you want to color the samples manually (order should be identical)
      -g\tgenes.txt:\tspecify that you want to color the samples by the avg. expression of genes specified
      --keepSamples:\tdo not filter samples based on their expression (if you are plotting bulk data)
      --keepGenes:\tdo not filter genes based on their normalized quantile expression")
  ags=c("raw.txt")
}
cat("Loading in the data...\n")
data <- as.matrix(read.delim(file = ags$`_`,header = TRUE, sep = "\t"))
cat("Removing duplicates and normalizing the raw counts...\n")
# remove duplicates
#data <-subset(data,!duplicated(data[,1]))
# remove row headers
if(dim(data)[2]>8 & !is.numeric(as.numeric(data[1,8])))
{
  for( i in 1:dim(data)[1]) 
  {
    data[i,8]= strsplit(data[i,8],"[|]")[[1]][1]
    if( i %% 1000==0)
      cat(".")
  }
  data2 <- data[,9:dim(data)[2]]
  data2= matrix(as.numeric(data2),nrow=dim(data2)[1],ncol=dim(data2)[2])
  rownames(data2)<-data[,8]
  colnames(data2)<-colnames(data)[9:length(colnames(data))]
} else {
  data2 <- data[,2:dim(data)[2]]
  data2= matrix(as.numeric(data2),nrow=dim(data2)[1],ncol=dim(data2)[2])
  rownames(data2)<-data[,1]
  colnames(data2)<-colnames(data)[2:length(colnames(data))]
}

#data2= as.data.frame(data2)
# add them as rownames
startingCells=dim(data2)[2]
startingGenes=dim(data2)[1]

if(length(ags$keepSamples)==0) 
{
  quantileThreshold=quantile(colSums(data2),probs=0.3)
  data3 <- data2[,colSums(data2)>quantileThreshold]
} else {
  data3 = data2
  quantileThreshold=-1
}
goodCells=dim(data3)[2]
cat(sprintf("Filtering cells in %s. Went from %g to %g samples which had counts greater than %g\n",ags[1],startingCells,goodCells,quantileThreshold))
data3 = data3 %*% diag(1000000/colSums(data3))
data3 = normalize.quantiles(data3)
colnames(data3)=colnames(data2)[colSums(data2)>quantileThreshold]
rownames(data3)=rownames(data2)
if(length(ags$keepGenes)==0) 
{
  data4 <- log(data3[apply(data3, 1, function(x) (sum(x>5)>dim(data3)[2]/10)&var(x)>0),]+1,base=2)
} else {
  data4 = log(data3+1,base=2)
}
for(i in 1:length(colnames(data4)))
{
  colnames(data4)[i]=strtrim(colnames(data4)[i],10)
}
goodGenes=dim(data4)[1]
cat(sprintf("Filtering genes in %s. Went from %g to %g genes which had more than 5 normalized counts in %g cells\n",ags[1],startingGenes,goodGenes,dim(data3)[2]/10))
set.seed(42)
#Get mean expression for all genes of interest
if(length(ags$g)>1)
{
  genes <- as.matrix(read.delim(file=ags$g,header = FALSE))
  expression<-colMeans(data4[rownames(data4) %in% genes,])
  my_paletteG <- colorRampPalette(c("black","red"))(n=max(1,max(expression)))
}
#get cluster color labels
if(length(ags$c)>0)
{
  clusters <- as.matrix(read.delim(file=ags$c,header = FALSE))
  if(length(clusters)>dim(data4)[2])
    clusters=clusters[colSums(data2)>quantileThreshold]
}

tsne<-Rtsne(t(data4), dims = 2, perplexity = 10,
            theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 500,
            verbose = TRUE, is_distance = FALSE)
#Color by count
sums<-colSums(data4)

my_palette <- colorRampPalette(c("black","red"))(n=max(sums)-min(sums))

label="TSNE colored by normalized count"
if(length(ags$g)>0)
{
  label=sprintf("TSNE colored by %s",ags$g)
  if(length(ags$c)>0)
  {
    label=sprintf("TSNE colored by %s and %s",ags$g,ags$c)
  }
} else {
  if(length(ags$c)>0)
  {
    label=sprintf("TSNE colored by %s",ags$c)
  }
}
for(i in 1:dim(tsne$Y)[1])
{
  tsne$bg[i]=my_palette[sums[i]-min(sums)+1]
  colorData=sums
  if(length(ags$g)>0)
  {
    tsne$bg[i]=my_paletteG[max(1,expression[i])]
    colorData=expression
    if(length(ags$c)>0)
    {
      tsne$fg[i]=clusters[i]
    }
  } else {
    if(length(ags$c)>0)
    {
      tsne$bg[i]=clusters[i]
    }
  }
}


png(sprintf("%s_tSNE.png",basename(ags$'_')), width=3000, height=3000,res=300,pointsize=8)
plot(tsne$Y,pch=21,col=tsne$fg,bg=tsne$bg,main=label,xlab="tSNE 1",ylab="tSNE 2")
text(tsne$Y[,1],tsne$Y[,2],cex=0.5,labels=colnames(data4))
if(length(ags$c)==0)
  color.bar(data=colorData,col=c("black","red"),xmin = min(tsne$Y[,1]),xmax=max(tsne$Y[,1]),y=max(tsne$Y[,2])+2)
dev.off()


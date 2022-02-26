

if (!require("gplots")) {
  install.packages("gplots", repos="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(gplots,quietly = T)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", repos ="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(RColorBrewer,quietly = T)
}
if (!require("igraph")) {
  install.packages("igraph", repos="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(igraph,quietly = T)
}
if (!require("visNetwork")) {
  install.packages("visNetwork", repos ="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(visNetwork,quietly = T)
}
if (!require("pkgmaker")) {
  install.packages("pkgmaker", repos ="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(pkgmaker,quietly = T)
}
if (!require("rjson")) {
  install.packages("rjson", repos ="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(rjson,quietly = T)
}
if (!require("data.table")) {
  install.packages("data.table", repos ="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(data.table,quietly = T)
}
if (!require("PythonInR")) {
  install.packages("PythonInR", repos ="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(PythonInR,quietly = T)
}
if (!require("parallel")) {
  install.packages("parallel", repos ="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(parallel,quietly = T)
}
if (!require("doParallel")) {
  install.packages("doParallel", repos ="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(doParallel,quietly = T)
}
if (!require("foreach")) {
  install.packages("foreach", repos ="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(foreach,quietly = T)
}
if (!require("WebGestaltR")) {
  install.packages("WebGestaltR", repos ="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(WebGestaltR,quietly = T)
}
if (!require("RCurl")) {
  install.packages("RCurl", repos ="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(RCurl,quietly = T)
}
options(stringsAsFactors = F)
#setwd("/Volumes/viper/analyses/maxsh/IGC/testing")
#ags<-c("diffExp.txt")
ags<- commandArgs(trailingOnly = TRUE)
if (length(ags)==0){
  cat("Usage: Enrichment single-column/HOMERPeaks/HOMERDiffExpResults [top N pathways, default=15]\n
      Will identify the input type and send results to WebGestaltR for enrichment analysis.
      Then, it will create an interactive network graph showing the overlapping genes between terms
      Results will go to html and folder with same prefix as your input file. 
      Thanks to Ling Huang for the interactive network script.\n")
  quit(status=1, save="no")
} 
split= strsplit(ags[1],"[.]")[[1]]
base=paste(split[1:length(split)-1],collapse = ".")
data <- read.delim(file = ags[1],header = TRUE, sep = "\t")
if(dim(data)[2]>15 & !is.numeric(data[1,16]))
{
  cat("Loading in genes that are next to your peaks!\n")
  ## This is a peak file, so lets just grab the symbols from column 14
  genes<-unique(as.matrix(data[,16]))
} else if (dim(data)[2]> 8 && is.numeric(data[1,9]) && !is.numeric(data[1,8])) {
  ## Looks like an fpkm table or differential expression list
  cat("Loading in a differential expression table or fpkm table\n")
  data2=as.matrix(data[,8])
  for( i in 1:dim(data)[1]) 
  {
    data2[i,1]=strsplit(data2[i,1],"[|]")[[1]][1]
  }
  genes = data2
} else {
  ## assume just that the first column is a list of symbols
  cat("Loading a column of gene symbols!\n")
  genes<-data[,1]
}
topN=15
if(length(ags)>1)
{
  topN=as.numeric(ags[2])
}
# let's query webgestaltR
if(identical(toupper(genes)[1],genes[1]))
{
  enrichment=WebGestaltR(enrichMethod="ORA", organism="hsapiens", 
            enrichDatabase="geneontology_Biological_Process_noRedundant",enrichDatabaseFile=NULL, 
            enrichDatabaseType=NULL,enrichDatabaseDescriptionFile=NULL,interestGeneFile=NULL, 
            interestGene=as.vector(genes),interestGeneType="genesymbol",collapseMethod="mean",referenceGeneFile=NULL,
            referenceGene=NULL,referenceGeneType=NULL,referenceSet="genome_protein-coding", minNum=5, maxNum=500,
            fdrMethod="BH",sigMethod="top",fdrThr=0.05,topThr=topN,dNum=20,perNum=1000,
            lNum=20,is.output=TRUE,outputDirectory=getwd(),projectName=paste(c(base,"_GO"),collapse=""),keepGSEAFolder=FALSE,
            methodType="R",dagColor="continuous",hostName="http://www.webgestalt.org/")
  enrichment2=WebGestaltR(enrichMethod="ORA", organism="hsapiens", 
                         enrichDatabase="pathway_Reactome",enrichDatabaseFile=NULL, 
                         enrichDatabaseType=NULL,enrichDatabaseDescriptionFile=NULL,interestGeneFile=NULL, 
                         interestGene=as.vector(genes),interestGeneType="genesymbol",collapseMethod="mean",referenceGeneFile=NULL,
                         referenceGene=NULL,referenceGeneType=NULL,referenceSet="genome_protein-coding", minNum=5, maxNum=500,
                         fdrMethod="BH",sigMethod="top",fdrThr=0.05,topThr=topN,dNum=20,perNum=1000,
                         lNum=20,is.output=TRUE,outputDirectory=getwd(),projectName=paste(c(base,"_Reactome"),collapse=""),keepGSEAFolder=FALSE,
                         methodType="R",dagColor="continuous",hostName="http://www.webgestalt.org/")
  enrichment3=WebGestaltR(enrichMethod="ORA", organism="hsapiens", 
                          enrichDatabase="pathway_Wikipathway",enrichDatabaseFile=NULL, 
                          enrichDatabaseType=NULL,enrichDatabaseDescriptionFile=NULL,interestGeneFile=NULL, 
                          interestGene=as.vector(genes),interestGeneType="genesymbol",collapseMethod="mean",referenceGeneFile=NULL,
                          referenceGene=NULL,referenceGeneType=NULL,referenceSet="genome_protein-coding", minNum=5, maxNum=500,
                          fdrMethod="BH",sigMethod="top",fdrThr=0.05,topThr=topN,dNum=20,perNum=1000,
                          lNum=20,is.output=TRUE,outputDirectory=getwd(),projectName=paste(c(base,"_Wikipathways"),collapse=""),keepGSEAFolder=FALSE,
                          methodType="R",dagColor="continuous",hostName="http://www.webgestalt.org/")
  enrichment4=WebGestaltR(enrichMethod="ORA", organism="hsapiens", 
                          enrichDatabase="drug_DrugBank",enrichDatabaseFile=NULL, 
                          enrichDatabaseType=NULL,enrichDatabaseDescriptionFile=NULL,interestGeneFile=NULL, 
                          interestGene=as.vector(genes),interestGeneType="genesymbol",collapseMethod="mean",referenceGeneFile=NULL,
                          referenceGene=NULL,referenceGeneType=NULL,referenceSet="genome_protein-coding", minNum=5, maxNum=500,
                          fdrMethod="BH",sigMethod="top",fdrThr=0.05,topThr=topN,dNum=20,perNum=1000,
                          lNum=20,is.output=TRUE,outputDirectory=getwd(),projectName=paste(c(base,"_DrugBank"),collapse=""),keepGSEAFolder=FALSE,
                          methodType="R",dagColor="continuous",hostName="http://www.webgestalt.org/")
  enrichment=merge(enrichment,enrichment2,all=T)
  
}else {
  enrichment=WebGestaltR(enrichMethod="ORA", organism="mus_musculus", 
                         enrichDatabase="geneontology_Biological_Process_noRedundant",enrichDatabaseFile=NULL, 
                         enrichDatabaseType=NULL,enrichDatabaseDescriptionFile=NULL,interestGeneFile=NULL, 
                         interestGene=as.vector(genes),interestGeneType="genesymbol",collapseMethod="mean",referenceGeneFile=NULL,
                         referenceGene=NULL,referenceGeneType=NULL,referenceSet="genome_protein-coding", minNum=10, maxNum=500,
                         fdrMethod="BH",sigMethod="top",fdrThr=0.05,topThr=topN,dNum=20,perNum=1000,
                         lNum=20,is.output=TRUE,outputDirectory=getwd(),projectName=paste(c(base,"_GO"),collapse=""),keepGSEAFolder=FALSE,
                         methodType="R",dagColor="continuous",hostName="http://www.webgestalt.org/")
}
# 1. make a subset of the result based on the cutoff
sorted=enrichment[order(enrichment$FDR),]
dfDEpath<-enrichment[1:topN,]
nNode<-nrow(dfDEpath)
# what is the smallest FDR except for 0?
#smallestFDR<-min(dfDEpath$FDR[dfDEpath$FDR>0])
dfDEpath$FDR<- -log10(dfDEpath$FDR)

# 2. create the dfNodes and dfEdges data frame
dfNodes<-data.frame(id="",description="",size="",numDEgenes="",FDR="")
dfEdges<-data.frame(from="",to="",overlap="",symbols="")

GeneList<-list()
SymbolList<-list()

for(i in seq(1,nNode)) {
  node<-dfDEpath$geneset[i]
  GeneList[[node]]<-strsplit(dfDEpath[i,10],split=";",fixed=T)[[1]]
  SymbolList[[node]]=strsplit(dfDEpath[i,11],split=";",fixed=T)[[1]]
}

dfNodes<-dfDEpath[,c(1,2,4,5,9,11)]
colnames(dfNodes)[1]<-"id"

for(i in seq(1,nNode-1)) {
  for(j in seq(min(i+1,nNode),nNode)) {
    #cat(paste0("i=",i," j=",j,"\n"))
    node1<-dfDEpath$geneset[i]
    node2<-dfDEpath$geneset[j]
    overlapGenes<-length(intersect(GeneList[[node1]],GeneList[[node2]]))
    if(overlapGenes > 0) {
      overlapSymbols<-intersect(SymbolList[[node1]],SymbolList[[node2]])
      edgeline=c(node1,node2,overlapGenes,paste(overlapSymbols,collapse = ","))
      names(edgeline)<-colnames(dfEdges)
      dfEdges<-rbind(dfEdges,edgeline)
    }
  }
}
dfEdges<-dfEdges[-1,]
dfEdges$overlap<-as.numeric(dfEdges$overlap)
dfNodes$size=5+as.numeric(dfNodes$O)
dfNodes$symbols=dfNodes$OverlapGene_UserID
#dfEdges<-subset(dfEdges,overlap>0) # only leave the connected edges


MapNodeColor<-function(dfNodes,discrete=F,limit=10,low="blue",mid="white",high="red") {
  # usually bins=limit*2+1
  bins=length(seq(0,limit,0.1))
  FDRpal<-grDevices::colorRampPalette(c(low,mid,high))(bins)
  # set the upper limit for the FDR value
  dfNodes$FDR <- sapply(dfNodes$FDR,function(x) min(limit,x)) 
  # add the 0 for lower boundary
  #mappedColor<-FDRpal[cut(c(0,dfNodes$FDR),breaks=bins)]
  mappedColor=c()
  for( i in 1:dim(dfNodes)[1])
  {
    mappedColor[i]<-FDRpal[as.integer(bins*dfNodes$FDR[i]/limit)]
  }
 #dfNodes$color<-as.character(mappedColor[-1]) # remove the first color
  # plot the color key
  jpeg("network_colorKey.jpeg",width=12,height=6,units="cm",res=600)
  plot(seq(0,limit,0.1),rep(0,bins),col=FDRpal,pch=15,cex=4,xlab="-log10FDR",ylab="")
  dev.off()
  return(dfNodes)
}
# map the color of nodes based on FDR values. 
# try to detect limit automatically
limit=min(10,max(dfNodes$FDR))
limit=round(limit,1)
dfNodes<-MapNodeColor(dfNodes,limit=limit)
dfEdges$width<-sqrt(dfEdges$overlap)



ConstructNetwork<-function(NetList) {
  # NetList is a list of 2, in the order of dfNodes and dfEdges
  dfNodes<-NetList[[1]]
  dfEdges<-NetList[[2]]
  g<-graph_from_data_frame(dfEdges,directed=F,vertices=dfNodes) 
  # sometimes the lable is too long...
  V(g)$label<-V(g)$description
  V(g)$label<-sapply(V(g)$label,function(x) gsub("_","\n",x))
  V(g)$label<-sapply(V(g)$label,function(x) gsub(" ","\n",x))
  V(g)$title<-V(g)$symbols
  E(g)$label<-E(g)$overlap
  E(g)$title<-E(g)$symbols
  visGraph<- visIgraph(g,smooth=T,idToLabel=F,physics = T) %>% 
    visPhysics(solver = "forceAtlas2Based", 
               forceAtlas2Based = list(gravitationalConstant = -25)) %>%
    #visIgraphLayout(layout = "layout_in_circle") %>%
    visEdges(color = list(color = "dimgrey", highlight = "black"),smooth=T) %>%
    visNodes(font=list(color="black",size=10)) %>%
    visOptions(highlightNearest=list(enabled=T,degree=2),selectedBy="label") 
  return(visGraph)
}

g<-ConstructNetwork(list(dfNodes=dfNodes,dfEdges=dfEdges,width="100%",height="100%"))
outfile<-base
outfile<-paste0(outfile,".html")
visSave(g,file=outfile)


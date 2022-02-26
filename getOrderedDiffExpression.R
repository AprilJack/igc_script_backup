if (!require("EBSeqHMM")) {
  source("http://bioconductor.org/biocLite.R")
  biocLite("EBSeqHMM")
  library(EBSeqHMM)
  library(EBSeq)
}
#args<-c("raw_tiny.txt","y","y","m","m","o","o")
args<- commandArgs(trailingOnly = TRUE)
if (length(args)<3){
  cat("Usage: getOrderedDiffExpression raw.txt T1 T1 T2 T2 T3 T3 ...\n
      Will output some analyses and significant trajectories for each gene\n")
  if(length(args)>0)
  {
    cat(sprintf("datafile:%s\n",args[1]))
    data<-as.matrix(read.delim(file = args[1],header = TRUE, sep = "\t"))
    print(colnames(data[9:dim(data)[2]]),humanReadable = TRUE) 
  }
}else{
  fdr=0.05 # for DE testing we want to control fdr to at most this value
  cat("Running ordered differential analysis on Homer raw expression count table [col1=refseq,col8=annot.,cols9-N=raw counts] given the ordering: \n")
  Conditions<-factor(args[2:length(args)],levels=unique(args[2:length(args)]))
  data<-as.matrix(read.delim(file = args[1],header = TRUE, sep = "\t"))
  for( i in 1:dim(data)[1]) 
  {
    data[i,8]= strsplit(data[i,8],"[|]")[[1]][1]
    if( i %% 1000==0)
      cat(".")
  }
  data2 = data[,9:dim(data)[2]]
  data2= matrix(as.numeric(data2),nrow=dim(data2)[1],ncol=dim(data2)[2])
  rownames(data2) <- data[,8]
  colnames(data2) <- colnames(data[,9:dim(data)[2]])
  CondLabels = as.matrix(colnames(data2))
  while(length(args)-1 < dim(CondLabels)[1])
  {
    args[length(args)+1]<-"-"
  }
  rownames(CondLabels)<-args[2:(1+min(length(args)-1,dim(CondLabels)[1]))]
  print(CondLabels)
  if(dim(data2)[2] == length(args)-1)
  {
    
    #cat("Making genes human-readable")
    #for( i in 1:dim(data)[1]) 
    #{
    #  symbol= strsplit(toString(data[i,8]),"[|]")[[1]][1]
    #  rownames(data2)[i]=sprintf("%s_%s",symbol,toString(data[i,1]))
    #  if( i %% 1000==0)
    #    cat(".")
    #}
    rownames(data2)<- data[,1]
    cat("\nRunning EBSeqHMM on your ordered dataset to find significant pathways across conditions...\n")
    Sizes <- MedianNorm(data2)
    EBSeqHMMGeneOut <- EBSeqHMMTest(Data=as.matrix(data2), sizeFactors=Sizes, Conditions=Conditions, UpdateRd = 10,FCV=c(1.5,1.75,2))
    GeneDECalls <- GetDECalls(EBSeqHMMGeneOut, FDR=fdr)
    #cat("Top 25 significant genes:\n")
    #head(GeneDECalls,n = 25)
    GeneConfCalls <- GetConfidentCalls(EBSeqHMMGeneOut, FDR=fdr,cutoff=0.5, OnlyDynamic=TRUE)    
    write.table(file=sprintf("%s_paths_%s.txt",substr(args[1],1,nchar(args[1])-4),paste(levels(Conditions),collapse='')),x=GeneConfCalls[[1]],sep="\t")  
    for( i in 1:length(GeneConfCalls$EachPath))
    {
      cat(sprintf("Path: %s had %d sig genes with FDR=%3.3f\n",GeneConfCalls$EachPath[i][[1]][2],dim(GeneConfCalls$EachPath[i][[1]])[1],fdr))
      dat=as.data.frame(GeneConfCalls$EachPath[i])
      if(GeneConfCalls$NumEach[i]>0)
      {
        colnames(dat)<-c("Path","Prob(Path)")
        write.table(file=sprintf("%s_%s_%s.txt",args[1],GeneConfCalls$EachPath[i][[1]][2],paste(levels(Conditions),collapse='')),dat,sep = "\t",col.names = NA,quote = FALSE)
      }
    }
  }else{
    cat(sprintf("Error: the number of columns in data: %d does not match the group label count of %d\n",dim(data2)[2],length(args)-1))
    cat("\nUsage: getOrderedDiffExpression raw.txt T1 T1 T2 T2 T3 T3 ...\n
    Will output a list of genes that follow significant dynamic paths to raw_paths.txt\n")
  }
}
print(.libPaths())
print(sessionInfo())
print(version)

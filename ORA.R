#WebGestalt R wrapper for GALAXY
library(WebGestaltR)
ags<- commandArgs(trailingOnly = TRUE)
if (length(ags)<4){
  cat("Usage: ORA gene_list organism outdir database [useTop10]\n
      gene_list should be a HOMER analyzeRepeats table or it will just use the first column\n");
  cat("\norganisms:\n")
  cat(listOrganism())
  cat("\ndatabases:\n")
  cat(listGeneSet()$name)
  cat("\n")
  quit("no")
}
data <- as.matrix(read.delim(file = ags[1],header = TRUE, sep = "\t",quote = "",blank.lines.skip = TRUE,stringsAsFactors = FALSE))
cat(sprintf("Read in a table containing %d rows and %d columns\n",dim(data)[1],dim(data)[2]))
geneType="genesymbol"
if(startsWith(x=colnames(data)[1], prefix="Transcript.RepeatID"))
{
  for( i in 1:dim(data)[1]) 
  {
    data[i,8]= strsplit(data[i,8],"[|]")[[1]][1]
    if( i %% 1000==0)
      cat(".")
  }
  genes=as.vector(data[,8])
  cat(sprintf("Assuming HOMER formated table which had %d gene symbols\n",length(genes)))
  head(genes)
} else {
  cat("Assuming just one column for ids and one row for sample names\n")
  genes=data[,1]
  head(genes)
  max=0
  geneType="genesymbol"
  #types=c("genesymbol","refseq_mrna","unigene")
  #write(sprintf("Checking for best id type among %d genes...\n",length(genes)),stdout())
  #for(type in types)
  #{
  #  map=idMapping(organism = ags[2], dataType = "list",inputGene = genes,sourceIdType = type,mappingOutput = FALSE,inputGeneFile=NULL)
  #  if(length(map)>1 && dim(map$mapped)[1] > max)
  #  {
  #    max= dim(map$mapped)[1]
  #    geneType=type
  #    write(sprintf("Mapped %d genes assuming %s\n",max,type),stdout())
  #  }
  #}
}
if(length(ags)==4) {
  WebGestaltR(organism = ags[2],enrichMethod = "ORA", enrichDatabase = ags[4],interestGene = genes,interestGeneType = geneType ,referenceSet = "genome_protein-coding",outputDirectory  = ags[3],projectName=sprintf("%s-%s",ags[1],ags[4]),minNum = 5)
} else {
  WebGestaltR(organism = ags[2],enrichMethod = "ORA", enrichDatabase = ags[4],interestGene = genes,interestGeneType = geneType ,referenceSet = "genome_protein-coding",outputDirectory  = ags[3],projectName=sprintf("%s-%s",ags[1],ags[4]),minNum = 5, sigMethod = "top")
} 


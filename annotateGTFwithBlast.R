ags<- commandArgs(trailingOnly = TRUE)
if (length(ags)<2){
  cat("Usage: annotateGTFwithBlast gtf blastOut\n")
  cat("blastOut should be -outFmt 10, where the first value is the gtf ID and the second value is the blast hit\n")
  quit(save = "no",status = 1)
}

dict=read.csv(ags[2],stringsAsFactors = FALSE,header = FALSE)
dict=dict[!duplicated(dict[,1]),]
gtf=read.delim(ags[1],sep = "\t",stringsAsFactors = FALSE,header = FALSE,quote = "")
write(sprintf("Processing %s using %s",ags[1],ags[2]),stderr())
for(i in 1:dim(dict)[1])
{
  hit=grep(dict[i,1],gtf$V9)
  for(j in 1:length(hit)){
    gtf$V9[hit[j]]=sprintf("%s blast_hit \"%s\";",gtf$V9[hit[j]],dict[i,2])
  }
  if( as.integer((i)*100/dim(dict)[1]) >as.integer((i-1)*100/dim(dict)[1]) ){
    write(as.integer((i)*100/dim(dict)[1]),stderr())
    
  }

}
write.table(gtf,row.names=FALSE,col.names = FALSE,sep = "\t",quote = FALSE)
cat("Done!\n")
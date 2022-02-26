#setwd("/gpfs/analyses/maxsh/IGC/TCGA/htseqCounts/")
ags<- commandArgs(trailingOnly = TRUE)
if(length(ags) <1)
{
  cat("Usage: grabTCGA manifestFile\n
      \nCreates a File_metadata.txt with metadata from TCGA for files in manifest.\n")
}
x=read.table(ags[1],header = T)
manifest_length= nrow(x)
id= toString(sprintf('"%s"', x$id))
Part1= '{"filters":{"op":"in","content":{"field":"files.file_id","value":[ '
Part2= '] }},"format":"TSV","fields":"file_id,file_name,cases.submitter_id,cases.case_id,data_category,data_type,cases.samples.tumor_descriptor,cases.samples.tissue_type,cases.samples.sample_type","size":'
Part3= paste0("\"",manifest_length, "\"", "}") #just change single quote ' to double quote "
Sentence= paste(Part1,id,Part2,Part3, collapse=" ")
write.table(Sentence,"Payload.txt",quote=F,col.names=F,row.names=F)
system('curl --request POST --header "Content-Type: application/json" --data @Payload.txt "https://api.gdc.cancer.gov/files" > File_metadata.txt')
meta=read.table("File_metadata.txt",header=TRUE,sep = "\t",stringsAsFactors = FALSE)

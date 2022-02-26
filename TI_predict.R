#Script for computing TI for a given dataset
args<- commandArgs(trailingOnly = TRUE)
if(length(args) < 1)
{
  cat("Usage: TI TPM.txt [model]\n")
  exit(1)
} 
model="PtenKO_rf.R"
if(length(args)>1)
{
  model=args[2]
}


#Convert ILMN ids to Gene Symbols!
if (!require("illuminaHumanv4.db")) {
  source("http://bioconductor.org/biocLite.R")
  biocLite("beadarray")
  biocLite("illuminaHumanv4.db")
  library(illuminaHumanv4.db)
}
#illuminaHumanv4()
cat("Please paste your ILMN ids and type Ctrl-D.\n");
ids <- readLines("stdin",n=99999);
as.matrix(unlist(mget(ids,illuminaHumanv4SYMBOL)))
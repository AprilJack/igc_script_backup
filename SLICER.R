if (!require("minimist")) {
  install.packages("minimist", repos ="http://cran.cnr.berkeley.edu/", dependencies= TRUE)
  library(minimist)
}

if (!require("SLICER")) {
  library("devtools")
  install_github("jw156605/SLICER")
  library(SLICER)
}

ags=minimist(c("raw_5-9-17_justAge.txt"))
if (length(ags$'_')==0){
  cat("Usage: SLICER raw.txt \nWill filter low-count samples, normalize samples, and filter out non-variable genes.\nThen will run SLICER on the data and output a trajectory and ordering of the samples")
  ags=c("raw.txt")
}
cat("Loading in the data...\n")
data <- as.matrix(read.delim(file = ags$`_`,header = TRUE, sep = "\t"))
cat("Removing duplicates and normalizing the raw counts...\n")
# remove duplicates
#data <-subset(data,!duplicated(data[,1]))
# remove row headers
if(dim(data)[2]>8 & !is.numeric(data[8]))
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

data3 = data2%*% diag(1000000/colSums(data2))
data3=data3[order(rowSums(data2)),]
colnames(data3)=colnames(data2)
data4=log(data3[1:min(500,dim(data3)[1]),]+1,base=2)
for(i in 1:length(colnames(data3)))
{
  colnames(data4)[i]=strtrim(colnames(data4)[i],10)
}
set.seed(42)
cat("Starting SLICER on the top 500 genes\n")
selectedk = select_k(t(data4), kmin=10, kmax=25,by=1)
cat(sprintf("Finished finding k which is %d\n",selectedk))
traj_lle = lle(t(data4), m=2, selectedk,iLLE = TRUE)$Y
cat(sprintf("Finished lle\n"))
traj_graph = conn_knn_graph(traj_lle,selectedk)
cat(sprintf("Finished conn knn graph calculations\n"))
ends = find_extreme_cells(traj_graph, traj_lle)
cat(sprintf("Found extreme cells\n"))
graph_process_distance(traj_graph,start=1)
cells_ordered = cell_order(traj_graph, start)
cat(sprintf("Finished cell order calculations... Done with everything\n"))

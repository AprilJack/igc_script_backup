#
# wgcna.R
# 20170714
# Shawn Driscoll
#
# wrapper to the WGCNA pipeline so i don't have to remember how to use it.
#

#
# example:
# 
# run 'wgcna_init()' followed by 'wgcna_choose_soft_threshold'
# select a power where the left plot seems to have stabilized. the lowest value
# possible.
# run 'wgcna_build_net' to do all the clustering magic.
#
# WGCNA expects genes are in columns and that we are 
# interested in clustering genes (not samples)
################################
#
# Edited by Max Shokhirev 2018
#
# So first load in your dataset such that you have rows of genes and columns with samples/cells
# Next, inverse the matrix to with data2=t(data)
# Then, plug into pickSoftThreshold as x
# Then, pick a threshold power (e.g. one that is stable), and run the wgcna_build_net with x and new power. 
#
################################
if (!require("WGCNA")) {
  source("http://bioconductor.org/biocLite.R") 
  biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore")) 
  install.packages("WGCNA")
}
wgcna_init <- function(p=4) {
	require(WGCNA)
	allowWGCNAThreads(nThreads=p)
}

wgcna_choose_soft_threshold <- function(X, max.power=20) {
	
	# Choose a set of soft-thresholding powers
	powers = c(c(1:10), seq(from = 12, to=20, by=2))
	# Call the network topology analysis function
	sft = pickSoftThreshold(X, powerVector = powers, verbose = 5)
	# Plot the results:
	sizeGrWindow(9, 5)
	par(mfrow = c(1,2));
	cex1 = 0.9;
	# Scale-free topology fit index as a function of the soft-thresholding power
	plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
			xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
			main = paste("Scale independence"));
	text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
			labels=powers,cex=cex1,col="red");
	# this line corresponds to using an R^2 cut-off of h
	abline(h=0.90,col="red")
	# Mean connectivity as a function of the soft-thresholding power
	plot(sft$fitIndices[,1], sft$fitIndices[,5],
			xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
			main = paste("Mean connectivity"))
	text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")	
	
}

wgcna_build_net <- function(X, power=6, TOMType="unsigned", minModuleSize=100, 
		reassignThreshold=0, mergeCutHeight=0.25, numericLabels=TRUE, 
		pamRespectsDendro=FALSE, saveTOMs=TRUE, saveTOMFileBase="wgcnaToms", 
		verbose=3, ...) {
	
	
	net <- blockwiseModules(X, power = power,
			TOMType = TOMType, minModuleSize = minModuleSize,
			reassignThreshold = reassignThreshold, mergeCutHeight = mergeCutHeight,
			numericLabels = numericLabels, pamRespectsDendro = pamRespectsDendro,
			saveTOMs = saveTOMs,
			saveTOMFileBase = saveTOMFileBase,
			verbose = verbose, ...)	
	
	return(net)
	
}

wgcna_groups <- function(net) {
	return(net$colors)
}

wgcna_plot_dendro <- function(net) {
	
	# open a graphics window
	sizeGrWindow(12, 9)
	# Convert labels to colors for plotting
	mergedColors = labels2colors(net$colors)
	# Plot the dendrogram and the module colors underneath
	plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
			"Module colors",
			dendroLabels = FALSE, hang = 0.03,
			addGuide = TRUE, guideHang = 0.05)	
	
}



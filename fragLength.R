#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

frag=read.csv(args[1],sep=' ')
frag=apply(as.matrix(frag), 2, as.numeric)
frag=frag[,c(2,1)]
frag[,2]=frag[,2]/sum(frag[,2])*100
colnames(frag)=c("Fragment Length","Frequency")
pdf(sprintf("%s.pdf",strtrim(args[1],nchar(args[1])-4)))
plot(frag,xlim = c(0,1000),type='l',main=args[1])
dev.off()

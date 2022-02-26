
ags<- commandArgs(trailingOnly = TRUE)
if(length(ags)==0)
{
  cat("Venn will draw a 2,3, or 4-way Venn diagram for a HOMER mergePeaks venn file (e.g. specify -venn) or a file with one column of values per condition.\n")
  cat("Usage: Venn file1.venn [file2.venn] [file3.venn] ...\n")
  cat("Usage: Venn Lists.txt [Lists2.txt] [Lists2.txt] ...\n")
  quit(save="no")
}
if (!require("VennDiagram")) {
   install.packages("VennDiagram",repos = "http://cran.us.r-project.org")
}
library("VennDiagram")
for( i in ags)
{
  data<-read.delim(file = i,header = TRUE,sep = "\t")
  data[is.na(data)]<-""
  if(colnames(data)[length(colnames(data))] == "Name" & colnames(data)[length(colnames(data))-1] == "Total")
  {
    dims=length(colnames(data))-2
    if(dims==2)
    {
      n2=0
      n1=0
      n12=0
      for(j in 1:dim(data)[1])
      {
        if(sum(data[j,1:2]==c("","X"))==2)
        {
          #n2
          n2=data[j,3]
        } else if(sum(data[j,1:2]==c("X",""))==2)
        {
          #n1
          n1=data[j,3]
        }
        else if(sum(data[j,1:2]==c("X","X"))==2)
        {
          #n12
          n12=data[j,3]
        }
      }
      area1=n1+n12;
      area2=n2+n12;
      png(file=sprintf("%s.png",i),width = 1000,height = 1000)
      draw.pairwise.venn(area1,area2,cat.just=list(c(0,1),c(1,0)),n12,category = colnames(data)[1:2],fill = c("red","blue"),alpha=c(0.5,0.5),euler.d = TRUE,ext.text = TRUE,col = c("red","blue"))
      dev.off()
    } else if (dims==3) {
      n3=0
      n2=0
      n1=0
      n12=0
      n13=0
      n23=0
      n123=0
      for(j in 1:dim(data)[1])
      {
        if(sum(data[j,1:3]==c("","","X"))==3){
          #n3
          n3=data[j,4]
        } else if(sum(data[j,1:3]==c("","X",""))==3){
          #n2
          n2=data[j,4]
        }else if(sum(data[j,1:3]==c("X","",""))==3){
          #n1
          n1=data[j,4]
        }else if(sum(data[j,1:3]==c("","X","X"))==3){
          #n23
          n23=data[j,4]
        }else if(sum(data[j,1:3]==c("X","","X"))==3){
          #n13
          n13=data[j,4]
        }else if(sum(data[j,1:3]==c("X","X",""))==3){
          #n12
          n12=data[j,4]
        }else if(sum(data[j,1:3]==c("X","X","X"))==3){
          #n123
          n123=data[j,4]
        }
      }
      area1=n1+n12+n13+n123;
      area2=n2+n12+n23+n123;
      area3=n3+n13+n23+n123;
      nn12=n12+n123
      nn13=n13+n123
      nn23=n23+n123
      nn123=n123
      png(file=sprintf("%s.png",i),width = 1000,height = 1000)
      draw.triple.venn(area1,area2,area3,cat.just=list(c(0,1),c(1,0),c(0.5,0)),n12 = nn12,n13 = nn13,n123=nn123,n23=nn23,category = colnames(data)[1:3],fill = c("red","blue","green"),alpha=c(0.5,0.5,0.5),euler.d = TRUE,ext.text = TRUE,col = c("red","blue","green"))
      dev.off()
    }  else if (dims==4) {
      n4=0
      n3=0
      n2=0
      n1=0
      n12=0
      n13=0
      n14=0
      n23=0
      n24=0
      n34=0
      n123=0
      n134=0
      n234=0
      n124=0
      n1234=0
      for(j in 1:dim(data)[1])
      {
        if(sum(data[j,1:4]==c("","","","X"))==4)
        {
          #n4
          n4=data[j,5]
        } else if(sum(data[j,1:4]==c("","","X",""))==4){
          #n3
          n3=data[j,5]
        } else if(sum(data[j,1:4]==c("","","X","X"))==4) {
          #n34
          n34=data[j,5]
        } else if(sum(data[j,1:4]==c("","X","",""))==4) {
          #n2
          n2=data[j,5]
        }else if(sum(data[j,1:4]==c("","X","","X"))==4) {
          #n24
          n24=data[j,5]
        }else if(sum(data[j,1:4]==c("","X","X",""))==4) {
          #n23
          n23=data[j,5]
        }else if(sum(data[j,1:4]==c("","X","X","X"))==4) {
          #n234
          n234=data[j,5]
        }else if(sum(data[j,1:4]==c("X","","",""))==4) {
          #n1
          n1=data[j,5]
        }else if(sum(data[j,1:4]==c("X","","","X"))==4) {
          #n14
          n14=data[j,5]
        }else if(sum(data[j,1:4]==c("X","","X",""))==4) {
          #n13
          n13=data[j,5]
        }else if(sum(data[j,1:4]==c("X","","X","X"))==4) {
          #n134
          n134=data[j,5]
        }else if(sum(data[j,1:4]==c("X","X","",""))==4) {
          #n12
          n12=data[j,5]
        }else if(sum(data[j,1:4]==c("X","X","","X"))==4) {
          #n124
          n124=data[j,5]
        }else if(sum(data[j,1:4]==c("X","X","X",""))==4) {
          #n123
          n123=data[j,5]
        }else if(sum(data[j,1:4]==c("X","X","X","X"))==4) {
          #n1234
          n1234=data[j,5]
        }
      }
      #n1+n2+n3+n4+n12+n13+n14+n23+n24+n34+n123+n234+n124+n134+n1234
      area1=n1+n12+n13+n14+n123+n124+n134+n1234
      area2=n2+n12+n23+n24+n123+n234+n124+n1234
      area3=n3+n13+n23+n34+n123+n234+n134+n1234
      area4=n4+n14+n24+n34+n234+n124+n134+n1234
      nn123=n123+n1234
      nn124=n124+n1234
      nn134=n134+n1234
      nn234=n234+n1234
      nn12=n12+n123+n124+n1234
      nn13=n13+n123+n134+n1234
      nn14=n14+n124+n134+n1234
      nn23=n23+n234+n123+n1234
      nn24=n24+n234+n124+n1234
      nn34=n34+n234+n134+n1234
      png(file=sprintf("%s.png",i),width = 1000,height = 1000)
      draw.quad.venn(area1=area1,area2=area2,area3=area3,area4=area4,cat.just=list(c(0,1),c(1,0),c(0.5,0),c(0.5,0)),n14=nn14,n12=nn12,n13=nn13,n23=nn23,n34=nn34,n24 = nn24,n123=nn123,n124=nn124,n134=nn134,n234=nn234,n1234=n1234,category = colnames(data)[1:4],fill = c("red","blue","green","cyan"),alpha=c(0.5,0.5,0.5,0.5),euler.d = TRUE,ext.text = TRUE,col = c("red","blue","green","cyan"))
      dev.off()
    }
  } else {
    #assume we just have lists
    cat("Processing some lists!\n")
    LS <- lapply(as.list(data), function(x) x[x != ""])
    names(LS)<-colnames(data)
    venn.plot <- venn.diagram(LS,file=sprintf("%s.png",i),imagetype = "png",main = sprintf("%s",i),fill=c("red","blue","green","yellow","cyan")[1:length(colnames(data))],alpha=c(0.5,0.5,0.5,0.5,0.5)[1:length(colnames(data))],euler.d = TRUE,ext.text = TRUE,margin=0.05)
    x=calculate.overlap(LS)
    x$a9=as.vector(x$a9)
    x$a14=as.vector(x$a14)
    x$a1=as.vector(x$a1)
    x$a2=as.vector(x$a2)
    x$a3=as.vector(x$a3)
    x$a4=as.vector(x$a4)
    x$a5=as.vector(x$a5)
    x$a6=as.vector(x$a6)
    x$a7=as.vector(x$a7)
    x$a8=as.vector(x$a8)
    x$a10=as.vector(x$a10)
    x$a11=as.vector(x$a11)
    x$a12=as.vector(x$a12)
    x$a13=as.vector(x$a13)
    mat=plyr::ldply(x, rbind)
    write.table(t(mat),file=sprintf("%s.txt",i),sep="\t",quote=FALSE)
  }
}

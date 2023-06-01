#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
if(length(args)<1) q()
libPath = tail(args,1)
if(nchar(libPath)>3){
  addPath <- unlist(strsplit(libPath,";"))
  addPath <- addPath[sapply(addPath,dir.exists)]
  .libPaths(c(addPath,.libPaths()))
}
suppressMessages(suppressWarnings(require(fgsea)))
suppressMessages(suppressWarnings(require(data.table)))
suppressMessages(suppressWarnings(require(ggplot2)))
options(bitmapType='cairo')
set.seed(0)
strCSV <- args[1]
strGMT <- args[2]
gsMin <- as.numeric(args[3])
gsMax <- as.numeric(args[4])
padjCut <- as.numeric(args[5])
upN <- as.numeric(args[6])
dnN <- as.numeric(args[7])
collapseF <- args[8]=='True'
strFun <- args[9]
fontsize <- as.numeric(args[10])
dpi <- as.numeric(args[11])

mtable <- read.csv(strCSV,as.is=T,check.names=F)
colnames(mtable) <- c('gene_name','logFC','pvalue','FDR')
mtable <- mtable[mtable$logFC!=0,]
gRank <- sort(setNames(mtable$logFC,mtable$gene_name))
GS <- gmtPathways(strGMT)

fgseaRes <- NULL
fgseaRes <- tryCatch(suppressMessages(suppressWarnings(fgsea(GS,gRank,nperm=1e5,minSize=gsMin,maxSize=gsMax))),error=function(e){return(NULL)})
if(!is.null(fgseaRes)) fwrite(fgseaRes,file=strCSV)#gsub("csv$","fgse.csv",strCSV)
if(sum(upN+dnN)==0 || is.null(fgseaRes) || nrow(fgseaRes)==0){
  cat("")
  q()
}

if(collapseF){
  collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj<padjCut],GS,gRank)
  fgseaRes <- fgseaRes[pathway %in% collapsedPathways$mainPathways]
}
topPathways <- c(fgseaRes[ES>0][head(order(pval),n=upN), pathway],
                 rev(fgseaRes[ES<0][head(order(pval), n=dnN), pathway]))
if(length(topPathways)==0){
  cat("")
  q()
}
pathwayW <- max(nchar(topPathways))/4.5
colwidths = c(pathwayW, 3.5, 1, 1.2, 1.2)
w = max(3,pathwayW/2)+4
h = max(6,6*length(topPathways)/20)

strImg <- gsub("csv$",strFun,strCSV)
f <- get(strFun)
if(sum(strFun%in%c('png','jpeg','tiff'))>0){
  f(strImg, width=w, height=h,units='in',res=dpi)
}else{
  f(strImg, width=w, height=h)
}
plotGseaTable(GS[topPathways], gRank, fgseaRes,gseaParam=0.5,colwidths=colwidths)
a <- dev.off()
fig = base64enc::dataURI(file = strImg)
cat(gsub("data:;base64,","",fig))
a <- file.remove(strImg)

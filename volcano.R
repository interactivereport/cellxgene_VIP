#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
if(length(args)<1) q()
libPath = tail(args,1)
if(nchar(libPath)>3){
  addPath <- unlist(strsplit(libPath,";"))
  addPath <- addPath[sapply(addPath,dir.exists)]
  .libPaths(c(addPath,.libPaths()))
}

library(ggplot2)
library(ggrepel)
library(ggrastr)
options(bitmapType='cairo')

strCSV = args[1]
genes <- unlist(strsplit(args[2],";"))
strFun <- args[3]
fontsize <- as.numeric(args[4])
dpi <- as.numeric(args[5])
logFCcut <- as.numeric(args[6])
strDN <- args[7]
strUP <- args[8]
FDR <- as.numeric(args[9])
logFC <- as.numeric(args[10])
labelSize <- as.numeric(args[11])
dotSize <- as.numeric(args[12])
ymin <- as.numeric(args[13]) 
ymax <- as.numeric(args[14]) 
rasterize <- ifelse(grepl(args[15],'Yes'), T, F)


mtable <- read.csv(strCSV,as.is=T,check.names=F)
colnames(mtable) <- c('gene_name','logFC','pvalue','FDR')
mtable$Top = ifelse(mtable$FDR >= FDR, "Not Sig", ifelse(mtable$logFC>logFC, "Up", ifelse(mtable$logFC< -logFC,"Down","Not Sig")))
Up <- length(which(mtable$Top=="Up"))
Down <- length(which(mtable$Top=="Down"))
# only label top 20 genes
genes <- unique(c(genes,mtable[1:20,"gene_name"]))
## remove genes with absolute logFC larger than the logFCcutoff
mtable <- mtable[abs(mtable$logFC)<logFCcut,]
#mtable[mtable$logFC>logFCcut,'logFC'] <- logFCcut
#mtable[mtable$logFC< -logFCcut,'logFC'] <- -logFCcut
mtable$FDR = ifelse(-log10(mtable$FDR)<ymax,mtable$FDR,1/10^ymax)

g <- ggplot(mtable, aes(x=logFC, y=-log10(FDR))) + ylim(c(ymin,ymax))

if (rasterize) {
  g <- g + geom_point_rast(aes(color = Top), size=dotSize, alpha=0.6,na.rm=TRUE)
} else {
  g <- g + geom_point(aes(color = Top), size=dotSize, alpha=0.6)
}

g <- g + theme_bw(base_size = 12) + xlab("log2(FC)") +
  scale_color_manual(values = c("Up"="#B41A29", "Not Sig"="grey", "Down"="#5B98E6")) +
  geom_text_repel(data = subset(mtable, gene_name %in% genes), aes(label = gene_name), size = labelSize, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"))
g <- g + theme(legend.position="none") + geom_hline(yintercept=-log10(FDR), linetype="dashed", color = "darkgreen")
xrange = ggplot_build(g)$layout$panel_scales_x[[1]]$range$range
yrange = ggplot_build(g)$layout$panel_scales_y[[1]]$range$range

g <- g + annotate(geom="text", x=ggplot_build(g)$layout$panel_params[[1]]$x.range[1]+1, y=-log10(FDR)+yrange[2]/80, label=paste0("FDR=",round(FDR,3)), color="darkgreen", size=labelSize, fontface="bold")
g <- g + annotate(geom="text", x=xrange[1]/2, y=ymin, label=paste("Up in",strDN,":",Down), color="#5B98E6", size=fontsize/3, fontface="bold")
g <- g + annotate(geom="text", x=xrange[2]/2, y=ymin, label=paste("Up in",strUP,":",Up), color="#B41A29", size=fontsize/3, fontface="bold")
g <- g + theme(panel.background = element_rect(fill = "transparent", color=NA),
               plot.background = element_rect(fill = "transparent", color = NA),
               text=element_text(size=fontsize))

strImg <- gsub("csv$",strFun,strCSV)
f <- get(strFun)
if(sum(strFun%in%c('png','jpeg','tiff'))>0){
  f(strImg, width=6.5, height=6.5,units='in',res=dpi)
}else{
  f(strImg, width=6.5, height=6.5)
}
print(g)
a <- dev.off()
fig = base64enc::dataURI(file = strImg)
cat(gsub("data:;base64,","",fig))
a <- file.remove(strImg)

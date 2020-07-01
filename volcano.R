#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
if(length(args)<1) q()
library(ggplot2)
library(ggrepel)
library(ggrastr)


strCSV = args[1]
genes <- unlist(strsplit(args[2],";"))
strFun <- args[3]
fontsize <- as.numeric(args[4])
dpi <- as.numeric(args[5])
mtable <- read.csv(strCSV,as.is=T,check.names=F)
colnames(mtable) <- c('gene_name','logFC','pvalue','FDR')
mtable$Top = ifelse(mtable$FDR >= 0.05, "Not Sig", ifelse(mtable$logFC>0, "Up", "Down"))
Up <- length(which(mtable$Top=="Up"))
Down <- length(which(mtable$Top=="Down"))
# only label top 20 genes
genes <- unique(c(genes,mtable[1:20,"gene_name"]))
g <- ggplot(mtable, aes(x=logFC, y=-log10(FDR))) +
  #geom_point(aes(color = Top), size=0.5, alpha=0.6) +
  geom_point_rast(aes(color = Top), size=0.5, alpha=0.6,na.rm=TRUE)+
  theme_bw(base_size = 12) + xlab("log2(FC)") +
  scale_color_manual(values = c("Up"="#B41A29", "Not Sig"="grey", "Down"="#5B98E6")) +
  geom_text_repel(data = subset(mtable, gene_name %in% genes), aes(label = gene_name), size = fontsize/4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"))
g <- g + theme(legend.position="none") + geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "darkgreen")
xrange = ggplot_build(g)$layout$panel_scales_x[[1]]$range$range
yrange = ggplot_build(g)$layout$panel_scales_y[[1]]$range$range

g <- g + annotate(geom="text", x=ggplot_build(g)$layout$panel_params[[1]]$x.range[1]+1, y=-log10(0.05)+yrange[2]/80, label="FDR=0.05", color="darkgreen", size=fontsize/4, fontface="bold")
g <- g + annotate(geom="text", x=xrange[1]/2, y=-5, label=paste("Down-regulated:",Down), color="#5B98E6", size=fontsize/3, fontface="bold")
g <- g + annotate(geom="text", x=xrange[2]/2, y=-5, label=paste("Up-regulated:",Up), color="#B41A29", size=fontsize/3, fontface="bold")
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

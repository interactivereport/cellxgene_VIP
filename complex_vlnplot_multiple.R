suppressMessages(suppressWarnings({
  require(ggplot2)
  require(reshape2)
  require(ggrastr)
  require(grid)
  require(ggbeeswarm)
}))
args = commandArgs(trailingOnly = TRUE)
if(length(args)<6){
  message("Error: too few input args, (",paste(args,collapse=", "),")")
}
libPath = tail(args,1)
if(nchar(libPath)>4){
  addPath <- unlist(strsplit(libPath,";"))
  addPath <- addPath[sapply(addPath,dir.exists)]
  .libPaths(c(addPath,.libPaths()))
}

strCSV <- args[1]
genes <- unlist(strsplit(args[2],","))
cluster <- args[3]
grp <- args[4]
imgW <- as.numeric(args[5])
imgH <- as.numeric(args[6])
#imgCol <- args[7]
i <- 7
strFun <- args[i]
fontsize <- as.numeric(args[i+1])
dpi <- as.numeric(args[i+2])
#swapAxes <- as.logical(args[i+3])
#rasterize <- ifelse(grepl(args[i+4],'Yes'), T, F)

gene_count <- as.data.frame(data.table::fread(strCSV))
cellID <- colnames(gene_count)[1]
cellN <- nrow(gene_count)
gene_count <- reshape2::melt(gene_count, id.vars = c(cellID,cluster,grp), measure.vars = genes,
                             variable.name = "Genes", value.name = "Expr")
gene_count$Genes<-factor(gene_count$Genes,levels = genes)
#
change_strip_background <- function(
  ggplt_obj,
  type = "top",
  strip.color=NULL
  ){
  g <- ggplot_gtable(ggplot_build(ggplt_obj))
  if(type == "top"){
    strip_both <- which(grepl('strip-t', g$layout$name))
    fills<-strip.color
    if(is.null(fills)){
    fills<- scales::hue_pal(l=90)(length(strip_both))
    }
  } else if(type=="right"){
    strip_both <- which(grepl('strip-r', g$layout$name))
    fills<-strip.color
    if(is.null(fills)){
      fills<- scales::hue_pal(l=90)(length(strip_both))
    }
  } else {
    strip_t <- which(grepl('strip-t', g$layout$name))
    strip_r <- which(grepl('strip-r', g$layout$name))
    strip_both<-c(strip_t, strip_r)
    fills<-strip.color
    if(is.null(fills)){
      fills<- c(scales::hue_pal(l=90)(length(strip_t)),scales::hue_pal(l=90)(length(strip_r)))
    }
  }
  k <- 1
  for (i in strip_both) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  g
}
## plot
alpha <- 0.01
strip.color <- NULL
font.size <- fontsize
pt.size <- 0.1
p<-ggplot(gene_count, aes_string(grp, "Expr", fill = grp)) +
  geom_violin(scale = 'width', adjust = 1, trim = TRUE, size=0.3, alpha=0.5, color="pink") +
  rasterise(geom_quasirandom(size=pt.size, alpha=alpha),dpi=dpi)+
  #geom_quasirandom(size=pt.size, alpha=alpha)+
  scale_y_continuous(expand = c(0, 0), position="left", labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(as.formula(paste(cluster,"~Genes")), scales =  'free_x') +
  xlab("") + ylab("")+ggtitle(paste(cellN,"cells"))+
  theme(panel.background = element_rect(fill = "white",colour = "black"),
        axis.title = element_text(size = font.size),
        axis.text.x = element_text(size = font.size, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size=(font.size)),
        strip.text = element_text( size = font.size-3),
        legend.title = element_blank(),
        legend.position = 'none')
g <- change_strip_background(p, type = 'both',  strip.color = strip.color)

strImg <- gsub("csv$",strFun,strCSV)
f <- get(strFun)
f(strImg, width=imgW, height=imgH,units='in',res=dpi)
grid::grid.draw(g)
a<-dev.off()
fig = base64enc::dataURI(file = strImg)
cat(gsub("data:;base64,","",fig))
a <- file.remove(strImg)

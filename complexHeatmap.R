#!/usr/bin/env Rscript
#devtools::install_version("rjson",version="0.2.20",repos="https://cran.us.r-project.org")
#BiocManager::install("ComplexHeatmap")
#options(stringsAsFactors=F)
# ./complexHeatmap.R /tmp/HEAT1661227613.365820.csv CD4,BTK,CSF2 Expression celltype.l1,celltype.l2 6 8 Reds png 11 150 F Yes ''
a<- suppressWarnings(suppressMessages({
require(dplyr)
graphics.off()
args = commandArgs(trailingOnly = TRUE)
if(length(args)<1) q()
libPath = tail(args,1)
if(nchar(libPath)>3){
  addPath <- unlist(strsplit(libPath,";"))
  addPath <- addPath[sapply(addPath,dir.exists)]
  .libPaths(c(addPath,.libPaths()))
}

library(ComplexHeatmap)
library(png)

strCSV <- args[1]
genes <- unlist(strsplit(args[2],","))
heatkey <- args[3]
clusterGrps <- unlist(strsplit(args[4],",")) #ifelse(is.na(as.logical(args[4])),T,as.logical(args[4]))
imgW <- as.numeric(args[5])
imgH <- as.numeric(args[6])
imgCol <- args[7]
#imgColRev <- ifelse(is.na(as.logical(args[8])),F,as.logical(args[8]))
legendRow <- as.numeric(args[8])
annFontsize <- as.numeric(args[9])
i <- 10
strFun <- args[i]
fontsize <- as.numeric(args[i+1])
dpi <- as.numeric(args[i+2])
#columnFormat <- args[i+3]
#rowFormat <- args[i+4]
#annoFormat <- args[i+5]
swapAxes <- as.logical(args[i+3])
rasterize <- ifelse(grepl(args[i+4],'Yes'), T, F)

imgColRev <- F
if(grepl("_r$",imgCol)){
  imgColRev <- T
  imgCol <- gsub("_r$","",imgCol)
}

D <- as.data.frame(data.table::fread(strCSV))
# if clusterRow is True, remove cells with sd=0 and value<1
clusterRow <- ifelse("Expression"%in%clusterGrps,T,F)
if(length(genes)<2) clusterRow <- F
if(clusterRow){
  mat <- as.matrix(D[,genes,drop=F])
  mat <- mat[ , colSums(is.na(mat)) == 0]
  matSD <- apply(mat,1,sd)
  matExp <- mat[,1]
  sel <- matSD>0 & matExp>1
  mat <- mat[sel,]
  D <- D[sel,]
}else{
  if(length(clusterGrps)>0){
    D <- D%>%arrange(!!!syms(clusterGrps))
    D <- D[,c(clusterGrps,setdiff(colnames(D),clusterGrps)),drop=F]
  }

  mat <- as.matrix(D[,genes,drop=F])
  mat <- mat[ , colSums(is.na(mat)) == 0]
}

## annotation
anno_df <- D[,!colnames(D)%in%genes,drop=F]
anno_df <- anno_df[,sapply(anno_df,function(x)return(length(unique(x))))<50,drop=F]
ann_col <- apply(anno_df,2,function(x){
  if(is.numeric(x)) {
    if(prod(range(x))<0)
      return(list(circlize::colorRamp2(c(min(x),0,max(x)),RColorBrewer::brewer.pal(name=imgCol,n=3))))
    else
      return(list(circlize::colorRamp2(c(min(x),median(x),max(x)),RColorBrewer::brewer.pal(name=imgCol,n=3))))
  }
  xN <- length(unique(x))
  if(xN<13)
    return(list(setNames(sample(RColorBrewer::brewer.pal(name="Set3",n=12),xN),unique(x))))
  return(list(setNames(scales::hue_pal()(xN),sort(unique(x)))))
})
if (!swapAxes) {
  rowAnno <- HeatmapAnnotation(df=anno_df,col=unlist(ann_col,recursive=F),
                         annotation_name_gp = grid::gpar(fontsize=max(3,fontsize=fontsize+annFontsize)), #eval(parse(text = paste0("grid::gpar(", annoFormat , ")"))),
                         annotation_legend_param=list(title_gp=gpar(fontsize = fontsize+annFontsize+2),
                                                      labels_gp=gpar(fontsize=fontsize+annFontsize+1)),which = "row")
} else {
  colAnno <- HeatmapAnnotation(df=anno_df,col=unlist(ann_col,recursive=F),
                         annotation_name_gp = grid::gpar(fontsize=max(3,fontsize=fontsize+annFontsize)), #eval(parse(text = paste0("grid::gpar(", annoFormat , ")"))),
                         annotation_legend_param=list(title_gp=gpar(fontsize = fontsize+annFontsize+2),
                                                      labels_gp=gpar(fontsize=fontsize+annFontsize+1),
                                                      nrow=legendRow))
  mat <- t(mat)
}
## color for heatmap
colorN <- 21 # please be an odd number
colorSet <- RColorBrewer::brewer.pal(name = imgCol, n = colorN)
if(length(colorSet)%%2==0){
 colorSet <- c(colorSet,tail(colorSet,1))
}
colorN <- length(colorSet)
if(imgColRev) colorSet <- rev(colorSet)
matRange <- range(mat)
if(prod(matRange)<0){
  if(max(abs(matRange))>30){
    colorBreak <- c(head(quantile(mat[mat<0],seq(0,1,length.out=ceiling(colorN/2))),-1),0,
                    quantile(mat[mat>0],seq(0,1,length.out=ceiling(colorN/2)))[-1])
  }else{
    colorBreak <- c(head(seq(matRange[1],0,length.out=ceiling(colorN/2)),-1),0,
                    seq(0,matRange[2],length.out=ceiling(colorN/2))[-1])
  }
}else{
  if(max(abs(matRange))>30){
    colorBreak <- quantile(mat[mat!=0],seq(0,1,length.out=colorN))
  }else{
    colorBreak <- seq(matRange[1],matRange[2],length.out=colorN)
  }

}
col_fun <- circlize::colorRamp2(colorBreak,colorSet)
}))

strImg <- gsub("csv$",strFun,strCSV)
f <- get(strFun)
if(sum(strFun%in%c('png','jpeg','tiff'))>0){
  SVGformat <- F
}else{
  SVGformat <- T
}

if (!swapAxes) {
  annotation_legend_side <- "right"
  p <- Heatmap(mat,name=heatkey,right_annotation=rowAnno,
#             col=col_fun,column_title = paste(nrow(mat),"cells"),
             col=col_fun,
             column_names_gp = grid::gpar(fontsize=max(2,fontsize=fontsize+annFontsize)),#eval(parse(text = paste0("grid::gpar(", columnFormat , ")"))),
             row_names_gp = grid::gpar(fontsize=max(2,fontsize=fontsize+annFontsize)),#eval(parse(text = paste0("grid::gpar(", rowFormat, ")"))),
             cluster_columns=F,cluster_rows=clusterRow,clustering_method_rows="ward.D",
             heatmap_legend_param=list(title_gp=gpar(fontsize = fontsize),
                                       labels_gp=gpar(fontsize=fontsize-1)))
} else {
#print(dim(mat))
#print(attributes(colAnno))
annotation_legend_side <- "top"
  p <- Heatmap(mat,name=heatkey,top_annotation=colAnno,
             col=col_fun,
             column_names_gp = grid::gpar(fontsize=max(2,fontsize=fontsize+annFontsize)),#eval(parse(text = paste0("grid::gpar(", rowFormat , ")"))),
             row_names_gp = grid::gpar(fontsize=max(2,fontsize=fontsize+annFontsize)),#eval(parse(text = paste0("grid::gpar(", columnFormat, ")"))),
             cluster_rows=F,cluster_columns=clusterRow,clustering_method_columns="ward.D",
             heatmap_legend_param=list(title_gp=gpar(fontsize = fontsize+annFontsize+2),
                                       labels_gp=gpar(fontsize=fontsize+annFontsize+1)))
}
if (SVGformat) {
  if (rasterize) {
    pngFile = paste0(strImg,".png")
    png(pngFile, width=imgW, height=imgH, units='in', res=dpi, bg='transparent')
    draw(p,annotation_legend_side=annotation_legend_side)

    for (comp in grid.ls(flatten=TRUE,print=FALSE,recursive=FALSE)[1]$name) {
      if (!grepl("GRID.rect",comp)) {
        grid.remove(comp)
      }
    }
    a <- dev.off()
    heatmap = readPNG(pngFile)
  }

  f(strImg, width=imgW, height=imgH)
  draw(p,annotation_legend_side=annotation_legend_side)

  if (rasterize) {
    for (comp in grid.ls(flatten=TRUE,print=FALSE,recursive=FALSE)[1]$name) {
      if (grepl("GRID.rect",comp)) {
        grid.remove(comp)
      }
    }
    grid.raster(heatmap)
    a <- file.remove(pngFile)
  }
} else {
  f(strImg, width=imgW, height=imgH,units='in',res=dpi)
  draw(p,annotation_legend_side=annotation_legend_side)
}

a <- dev.off()
fig = base64enc::dataURI(file = strImg)
cat(gsub("data:;base64,","",fig))
a <- file.remove(strImg)

#!/usr/bin/env Rscript
#devtools::install_version("rjson",version="0.2.20",repos="https://cran.us.r-project.org")
#BiocManager::install("ComplexHeatmap")
#options(stringsAsFactors=F)
# ./complexHeatmap.R HEAT1650501335.568982.csv Cd9,Neurod2,Plp1,Mbp,Mag Expression T 6 8 RdBu pdf 8 300 ''
a<- suppressWarnings(suppressMessages({
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

strCSV <- args[1]
genes <- unlist(strsplit(args[2],","))
heatkey <- args[3]
clusterRow <- ifelse(is.na(as.logical(args[4])),T,as.logical(args[4]))
imgW <- as.numeric(args[5])
imgH <- as.numeric(args[6])
imgCol <- args[7]
#imgColRev <- ifelse(is.na(as.logical(args[8])),F,as.logical(args[8]))
i <- 8
strFun <- args[i]
fontsize <- as.numeric(args[i+1])
dpi <- as.numeric(args[i+2])

imgColRev <- F
if(grepl("_r$",imgCol)){
  imgColRev <- T
  imgCol <- gsub("_r$","",imgCol)
}

D <- as.data.frame(data.table::fread(strCSV))
mat <- as.matrix(D[,genes,drop=F])
# if clusterRow is True, remove cells with sd=0 and value<1
if(ncol(mat)==1) clusterRow <- F
if(clusterRow){
  matSD <- apply(mat,1,sd)
  matExp <- mat[,1]
  sel <- matSD>0 & matExp>1
  mat <- mat[sel,]
  D <- D[sel,]
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
rowAnno <- rowAnnotation(df=anno_df,col=unlist(ann_col,recursive=F),
                         annotation_legend_param=list(title_gp=gpar(fontsize = fontsize),
                                                      labels_gp=gpar(fontsize=fontsize-1)))
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
  f(strImg, width=imgW, height=imgH,units='in',res=dpi)
  use_raster <- F
}else{
  f(strImg, width=imgW, height=imgH)
  use_raster <- T
}

p <- Heatmap(mat,name=heatkey,right_annotation=rowAnno,
             col=col_fun,use_raster=use_raster,column_title = paste(nrow(mat),"cells"),
             cluster_columns=F,cluster_rows=clusterRow,clustering_method_rows="ward.D",
             heatmap_legend_param=list(title_gp=gpar(fontsize = fontsize),
                                       labels_gp=gpar(fontsize=fontsize-1)))
draw(p)
a <- dev.off()
fig = base64enc::dataURI(file = strImg)
cat(gsub("data:;base64,","",fig))
a <- file.remove(strImg)

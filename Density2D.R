#!/usr/bin/env Rscript
##Density2D.R
##
args = commandArgs(trailingOnly = TRUE)
if(length(args)<1) q()
if(!require(ggplot2,quietly=T) || !require(MASS,quietly=T) || !require(dplyr,quietly=T,warn.conflicts=F)) stop("ggplot2, dplyr and MASS are required R package!")

get_density = function(x, y, ...) {
  dens <- kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

strCSV = args[1]
strFun <- args[2]
minExpr <- as.numeric(args[3])
colMap <- switch(args[4],"B",magma="A",inferno="B",plasma="C",viridis="D")

expr <- read.csv(strCSV,row.names=1,as.is=T,check.names=F)
#minExpr <- apply(expr,2,min)
expr <- expr[apply(expr,1,function(x)return(sum(x>minExpr)))>0,]
if(nrow(expr)<50) stop("Less than 50 cells expression above minimal value in at least one of two genes!")
if(MASS::bandwidth.nrd(expr[,1])<=0 || MASS::bandwidth.nrd(expr[,2])<=0)
  stop(paste("The difference between the first and second quantile of the expression (min=",min(apply(expr,2,min)),") are NOT positive for at least one gene in selected cells! Try to increase the expression cutoff (",minExpr,").",sep=""))
D <- expr %>% dplyr::mutate(density=get_density(expr[,1],expr[,2],n=min(100,nrow(expr)/10)))

p <- ggplot(D, aes_string(x = colnames(D)[1], y = colnames(D)[2], color = 'density')) +
  geom_point(size = 0.4) + theme_classic() + ggtitle(paste("Density on",nrow(expr),"cells")) +
  geom_vline(xintercept = 0, color = "red", linetype = 2) + 
  geom_hline(yintercept = 0, color = "red", linetype = 2) + 
  theme(axis.text = element_text(face = "bold",size = 12)) +
  viridis::scale_color_viridis(option = colMap) +  
  scale_shape_identity() 

dpi <- 1
if(sum(strFun%in%c('png','jpeg','tiff'))>0) dpi <- 150
strImg <- gsub("csv$",strFun,strCSV)
f <- get(strFun)
f(strImg, width=3*dpi, height=3*dpi)
print(p)
a <- dev.off()
fig = base64enc::dataURI(file = strImg)
cat(gsub("data:;base64,","",fig))
a <- file.remove(strImg)


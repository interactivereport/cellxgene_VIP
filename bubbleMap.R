#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
if(length(args)<1) q()
libPath = tail(args,1)
if(nchar(libPath)>3){
  addPath <- unlist(strsplit(libPath,";"))
  addPath <- addPath[sapply(addPath,dir.exists)]
  .libPaths(c(addPath,.libPaths()))
}
require(tidyverse)
require(ggpubr)
options(bitmapType='cairo')

strCSV = args[1]# csv [gene,tag,log2fc,pval,qval]
strFun <- args[2]
fontsize <- as.numeric(args[3])
dpi <- as.numeric(args[4])
figScale <- as.numeric(args[5])#0.5

X <- read.csv(strCSV,check.names=F)

container = data.frame(tag = sort(factor(rep(levels(X$tag),nlevels(X$gene)),levels=unique(X$tag))),
                       gene = rep(levels(X$gene),nlevels(X$tag)))
#Make graph format
D <- container %>% left_join(X) %>%
  mutate(abslog2FC = abs(log2fc)) %>%
  mutate(abslog2FC.bin = cut(abslog2FC, breaks=c(-Inf,log2(1.2),log2(1.5),Inf),
                             labels = c("abs(log2FC)<=log2(1.2)","abs(log2FC)<=log2(1.5)","abs(log2FC)>log2(1.5)"))) %>%
  mutate(abslog2FC.bin = case_when(is.na(abslog2FC.bin) ~ "Not Measured",# & is.na(Warning)
                                   #is.na(abslog2FC.bin) & !is.na(Warning) ~ "glmmTMB Warning",
                                   TRUE ~ as.character(abslog2FC.bin) )) %>%
  mutate(FDR.bin = cut(qval, breaks=c(-Inf,0.01,0.05,0.1,1),
                       labels = c("FDR<=0.01","FDR<=0.05","FDR<=0.1","FDR>0.1"))) %>%
  mutate(FDR.bin = case_when(is.na(FDR.bin) ~ "Not Measured",#& is.na(Warning)
                             #is.na(FDR.bin) & !is.na(Warning) ~ "glmmTMB Warning",
                             TRUE ~ as.character(FDR.bin) )) %>%
  mutate(abslog2FC.bin = factor(abslog2FC.bin, levels = c("abs(log2FC)>log2(1.5)","abs(log2FC)<=log2(1.5)","abs(log2FC)<=log2(1.2)","Not Measured") )) %>% #,"glmmTMB Warning"
  mutate(log2FC.sign = case_when(sign(log2fc) == 1 ~ "Positive",
                                 sign(log2fc) == -1 ~ "Negative")) %>%
  mutate(log2FC.sign = case_when(is.na(log2FC.sign)  ~ "Not Measured",#& is.na(Warning)
                                 #is.na(log2FC.sign) & !is.na(Warning) ~ "glmmTMB Warning",
                                 TRUE ~ as.character(log2FC.sign) )) %>%
  mutate(FDR.bin = as.factor(FDR.bin),
         log2FC.sign = as.factor(log2FC.sign))

#Do graph output
g <- ggplot(D, aes(y=gene, x=tag, shape=log2FC.sign, fill=FDR.bin, size=abslog2FC.bin)) +
  scale_size_manual(name = "abs(log2FC) bin", values = c("abs(log2FC)>log2(1.5)"=18,"abs(log2FC)<=log2(1.5)"=12,"abs(log2FC)<=log2(1.2)"=6,"Not Measured" = 2)) +#
  scale_shape_manual(name= "log2FC sign",values = c("Positive" = 24, "Negative" = 25, "Not Measured" = 32 ), breaks = c("Positive","Negative")) +#
  scale_fill_manual(name = "FDR bin", values = c("FDR<=0.01" = "#FC8D59", "FDR<=0.05" = "#FED99D",
                                                 "FDR<=0.1" = "#DAE9C8", "FDR>0.1" = "#91BFDB", "Not Measured" = "transparent"),#
                    breaks = c("FDR<=0.01", "FDR<=0.05",
                               "FDR<=0.1", "FDR>0.1") ) +
  geom_point(shape=21, color="black") +
  geom_point(size=3, fill="gray50") +
  theme(panel.background = element_blank(),
        axis.title = element_text(size=rel(1.25)),
        axis.text.x = element_text(angle=90,hjust = 0.05,size=rel(1.25)),
        axis.text.x.top = element_text(vjust = 0.5),
        axis.text.y = element_text(size=rel(1.25)),
        panel.border = element_rect(color="black",fill=NA,size=1),
        legend.key = element_rect(fill=NA),
        plot.margin=unit(c(0.5, 11, 0.5, 0.5), "lines"),
        legend.justification = c(0,0),
        legend.position=c(1,0),
        text=element_text(size=fontsize)) +
  scale_x_discrete(position = "top") +
  labs(y = "Gene Symbol",x="") +
  guides(fill = guide_legend(override.aes = list(size=3)))

strImg <- gsub("csv$",strFun,strCSV)
f <- get(strFun)
if(sum(strFun%in%c('png','jpeg','tiff'))>0){
  f(strImg, width=5+round(nlevels(D$tag)*figScale,1), height=max(6,2+round(nlevels(D$gene)*figScale,1)),units='in',res=dpi)
}else{
  f(strImg, width=5+round(nlevels(D$tag)*figScale,1), height=max(6,2+round(nlevels(D$gene)*figScale,1)))
}
print(g)
a <- dev.off()
fig = base64enc::dataURI(file = strImg)
cat(gsub("data:;base64,","",fig))
a <- file.remove(strImg)

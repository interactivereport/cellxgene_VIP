%load_ext rpy2.ipython
%Rdevice svg

%%R -w 8 -h 8
library(Seurat, lib.loc = '/usr/lib64/R/library/Seurat_3.0.2/')
suppressPackageStartupMessages(library(SeuratData))
library(ggplot2)
library(patchwork)
data("pbmc3k.final")
pbmc3k.final$groups <- sample(c("group1", "group2"), size = ncol(pbmc3k.final), replace = TRUE)
features <- c("LYZ", "CCL5", "IL32", "PTPRCAP", "FCGR3A", "PF4")
#
p<-RidgePlot(pbmc3k.final, features = features, ncol = 2)
print(p)

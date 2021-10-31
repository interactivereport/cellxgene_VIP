#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
if(length(args)<1) q()
libPath = tail(args,1)
if(nchar(libPath)>3){
  addPath <- unlist(strsplit(libPath,";"))
  addPath <- addPath[sapply(addPath,dir.exists)]
  .libPaths(c(addPath,.libPaths()))
}
loadPackages <- function(){
  require(ggplot2)
  require(dplyr)
  require(RColorBrewer)
  require(glue)
  require(ggpubr)
  require(gridExtra)
}
suppressMessages(suppressWarnings(loadPackages()))
options(bitmapType='cairo')

violinPlot <- function(X,cutoff=0){
  #X[cells,(gene.name,grouping,subgrouping)]
  gene.name <- colnames(X)[1]
  colnames(X)[1] <- 'count'
  color.type <- fill.type <- grouping <- colnames(X)[2]
  paired.type <- "Unpaired Comparison"
  if(ncol(X)>2){
    color.type <- fill.type <- colnames(X)[3]
    paired.type <- "Paired Comparison"
  }
  X <- cbind(X,ordering=as.numeric(X[,grouping]))

  per.obj <- X %>% group_by(!!sym(grouping)) %>%
    mutate(n.per = sum(count>cutoff) / length(count),
           average.count = sum(count, na.rm=TRUE) / length(count),
           gene.id = gene.name) %>%
    dplyr::select(!!grouping, n.per, average.count, ordering, gene.id) %>% unique() %>% as.data.frame() %>%
    mutate(!!sym(grouping) := reorder(!!sym(grouping), ordering) ) %>% #keep the reorder
    mutate(n.per.bin = cut(n.per, breaks=c(-Inf,0.2,0.4,0.6,0.8,Inf),
                           labels = c("<20%","<40%","<60%","<80%",">80%"))
    ) #n.per will be used as a size scale, but continuous size scale is not informative, will bin them to a discrete scale

  # The ggplot becomes generalized

  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }

  fill.expand <- gg_color_hue( length( levels(as.factor(X[[fill.type]])) ) )

  #Point color needs additional levels based on selection to use alternative palette
  colors.expand <- colorRampPalette(brewer.pal(8, "Dark2"))( length( levels(as.factor(X[[color.type]])) ) )

  if(paired.type == "Unpaired Comparison"){

    general.violin <- ggplot(X, aes_string(x = grouping, y= "count" )) +
      geom_violin(aes_string(fill = fill.type), draw_quantiles =  0.5,  alpha=0.2) +
      geom_jitter(aes_string(colour = color.type), size = 0.2,position = position_jitter(0.2)) +
      scale_color_manual(values = colors.expand, name = glue("{color.type}.color")) +
      scale_fill_manual(values = fill.expand, name = glue("{fill.type}.fill") ) +
      stat_summary(
        fun=mean,fun.min=mean,fun.max=mean,
        #fun.y = mean, fun.ymin = mean, fun.ymax = mean,
        geom = "point",pch=8, color="red",size=2.25,
        position = position_dodge(width = .25)
      ) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, size=rel(1.25)),
            axis.text.y = element_text(size = rel(1.25)),
            axis.title = element_text(size = rel(1.25)),
            legend.text = element_text(size=rel(1.25)),
            legend.title = element_text(size=rel(1.25)),
            legend.margin = margin(r=100) #Not ideal but provides some right margin pad, may want to make generic legend label
      ) +
      labs(y = "log2(Normalized Counts + 1)"
           #title = paste0(gene.name," in ",cluster,". The horizontal line in the violin plot represents the median")
      ) +
      guides(colour = guide_legend(title=glue("{color.type} "), override.aes = list(size=5) ), fill = guide_legend(title=glue("{fill.type}")) )

  }

  if(paired.type == "Paired Comparison"){

    general.violin <- ggplot(X, aes_string(x = grouping, y= "count" )) +
      geom_violin(aes_string(fill = fill.type), draw_quantiles =  0.5,  alpha=0.2, position = position_dodge(0.75) ) +
      geom_point(aes_string(colour = color.type), size = 0.2, position = position_jitterdodge(seed=1, dodge.width = 0.75) ) +
      #geom_jitter(aes_string(colour = color.type), size = 0.2, position = position_jitter(0.2)) +
      scale_color_manual(values = colors.expand, name = glue("{color.type}.color")) +
      scale_fill_manual(values = fill.expand, name = glue("{fill.type}.fill") ) +
      stat_summary(
        aes_string(group=color.type),
        fun=mean,
        fun.min=mean,
        fun.max=mean,
        #fun.y = mean, fun.ymin = mean, fun.ymax = mean,
        geom = "point",pch=8, color="red",size=2.25,
        position = position_dodge(width = 0.75)
      ) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, size=rel(1.25)),
            axis.text.y = element_text(size = rel(1.25)),
            axis.title = element_text(size = rel(1.25)),
            legend.text = element_text(size=rel(1.25)),
            legend.title = element_text(size=rel(1.25)),
            legend.margin = margin(r=100) #Not ideal but provides some right margin pad, may want to make generic legend label
      ) +
      labs(y = "log2(Normalized Counts + 1)"
           #title = paste0(gene.name," in ",cluster,". The horizontal line in the violin plot represents the median")
      ) +
      guides(colour = guide_legend(title=glue("{color.type} "), override.aes = list(size=5) ), fill = guide_legend(title=glue("{fill.type}")) )

  }

  #Make a ggplot for the expression dots
  per.dots <- ggplot(per.obj, aes_string(x = grouping, y = "gene.id", color = "average.count", size = "n.per.bin") ) +
    geom_point() +
    scale_size_manual(name = "Percent \nExpr.\n(binned)", values = c("<20%"=2,"<40%"=6,"<60%"=10,"<80%"=14,">80%"=18)) +
    theme_bw() +
    labs(y="Gene Name", color="Avg.\nNorm. \nCount") +
    theme(legend.position = "right",
          legend.margin = margin(r=100),
          legend.text = element_text(size=rel(1.25)),
          legend.title = element_text(size=rel(1.25)),
          axis.text.x = element_text(angle = 90, size = rel(1.25)),
          axis.title.x = element_blank(),
          axis.line = element_blank(),
          axis.text.y = element_text(angle=90, hjust=0.5, vjust=0.5, size=rel(1.25)),
          axis.title.y = element_text(size=rel(1.25)),
          plot.margin = margin(t = 25)
    )

  #Combine the expression dots and the general violin plot
  annotate_figure(
    ggarrange( ggarrange(per.dots, general.violin, nrow=2, align="v", heights = c(1.5,4), legend="none"),
               ggarrange(NULL, ncol=1),
               widths=c(4,0.75),
               common.legend = TRUE,
               legend = "right",
               legend.grob = arrangeGrob(grobs = list(get_legend(per.dots),get_legend(general.violin)))
    ),
    top = text_grob(""),
    fig.lab = paste0(gene.name,". The horizontal line in the violin plot represents the median.\n      The red star represents the mean."), fig.lab.size=14
  )

}

strCSV = args[1]
expCut <- as.numeric(args[2])
strFun <- args[3]
fontsize <- as.numeric(args[4])
dpi <- as.numeric(args[5])
## process
mtable <- read.csv(strCSV,check.names=F)
## remove 0 on the violin plot
#mtable <- mtable[mtable[,1]>=expCut,]

g <- violinPlot(mtable,expCut)
## plot
grpN <- compN <- nlevels(mtable[,2])
if(ncol(mtable)>2){
  grpN <- length(unique(apply(mtable,1,function(x)return(paste(x[-1],collapse="::")))))
  compN <- nlevels(mtable[,3])
}
w <- max(8,2+round(2*grpN/3,1))
h <- max(8,compN/2)
strImg <- gsub("csv$",strFun,strCSV)
f <- get(strFun)
if(sum(strFun%in%c('png','jpeg','tiff'))>0){
  f(strImg, width=w, height=h,units='in',res=dpi)
}else{
  f(strImg, width=w, height=h)
}
print(g)
a <- dev.off()
fig = base64enc::dataURI(file = strImg)
cat(gsub("data:;base64,","",fig))
a <- file.remove(strImg)

#./violin.R /share/oyoung/violin.csv Gene1 pdf 10 300

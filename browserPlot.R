#!/usr/bin/env Rscript
load <- function(libPath){
    if(nchar(libPath)>3){
      addPath <- unlist(strsplit(libPath,";"))
      addPath <- addPath[sapply(addPath,dir.exists)]
      .libPaths(c(addPath,.libPaths()))
    }
    require(ggplot2)
    require(reshape2)
    require(ggforce)
    require(patchwork)
    require(rtracklayer)
}

main <- function(){
    ## process input ----
    args <- commandArgs(trailingOnly=TRUE)
    suppressMessages(suppressWarnings(load(tail(args,1))))
    strPath <- paste0(normalizePath(args[1]),"/")
    region <- args[2]
    bwList <- unlist(strsplit(args[3],","))
    extend.upstream <- as.numeric(args[4])
    extend.downstream <- as.numeric(args[5])
    strExp <- args[6] # have to be provided, used to save temp figures
    expCutoff <- as.numeric(args[7])
    strFun <- args[8]
    fontsize <- as.numeric(args[9])
    dpi <- as.numeric(args[10])

    ## obtain the region -----
    annotations <- NULL
    if(file.exists(paste0(strPath,"annotation.rds"))) annotations <- readRDS(paste0(strPath,"annotation.rds"))
    region <- customFindRegion(
        region = region,
        annotations = annotations,
        extend.upstream = extend.upstream,
        extend.downstream = extend.downstream
    )

    ## plot the bigwig -----    
    bwList <- bwList[bwList%in%list.files(strPath,"bw$")]
    strBW <- paste0(strPath,bwList)
    AllPlots <- suppressWarnings(suppressMessages(customBigwigTrack(strBW,region,fontsize)))
    #h <- rep(1,length(AllPlots))
    h <- length(strBW)
    w <- 6

    ## plot gene expression ----
    strCluster <- paste0(strPath,"/bw.cluster")
    expPlot <- customExpressionPlot(strExp,strCluster,bwList,expCutoff,fontsize)
    if(!is.null(expPlot)){
      w <- c(w,nlevels(expPlot$data$gene))
    }

    ## plot gene annotation -----
    if(!is.null(annotations)){
      p <- customAnnotationPlot(annotations,region,fontsize)#suppressWarnings(suppressMessages())
      if(!is.null(p)){
          AllPlots <- c(AllPlots,list(p))
          h <- c(h,0.7)
      }
    }
    strAnno <- paste0(strPath,"/annotation.rds")
    if(file.exists(strAnno)){

    }

    ## plot peaks -----
    strPeaks <- paste0(strPath,"/peaks.rds")
    if(file.exists(strAnno)){
        p <- suppressWarnings(suppressMessages(customPeakPlot(readRDS(strPeaks),region,fontsize)))
        if(!is.null(p)){
            AllPlots <- c(AllPlots,list(p))
            h <- c(h,0.4)
        }
    }

    ## plot links -----
    strLinks <- paste0(strPath,"/links.rds")
    if(file.exists(strLinks)){
        p <- suppressWarnings(suppressMessages(customLinkPlot(readRDS(strLinks),region,fontsize)))
        if(!is.null(p)){
            AllPlots <- c(AllPlots,list(p))
            h <- c(h,0.7)
        }
    }
    # save all plots -----
    strImg <- gsub("csv$",strFun,strExp)
    f <- get(strFun)
    if(sum(strFun%in%c('png','jpeg','tiff'))>0){
        f(strImg, width=sum(w), height=max(4,sum(h)),units='in',res=dpi)
    }else{
        f(strImg, width=sum(w), height=max(4,sum(h)))
    }
    suppressWarnings(suppressMessages(print(CombineTracks(plotlist=AllPlots,heights=h,
                                                          expression.plot=expPlot,widths=w))))
    a <- dev.off()
    fig = base64enc::dataURI(file = strImg)
    cat(gsub("data:;base64,","",fig))
    a <- file.remove(strImg)

}
customFindRegion <- function(region,annotations=NULL,sep = c("-", "-"),extend.upstream = 0,extend.downstream = 0) {
    if (!is(object = region, class2 = "GRanges")) {
        # first try to convert to coordinates, if not lookup gene
        region <- tryCatch(
            expr = suppressWarnings(
                expr = StringToGRanges(regions = region, sep = sep)
            ),
            error = function(x) {
                if(is.null(annotations)) stop("Annotation not found")
                region <- customLookupGeneCoords(annotations,region)
                return(region)
            }
        )
        if (is.null(x = region)) {
            stop("Gene/Region is not found")
        }
    }
    region <- suppressWarnings(expr = Extend(
        x = region,
        upstream = extend.upstream,
        downstream = extend.downstream
    )
    )
    return(region)
}
StringToGRanges <- function(regions, sep = c("-", "-"), ...) {
  ranges.df <- data.frame(ranges = regions)
  ranges.df <- tidyr::separate(
    data = ranges.df,
    col = "ranges",
    sep = paste0(sep[[1]], "|", sep[[2]]),
    into = c("chr", "start", "end")
  )
  granges <- makeGRangesFromDataFrame(df = ranges.df, ...)
  return(granges)
}
customLookupGeneCoords <- function(annotations, gene) {
    isgene <- annotations$gene_name == gene
    isgene <- !is.na(x = isgene) & isgene
    annot.sub <- annotations[isgene]
    if (length(x = annot.sub) == 0) {
        return(NULL)
    } else {
        gr <- GRanges(seqnames = as.character(x = GenomicRanges::seqnames(x = annot.sub))[[1]],
                      ranges = IRanges::IRanges(start = min(IRanges::start(x = annot.sub)),
                                       end = max(IRanges::end(x = annot.sub))))
        return(gr)
    }
}
customExpressionPlot <- function(strExp,strCluster,bwList,geneCutoff=-0.1,fontsize=9) {
  if(!file.exists(strExp)) return(NULL)
  Exp <- read.csv(strExp,row.names=1,as.is=T)
  clusterInfo <- read.table(strCluster,row.names=1,header=T,as.is=T,sep="\t")[bwList,,drop=F]
  grp <- intersect(colnames(Exp),colnames(clusterInfo))
  if(length(grp)==0) return(NULL)
  features <- colnames(Exp)[!colnames(Exp)%in%grp][-1]
  # get data
  df <- NULL
  for(oneBW in bwList){
    selC <- rep(T,nrow(Exp))
    for(oneG in grp){
        if(nchar(clusterInfo[oneBW,oneG])>0){
            selC <- selC & Exp[,oneG]==clusterInfo[oneBW,oneG]
        }
    }
    if(sum(selC)==0) stop(paste("No cell for",oneBW))
    oneD = setNames(melt(Exp[selC,features,drop=F],id.vars=NULL),c("gene","expression"))
    oneD$group <- oneBW
    df <- rbind(df,oneD)
  }
  df$group <- factor(df$group,levels=bwList)
  # set the color
  colors_all <- setNames(gg_color_hue(length(bwList)),bwList)

  p.list <- list()
  for (i in seq_along(along.with = features)) {
    df.use.all <- df[df$gene == features[i], ]
    df.use <- df.use.all[df.use.all$expression>geneCutoff,,drop=F]
    df.rate <- round(100*table(df.use$group)/table(df.use.all$group),1)
    df.rate <- data.frame(group=names(df.rate),
                          text=paste0(df.rate,"%"),
                          xpos=c(Inf),ypos=c(Inf))
    p <- ggplot(data = df.use, aes(x = expression, y = gene, fill = group)) +
      geom_violin(size = 1/4) +
      facet_wrap(~group, ncol = 1, strip.position = "right") +
      theme_classic() +
      scale_y_discrete(position = "top") +
      scale_x_continuous(position = "bottom", limits = c(0, NA)) +
      geom_text(data=df.rate,aes(xpos,ypos,label=text,color=group),hjust="right",vjust=1.1, size=fontsize-7)+
      scale_fill_manual(values = colors_all)+
      theme(
        axis.text.x.top = element_text(size=fontsize, face="bold", angle=0), # Gene names @ top
        axis.text.x.bottom = element_text(size=fontsize-1, angle=60, vjust=0.6),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none"
      )
    p.list[[i]] <- p
  }
  p <- patchwork::wrap_plots(p.list, ncol = length(x = p.list))
  return(p)
}
customBigwigTrack <- function(strBW,region,fontsize) {
    namePos <- data.frame(xpos=c(-Inf),ypos=c(Inf),text="",stringsAsFactors=F)
    nameBW <- setNames(sapply(strBW,function(x)return(trimws(gsub("\\.bw$","",basename(x))))),
                       strBW)
    bwCol <- setNames(gg_color_hue(length(strBW)),strBW)

    # each bigwig a plot -----
    if(F){
      AllPlots <- list()
      yLim <- 0
      yLab <- setNames(rep("",length(strBW)),strBW)
      yLab[round(length(yLab)/2)] <- "Coverage"
      for(one in strBW){
          namePos[1,"text"] <- nameBW[one]
          p <- BigwigTrack(region,one,y_label=yLab[one])+
              geom_text(data=namePos,aes(xpos,ypos,label=text), hjust="left", vjust="top", size=fontsize-4,color=bwCol[one])+
              geom_area(fill=bwCol[one])
          yLim <- max(c(yLim,layer_scales(p)$y$get_limits()))
          AllPlots <- c(AllPlots,list(p))
      }
      yLim <- round(ceiling(yLim/10))*10
      for(i in 1:length(AllPlots)){
          AllPlots[[i]] <- AllPlots[[i]]+
              scale_y_continuous(limits=c(0,yLim*1.05),breaks=c(0,yLim))+
              ##theme_classic()+
              theme(plot.margin=margin(0,0,0,0),
                    axis.text.y=element_text(face="italic"),
                    axis.title.y=element_text(size=fontsize))
      }
    }
    ## all bigwig in one plot -------
    if(T){
      coverages <- c()
      for(one in strBW){
        oneWig <- getBigwigCoverage(region,one)
        coverages <- rbind(coverages,cbind(oneWig,group=basename(one)))
      }
      bwCol <- setNames(gg_color_hue(length(strBW)),basename(strBW))
      xlabel <- paste0(as.character(x = GenomicRanges::seqnames(x = region)), " position (bp)")
      AllPlots <- list(customCoverageTrack(coverages,xlabel,bwCol,fontsize=fontsize))
    }
    return(AllPlots)
}
getBigwigCoverage <- function(region,bigwig,smooth = 200,max.downsample = 3000,downsample.rate = 0.1) {
    if (.Platform$OS.type == "windows") {
        message("BigwigTrack not supported on Windows")
        return(NULL)
    }
    if (!requireNamespace("rtracklayer", quietly = TRUE)) {
        message("Please install rtracklayer. http://www.bioconductor.org/packages/rtracklayer/")
        return(NULL)
    }
    if (!inherits(x = region, what = "GRanges")) {
        stop("region should be a GRanges object")
    }
    region_data <- rtracklayer::import(
        con = bigwig,
        which = region,
        as = "NumericList"
    )[[1]]
    if (!is.null(x = smooth)) {
        region_data <- RcppRoll::roll_mean(x = region_data, n = smooth, fill = 0L)
    }
    region_data <- data.frame(
        position = start(x = region):end(x = region),
        coverage = region_data,
        stringsAsFactors = FALSE
    )
    window.size = width(x = region)
    sampling <- min(max.downsample, window.size * downsample.rate)
    coverages <- dplyr::slice_sample(.data = region_data, n = sampling)
    return(coverages)
}
customCoverageTrack <- function(coverages,xlabel,colors_all=NULL,fontsize=12,region.highlight = NULL) {
  bwName <- levels(coverages$group)
  bwName <- data.frame(xpos=c(-Inf),ypos=c(Inf),text=paste0(" ",gsub("\\.bw$","",bwName)),group=bwName)
  yMax <- round(ceiling(max(coverages$coverage)/10))*10
  p <- ggplot(
    data = coverages,
    mapping = aes(x = position, y = coverage, fill = group)
  ) +
    geom_area(stat = "identity") +
    geom_hline(yintercept = 0, size = 0.1) +
    facet_wrap(facets = ~group, strip.position = "left", ncol = 1) +
    xlab(label = xlabel) +
    ylab(label = paste0("Normalized accessibility")) +
    geom_text(data=bwName,aes(xpos,ypos,label=text,color=group), hjust="left",vjust=1.1, size=fontsize-7)+
    theme_browser(legend = FALSE) +
    scale_y_continuous(limits=c(0,yMax*1.05),breaks=c(yMax))+
    theme(panel.spacing.y = unit(x = 0, units = "line"),
          axis.text.y=element_text(face="italic"),
          axis.title.y=element_text(size=fontsize,face="bold"),
          strip.text.y.left=element_blank(),
          strip.background = element_blank(),
          legend.position="none")
  if (!is.null(colors_all)) {
    p <- p + scale_fill_manual(values = colors_all) +
              scale_color_manual(values = colors_all)
  }
  if (!is.null(x = region.highlight)) {
    if (!inherits(x = region.highlight, what = "GRanges")) {
      warning("region.highlight must be a GRanges object")
    } else {
      md <- mcols(x = region.highlight)
      if ("color" %in% colnames(x = md)) {
        color.use <- md$color
      } else {
        color.use <- rep(x = "grey", length(x = region.highlight))
      }
      df <- data.frame(
        "start" = start(x = region.highlight),
        "end" = end(x = region.highlight),
        "color" = color.use
      )
      df$start <- ifelse(
        test = df$start < start.pos,
        yes = start.pos,
        no = df$start
      )
      df$end <- ifelse(
        test = df$end > end.pos,
        yes = end.pos,
        no = df$end
      )
      p <- p +
        geom_rect(
          data = df,
          inherit.aes = FALSE,
          aes_string(
            xmin = "start",
            xmax = "end",
            ymin = 0,
            ymax = ymax),
          fill = rep(x = df$color, length(x = unique(x = coverages$group))),
          color = "transparent",
          alpha = 0.2
        )
    }
  }
  return(p)
}

customAnnotationPlot <- function(annotation, region, fontsize) {
    if (is.null(x = annotation)) {
        return(NULL)
    }
    if (!inherits(x = region, what = "GRanges")) {
        region <- StringToGRanges(regions = region)
    }
    start.pos <- IRanges::start(x = region)
    end.pos <- IRanges::end(x = region)
    chromosome <- GenomicRanges::seqnames(x = region)

    # get names of genes that overlap region, then subset to include only those
    # genes. This avoids truncating the gene if it runs outside the region
    annotation.subset <- IRanges::subsetByOverlaps(x = annotation, ranges = region)
    genes.keep <- unique(x = annotation.subset$gene_name)
    annotation.subset <- annotation[
        fastmatch::fmatch(x = annotation$gene_name, table = genes.keep, nomatch = 0L) > 0L
        ]
    #print(annotation.subset)
    if (length(x = annotation.subset) == 0) {
        # make empty plot
        p <- ggplot(data = data.frame())
        y_limit <- c(0, 1)
    } else {
        annotation_df_list <- reformat_annotations(
            annotation = annotation.subset,
            start.pos = start.pos,
            end.pos = end.pos
        )
        p <- ggplot() +
            # exons
            geom_segment(
                data = annotation_df_list$exons,
                mapping = aes_string(
                    x = "start",
                    y = annotation_df_list$exons$dodge,
                    xend = "end",
                    yend = annotation_df_list$exons$dodge,
                    color = "strand"
                ),
                show.legend = FALSE,
                size = 5
            ) +
            # gene body
            geom_segment(
                data = annotation_df_list$labels,
                mapping = aes_string(
                    x = "start",
                    y = annotation_df_list$labels$dodge,
                    xend = "end",
                    yend = annotation_df_list$labels$dodge,
                    color = "strand"
                ),
                show.legend = FALSE,
                size = 1/2
            )
        if (nrow(x = annotation_df_list$plus) > 0) {
            # forward strand arrows
            p <- p + geom_segment(
                data = annotation_df_list$plus,
                mapping = aes_string(
                    x = "start",
                    y = annotation_df_list$plus$dodge,
                    xend = "end",
                    yend = annotation_df_list$plus$dodge,
                    color = "strand"
                ),
                arrow = arrow(
                    ends = "last",
                    type = "open",
                    angle = 45,
                    length = unit(x = 0.05, units = "inches")
                ),
                show.legend = FALSE,
                size = 1/2
            )
        }
        if (nrow(x = annotation_df_list$minus) > 0) {
            # reverse strand arrows
            p <- p + geom_segment(
                data = annotation_df_list$minus,
                mapping = aes_string(
                    x = "start",
                    y = annotation_df_list$minus$dodge,
                    xend = "end",
                    yend = annotation_df_list$minus$dodge,
                    color = "strand"
                ),
                arrow = arrow(
                    ends = "first",
                    type = "open",
                    angle = 45,
                    length = unit(x = 0.05, units = "inches")
                ),
                show.legend = FALSE,
                size = 1/2
            )
        }
        # label genes
        n_stack <- max(annotation_df_list$labels$dodge)
        annotation_df_list$labels$dodge <- annotation_df_list$labels$dodge + (n_stack * 0.2)
        p <- p + geom_text(
            data = annotation_df_list$labels,
            mapping = aes_string(x = "position", y = "dodge", label = "gene_name"),
            nudge_y = 0.3,
            size = 3
        )
        y_limit <- c(0.9, n_stack + (n_stack * 0.5))
    }
    p <- p +
        theme_classic() +
        ylab("Genes") +
        xlab(label = paste0(chromosome, " position (kb)")) +
        xlim(start.pos, end.pos) +
        ylim(y_limit) +
        theme(
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.title=element_text(size=fontsize, face="bold")
        ) +
        scale_color_manual(values = c("darkblue", "darkgreen"))
    return(p)
}
customPeakPlot <- function(peaks,region,fontsize,group.by = NULL,color = "dimgrey") {
    if (!inherits(x = region, what = "GRanges")) {
        region <- StringToGRanges(regions = region)
    }
    if (is.null(x = peaks)) {
        return(NULL)
    }
    # subset to covered range
    peak.intersect <- IRanges::subsetByOverlaps(x = peaks, ranges = region)
    peak.df <- as.data.frame(x = peak.intersect)
    start.pos <- IRanges::start(x = region)
    end.pos <- IRanges::end(x = region)
    chromosome <- GenomicRanges::seqnames(x = region)

    if (nrow(x = peak.df) > 0) {
        if (!is.null(x = group.by)) {
            if (!(group.by %in% colnames(x = peak.df))) {
                warning("Requested grouping variable not found")
                group.by <- NULL
            }
        }
        peak.df$start[peak.df$start < start.pos] <- start.pos
        peak.df$end[peak.df$end > end.pos] <- end.pos
        peak.plot <- ggplot(
            data = peak.df,
            aes_string(color = SetIfNull(x = group.by, y = "color"))
        ) +
            geom_segment(aes(x = start, y = 0, xend = end, yend = 0),
                         size = 2,
                         data = peak.df)
    } else {
        # no peaks present in region, make empty panel
        peak.plot <- ggplot(data = peak.df)
    }
    peak.plot <- peak.plot + theme_classic() +
        ylab(label = "Peaks") +
        theme(axis.ticks.y = element_blank(),
              axis.text.y = element_blank(),
              axis.title=element_text(size=fontsize, face="bold")) +
        xlab(label = paste0(chromosome, " position (bp)")) +
        xlim(c(start.pos, end.pos))
    if (is.null(x = group.by)) {
        # remove legend, change color
        peak.plot <- peak.plot +
            scale_color_manual(values = color) +
            theme(legend.position = "none")
    }
    return(peak.plot)
}
customLinkPlot <- function(links,region, fontsize, object=NULL, min.cutoff = 0) {
    if (!inherits(x = region, what = "GRanges")) {
        region <- StringToGRanges(regions = region)
    }
    chromosome <- GenomicRanges::seqnames(x = region)

    # extract link information
    if(!is.null(object)) links <- Links(object = object)
    if(is.null(links)) stop("Either provide seruat object contains links or Link object itself!")
    # if links not set, return NULL
    if (length(x = links) == 0) {
        return(NULL)
    }

    # subset to those in region
    links.keep <- IRanges::subsetByOverlaps(x = links, ranges = region)

    # filter out links below threshold
    link.df <- as.data.frame(x = links.keep)
    link.df <- link.df[abs(x = link.df$score) > min.cutoff, ]

    # remove links outside region
    link.df <- link.df[link.df$start >= IRanges::start(x = region) & link.df$end <= IRanges::end(x = region), ]

    # plot
    if (nrow(x = link.df) > 0) {
        # convert to format for geom_bezier
        link.df$group <- seq_len(length.out = nrow(x = link.df))
        df <- data.frame(
            x = c(link.df$start,
                  (link.df$start + link.df$end) / 2,
                  link.df$end),
            y = c(rep(x = 0, nrow(x = link.df)),
                  rep(x = -1, nrow(x = link.df)),
                  rep(x = 0, nrow(x = link.df))),
            group = rep(x = link.df$group, 3),
            score = rep(link.df$score, 3)
        )
        p <- ggplot(data = df) +
            geom_bezier(
                mapping = aes_string(x = "x", y = "y", group = "group", color = "score")
            ) +
            geom_hline(yintercept = 0, color = 'grey') +
            scale_color_gradient2(low = "red", mid = "grey", high = "blue")
    } else {
        p <- ggplot(data = link.df)
    }
    p <- p +
        theme_classic() +
        theme(axis.ticks.y = element_blank(),
              axis.text.y = element_blank(),
              axis.title=element_text(size=fontsize, face="bold")) +
        ylab("Links") +
        xlab(label = paste0(chromosome, " position (bp)")) +
        xlim(c(IRanges::start(x = region), IRanges::end(x = region)))
    return(p)
}
BigwigTrack <- function(region,bigwig,smooth = 200,type = "coverage",y_label = "Score",max.downsample = 3000,downsample.rate = 0.1) {
    possible_types <- c("line", "heatmap", "coverage")
    if (!(type %in% possible_types)) {
        stop(
            "Invalid type requested. Choose ",
            paste(possible_types, collapse = ", ")
        )
    }
    if (.Platform$OS.type == "windows") {
        message("BigwigTrack not supported on Windows")
        return(NULL)
    }
    if (!requireNamespace("rtracklayer", quietly = TRUE)) {
        message("Please install rtracklayer. http://www.bioconductor.org/packages/rtracklayer/")
        return(NULL)
    }
    if (!inherits(x = region, what = "GRanges")) {
        stop("region should be a GRanges object")
    }
    region_data <- rtracklayer::import(
        con = bigwig,
        which = region,
        as = "NumericList"
    )[[1]]
    if (!is.null(x = smooth)) {
        region_data <- RcppRoll::roll_mean(x = region_data, n = smooth, fill = 0L)
    }
    region_data <- data.frame(
        position = start(x = region):end(x = region),
        score = region_data,
        stringsAsFactors = FALSE
    )
    window.size = width(x = region)
    sampling <- min(max.downsample, window.size * downsample.rate)
    coverages <- dplyr::slice_sample(.data = region_data, n = sampling)
    p <- ggplot(
        data = coverages,
        mapping = aes_string(x = "position", y = "score")
    )
    if (type == "line") {
        p <- p + geom_line()
    } else if (type == "heatmap") {
        # different downsampling needed for heatmap
        # cut into n bins and average within each bin
        region_data$bin <- floor(x = region_data$position / smooth)
        region_data <- dplyr::group_by(region_data, bin)
        region_data <- dplyr::mutate(region_data, score = mean(x = score))
        region_data <- dplyr::ungroup(region_data)
        region_data <- unique(x = region_data[, c("bin", "score")])
        p <- ggplot(
            data = region_data,
            mapping = aes_string(x = "bin", y = 1, fill = "score")
        ) + geom_tile() + scale_fill_viridis_c()
    } else if (type == "coverage") {
        p <- p + geom_area()
    }
    chromosome <- as.character(x = GenomicRanges::seqnames(x = region))
    p <- p + theme_browser() +
        xlab(label = paste0(chromosome, " position (bp)")) +
        ylab(label = y_label)
    return(p)
}
CombineTracks <- function(plotlist,expression.plot = NULL,heights = NULL,widths = NULL) {
    # remove any that are NULL
    nullplots <- sapply(X = plotlist, FUN = is.null)
    plotlist <- plotlist[!nullplots]
    heights <- heights[!nullplots]

    if (length(x = plotlist) == 1) {
        return(plotlist[[1]])
    }

    # remove x-axis from all but last plot
    for (i in 1:(length(x = plotlist) - 1)) {
        plotlist[[i]] <- plotlist[[i]] + theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.line.x.bottom = element_blank(),
            axis.ticks.x.bottom = element_blank()
        )
    }

    # combine plots
    if (is.null(x = heights)) {
        # set height of first element to 10x more than other elements
        n.plots <- length(x = plotlist)
        heights <- c(8, rep(1, n.plots - 1))
    } else {
        if (length(x = heights) != length(x = plotlist)) {
            stop("Relative height must be supplied for each plot")
        }
    }
    if (!is.null(x = expression.plot)) {
        # align expression plot with the first element in plot list
        p <- (plotlist[[1]] + expression.plot) +
            patchwork::plot_layout(widths = widths)

        n <- length(x = plotlist)
        if(n>1){
          heights.2 <- heights[2:n]
          p2 <- patchwork::wrap_plots(plotlist[2:n], ncol = 1, heights = heights.2)

          p <- p + p2 + patchwork::guide_area() + patchwork::plot_layout(
              ncol = 2, heights = c(heights[[1]], sum(heights.2)),
              guides = "collect")
        }
    } else {
        p <- patchwork::wrap_plots(plotlist, ncol = 1, heights = heights)
    }
    return(p)
}

gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}
split_body <- function(df, width = 1000) {
    wd <- df$end - df$start
    nbreak <- wd / width
    if (nbreak > 1) {
        steps <- 0:(nbreak)
        starts <- (width * steps) + df$start
        starts[starts > df$end] <- NULL
    } else {
        starts <- df$end
    }
    breaks <- data.frame(
        seqnames = df$seqnames[[1]],
        start = starts,
        end = starts + 1,
        strand = df$strand[[1]],
        gene_name = df$gene_name[[1]],
        gene_biotype = df$gene_biotype[[1]],
        type = "arrow"
    )
    return(breaks)
}
theme_browser <- function(..., legend = TRUE) {
    browser.theme <- theme_classic() +
        theme(
            axis.text.y = element_blank(),
            strip.background = element_blank(),
            strip.text.y.left = element_text(angle = 0)
        )
    if (!legend) {
        browser.theme <- browser.theme +
            theme(
                legend.position = "none"
            )
    }
    return(browser.theme)
}
Extend <- function(x,upstream = 0,downstream = 0,from.midpoint = FALSE) {
    if (any(GenomicRanges::strand(x = x) == "*")) {
        warning("'*' ranges were treated as '+'")
    }
    on_plus <- GenomicRanges::strand(x = x) == "+" | GenomicRanges::strand(x = x) == "*"
    if (from.midpoint) {
        midpoints <- GenomicRanges::start(x = x) + (GenomicRanges::width(x = x) / 2)
        new_start <- midpoints - ifelse(
            test = on_plus, yes = upstream, no = downstream
        )
        new_end <- midpoints + ifelse(
            test = on_plus, yes = downstream, no = upstream
        )
    } else {
        new_start <- GenomicRanges::start(x = x) - ifelse(
            test = on_plus, yes = upstream, no = downstream
        )
        new_end <- end(x = x) + ifelse(
            test = on_plus, yes = downstream, no = upstream
        )
    }
    IRanges::ranges(x = x) <- IRanges::IRanges(start = new_start, end = new_end)
    x <- GenomicRanges::trim(x = x)
    return(x)
}
SetIfNull <- function(x, y) {
    if (is.null(x = x)) {
        return(y)
    } else {
        return(x)
    }
}
record_overlapping <- function(annotation, min.gapwidth = 1000) {
    # convert back to granges
    annotation.stash <- annotation
    annotation$strand <- "*"
    gr <- GenomicRanges::makeGRangesFromDataFrame(
        df = annotation[annotation$type == "body", ], keep.extra.columns = TRUE
    )
    # work out which ranges overlap
    collapsed <- GenomicRanges::reduce(
        x = gr, with.revmap = TRUE, min.gapwidth = min.gapwidth
    )$revmap
    idx <- seq_along(gr)
    for (i in seq_along(collapsed)) {
        mrg <- collapsed[[i]]
        for (j in seq_along(mrg)) {
            idx[[mrg[[j]]]] <- j
        }
    }
    names(x = idx) <- gr$gene_name
    return(idx)
}
reformat_annotations <- function(annotation,start.pos,end.pos) {
    annotation <- annotation[annotation$type == "exon"]
    exons <- as.data.frame(x = annotation,row.names=NULL)
    annotation <- split(
        x = annotation,
        f = annotation$gene_name
    )
    annotation <- lapply(X = annotation, FUN = as.data.frame,row.names=NULL)

    # add gene total start / end
    gene_bodies <- list()
    for (i in seq_along(annotation)) {
        df <- data.frame(
            seqnames = annotation[[i]]$seqnames[[1]],
            start = min(annotation[[i]]$start),
            end = max(annotation[[i]]$end),
            strand = annotation[[i]]$strand[[1]],
            gene_name = annotation[[i]]$gene_name[[1]],
            gene_biotype = annotation[[i]]$gene_biotype[[1]],
            type = "body"
        )
        # trim any that extend beyond region
        df$start <- ifelse(
            test = df$start < start.pos,
            yes = start.pos,
            no = df$start
        )
        df$end <- ifelse(
            test = df$end > end.pos,
            yes = end.pos,
            no = df$end
        )
        breaks <- split_body(df = df)
        df <- rbind(df, breaks)
        gene_bodies[[i]] <- df
    }
    gene_bodies <- do.call(what = rbind, args = gene_bodies)

    # record if genes overlap
    overlap_idx <- record_overlapping(annotation = gene_bodies, min.gapwidth = 1000)
    gene_bodies$dodge <- overlap_idx[gene_bodies$gene_name]
    exons$dodge <- overlap_idx[exons$gene_name]

    label_df <- gene_bodies[gene_bodies$type == "body", ]
    label_df$width <- label_df$end - label_df$start
    label_df$position <- label_df$start + (label_df$width / 2)

    onplus <- gene_bodies[gene_bodies$strand %in% c("*", "+"), ]
    onminus <- gene_bodies[gene_bodies$strand == "-", ]

    return(
        list(
            "labels" = label_df,
            "exons" = exons,
            "plus" = onplus,
            "minus" = onminus
        )
    )
}

main()

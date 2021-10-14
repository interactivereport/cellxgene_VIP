#!/usr/bin/env Rscript
load <- function(libPath){
    if(nchar(libPath)>3){
      addPath <- unlist(strsplit(libPath,";"))
      addPath <- addPath[sapply(addPath,dir.exists)]
      .libPaths(c(addPath,.libPaths()))
    }
    require(Signac)
    require(ggplot2)
    require(ggforce)
    require(rtracklayer)
}

main <- function(){
    ## process input ----
    args <- commandArgs(trailingOnly=TRUE)
    suppressMessages(suppressWarnings(load(tail(args,1))))
    strPath <- args[1]
    region <- args[2]
    extend.upstream <- as.numeric(args[3])
    extend.downstream <- as.numeric(args[4])

    strExp <- args[5] # have to be provided, used to save temp figures
    strFun <- args[6]
    fontsize <- as.numeric(args[7])
    dpi <- as.numeric(args[8])

    yLabSize <- 9
    ## obtain the region -----
    annotations <- NULL
    if(file.exists(paste0(strPath,"/annotation.rds"))) annotations <- readRDS(paste0(strPath,"/annotation.rds"))
    region <- customFindRegion(
        region = region,
        annotations = annotations,
        extend.upstream = extend.upstream,
        extend.downstream = extend.downstream
    )

    ## plot the bigwig -----
    strBW <- list.files(strPath,"bw$",full.names=T)
    AllPlots <- suppressWarnings(suppressMessages(customBigwigTrack(strBW,region,fontsize,yLabSize)))
    h <- rep(1,length(AllPlots))

    ## plot gene annotation -----
    strAnno <- paste0(strPath,"/annotation.rds")
    if(file.exists(strAnno)){
        p <- suppressWarnings(suppressMessages(customAnnotationPlot(readRDS(strAnno),region,yLabSize)))
        if(!is.null(p)){
            AllPlots <- c(AllPlots,list(p))
            h <- c(h,0.4)
        }
    }

    ## plot peaks -----
    strPeaks <- paste0(strPath,"/peaks.rds")
    if(file.exists(strAnno)){
        p <- suppressWarnings(suppressMessages(customPeakPlot(readRDS(strPeaks),region,yLabSize)))
        if(!is.null(p)){
            AllPlots <- c(AllPlots,list(p))
            h <- c(h,0.2)
        }
    }

    ## plot links -----
    strLinks <- paste0(strPath,"/links.rds")
    if(file.exists(strLinks)){
        p <- suppressWarnings(suppressMessages(customLinkPlot(readRDS(strLinks),region,yLabSize)))
        if(!is.null(p)){
            AllPlots <- c(AllPlots,list(p))
            h <- c(h,0.7)
        }
    }
    # save all plots -----
    strImg <- gsub("csv$",strFun,strExp)
    f <- get(strFun)
    if(sum(strFun%in%c('png','jpeg','tiff'))>0){
        f(strImg, width=8, height=max(4,sum(h)),units='in',res=dpi)
    }else{
        f(strImg, width=8, height=max(4,sum(h)))
    }
    suppressWarnings(suppressMessages(print(CombineTracks(plotlist=AllPlots,heights=h))))
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
            stop("Gene not found")
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
customLookupGeneCoords <- function(annotations, gene) {
    isgene <- annotations$gene_name == gene
    isgene <- !is.na(x = isgene) & isgene
    annot.sub <- annotations[isgene]
    if (length(x = annot.sub) == 0) {
        return(NULL)
    } else {
        gr <- GRanges(seqnames = as.character(x = seqnames(x = annot.sub))[[1]],
                      ranges = IRanges::IRanges(start = min(IRanges::start(x = annot.sub)),
                                       end = max(IRanges::end(x = annot.sub))))
        return(gr)
    }
}
customBigwigTrack <- function(strBW,region,fontsize,yLabSize) {
    namePos <- data.frame(xpos=c(-Inf),ypos=c(Inf),text="",stringsAsFactors=F)
    nameBW <- setNames(sapply(strBW,function(x)return(trimws(gsub("\\.bw$","",basename(x))))),
                       strBW)
    bwCol <- setNames(gg_color_hue(length(strBW)),strBW)

    AllPlots <- list()
    yLim <- 0
    yLab <- setNames(rep("",length(strBW)),strBW)
    yLab[round(length(yLab)/2)] <- "Coverage"
    for(one in strBW){
        namePos[1,"text"] <- nameBW[one]
        p <- BigwigTrack(region,one,y_label=yLab[one])+
            geom_area(fill=bwCol[one])+
            geom_text(data=namePos,aes(xpos,ypos,label=text),
                      hjust="left",vjust="top",
                      size=fontsize,color=bwCol[one])
        yLim <- max(c(yLim,layer_scales(p)$y$get_limits()))
        AllPlots <- c(AllPlots,list(p))
    }
    yLim <- round(ceiling(yLim/10))*10
    saveRDS(AllPlots,file="bw.rds")
    for(i in 1:length(AllPlots)){
        AllPlots[[i]] <- AllPlots[[i]]+
            scale_y_continuous(limits=c(0,yLim*1.05),breaks=c(0,yLim))+
            ##theme_classic()+
            theme(plot.margin=margin(0,0,0,0,"inches"),
                  axis.text.y=element_text(face="italic"),
                  axis.title.y=element_text(size=yLabSize))
    }
    return(AllPlots)
}
customAnnotationPlot <- function(annotation, region, yLabSize) {
    if (is.null(x = annotation)) {
        return(NULL)
    }
    if (!inherits(x = region, what = "GRanges")) {
        region <- StringToGRanges(regions = region)
    }
    start.pos <- IRanges::start(x = region)
    end.pos <- IRanges::end(x = region)
    chromosome <- seqnames(x = region)

    # get names of genes that overlap region, then subset to include only those
    # genes. This avoids truncating the gene if it runs outside the region
    annotation.subset <- subsetByOverlaps(x = annotation, ranges = region)
    genes.keep <- unique(x = annotation.subset$gene_name)
    annotation.subset <- annotation[
        fastmatch::fmatch(x = annotation$gene_name, table = genes.keep, nomatch = 0L) > 0L
        ]

    if (length(x = annotation.subset) == 0) {
        # make empty plot
        p <- ggplot(data = data.frame())
        y_limit <- c(0, 1)
    } else {
        annotation_df_list <- Signac:::reformat_annotations(
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
            axis.title=element_text(size=yLabSize)
        ) +
        scale_color_manual(values = c("darkblue", "darkgreen"))
    return(p)
}
customPeakPlot <- function(peaks,region,yLabSize,group.by = NULL,color = "dimgrey") {
    if (!inherits(x = region, what = "GRanges")) {
        region <- StringToGRanges(regions = region)
    }
    if (is.null(x = peaks)) {
        return(NULL)
    }
    # subset to covered range
    peak.intersect <- subsetByOverlaps(x = peaks, ranges = region)
    peak.df <- as.data.frame(x = peak.intersect)
    start.pos <- IRanges::start(x = region)
    end.pos <- IRanges::end(x = region)
    chromosome <- seqnames(x = region)

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
            aes_string(color = Signac:::SetIfNull(x = group.by, y = "color"))
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
              axis.title=element_text(size=yLabSize)) +
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
customLinkPlot <- function(links,region, yLabSize, object=NULL, min.cutoff = 0) {
    if (!inherits(x = region, what = "GRanges")) {
        region <- StringToGRanges(regions = region)
    }
    chromosome <- seqnames(x = region)

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
              axis.title.y=element_text(size=yLabSize)) +
        ylab("Links") +
        xlab(label = paste0(chromosome, " position (bp)")) +
        xlim(c(IRanges::start(x = region), IRanges::end(x = region)))
    return(p)
}
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}
main()

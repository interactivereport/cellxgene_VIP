#------- Load needed packages ---------------
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(glue))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(edgeR))
suppressMessages(library(DESeq2))
suppressMessages(library(apeglm))
suppressMessages(library(glmmTMB))
suppressMessages(library(BiocParallel))
suppressMessages(library(Seurat))
suppressMessages(library(MAST))
suppressMessages(library(R6))
suppressMessages(library(Matrix))
suppressMessages(library(slam))
suppressMessages(library(foreach))
suppressMessages(library(doMC))
suppressMessages(library(biomaRt))
suppressMessages(library(emmeans))
suppressMessages(library(nebula))

# ---------------------
# make helper functions
# ---------------------

my_CV <- function(x) sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)*100
my_FC_row <- function(arow, Gone, Gtwo) mean(arow[Gone], na.rm = TRUE)/mean(arow[Gtwo], na.rm = TRUE)
my_FC <- function(some_data,Gone,Gtwo) apply(some_data,1,function(x) my_FC_row(x,Gone,Gtwo))

#------- Function added for sparse matrix calculation ---------------

apply_sparseMatrix <- function(X, MARGIN, FUN, ...) {
    require(slam)
    stopifnot(MARGIN %in% c(1, 2))
    X2 <- as.simple_triplet_matrix(X) # X could be a dense matrix or sparse matrix (dgCMatrix/dgTMatrix).
    if (MARGIN == 1) {
        res <- rowapply_simple_triplet_matrix(X2, FUN, ...)
    } else {
        res <- colapply_simple_triplet_matrix(X2, FUN, ...)
    }
    return(res)
}

sweep_sparseMatrix <- function(x, margin, stats, fun = "*") {
    stopifnot(margin %in% c(1, 2))
    stopifnot(is.vector(stats))
    f <- match.fun(fun)
    x <- as(x, "dgTMatrix")
    if (margin == 1) {
        stopifnot(nrow(x) == length(stats))
        idx <- x@i + 1
    } else {
        stopifnot(ncol(x) == length(stats))
        idx <- x@j + 1
    }
    x@x <- f(x@x, stats[idx])
    return(x)
}

# R6 Class
# -------------------------------------------------------------------
BiostatsSingleCell =
    R6Class("BiostatsSingleCell",
            private = list(
                counts = NULL,
                norm_counts = NULL,
                cells = NULL,
                genes = NULL,
                meta_data = NULL,
                cluster_col = NULL,
                sampleId_col = NULL,
                treatment_col = NULL,
                MT_rows = NULL,
                lib_sizes = NULL,
                av_lib = NULL,
                pseudoBulk_counts = NULL,
                pseudoBulk_meta = NULL,
                norm_pseudoBulk_counts = NULL,
                mode = NULL,
                de_results = list(),
                filter_info = list(),
                de_method = NULL
            ),
            public = list(
                # --------------------------------
                # initialize function...
                #
                # option 1:
                # Read single cell data into memory from count_file (rows are genes, columns are cells,
                #                                                    gene names in column 1, cell names in column names)
                # Read meta data (sample id, treatment, and cluster information) from meta_file,
                # and column names given by sampleId_col, treatment_col, cluster_col
                #
                # option 2:
                # Initialize with count table (sparse matrix, rows are genes, columns are cells,
                #                              gene names in row names, cell names in column names)
                # and sample id, treatment, and cluster information in meta_data data frame
                # (column names given by sampleId_col, treatment_col, cluster_col)
                #
                # option 3:
                # Read single cell data (as a Seurat object) into memory from a rds_file
                # assuming there is a meta.data slot that includes sample id, treatment, and cluster information
                # and column names given by sampleId_col, treatment_col, cluster_col
                # ---------------------------------
                initialize = function(count_file = NULL, meta_file = NULL, sampleId_col = NULL, treatment_col = NULL,
                                      cluster_col = NULL, count_data = NULL, meta_data = NULL, rds_file = NULL, mode = NULL)
                {
                    method = 0
                    if (!is.null(count_file) & !is.null(meta_file) & !is.null(sampleId_col) & !is.null(treatment_col) &
                        !is.null(cluster_col) & is.null(count_data) & is.null(meta_data) & is.null(rds_file)) method = 1
                    if (is.null(count_file) & is.null(meta_file) & !is.null(sampleId_col) & !is.null(treatment_col) &
                        !is.null(cluster_col) & !is.null(count_data) & !is.null(meta_data) & is.null(rds_file)) method = 2
                    if (is.null(count_file) & is.null(meta_file) & !is.null(sampleId_col) & !is.null(treatment_col) &
                        !is.null(cluster_col) & is.null(count_data) & is.null(meta_data) & !is.null(rds_file)) method = 3
                    
                    stopifnot(method == 1 | method == 2 | method == 3)
                    
                    if (method == 1)
                    {
                        stopifnot(is.character(count_file))
                        stopifnot(is.character(meta_file))
                        stopifnot(is.character(sampleId_col))
                        stopifnot(is.character(treatment_col))
                        stopifnot(is.character(cluster_col))
                    }
                    if (method == 2)
                    {
                        stopifnot(inherits(count_data, "Matrix"))
                        stopifnot(is.data.frame(meta_data))
                        stopifnot(is.character(sampleId_col))
                        stopifnot(is.character(treatment_col))
                        stopifnot(is.character(cluster_col))
                    }
                    if (method == 3)
                    {
                        stopifnot(is.character(rds_file))
                        stopifnot(is.character(sampleId_col))
                        stopifnot(is.character(treatment_col))
                        stopifnot(is.character(cluster_col))
                    }
                    
                    if (method == 1) {
                        my_data = fread(count_file, data.table = FALSE)
                        rownames(my_data) <- my_data[[1]]
                        my_data[, 1] <- NULL
                        private$counts = as(as.matrix(my_data), "dgCMatrix")
                    }
                    if (method == 2) private$counts = count_data
                    if (method == 3) {
                        my_data = readRDS(rds_file)
                        private$counts = my_data@assays$RNA@counts
                    }
                    private$cells = colnames(private$counts)
                    private$genes = rownames(private$counts)
                    
                    cat("Dimensions of counts data: \n")
                    print(dim(private$counts))
                    cat("\n")
                    
                    # load the meta data
                    if (method == 1)
                    {
                        private$meta_data = fread(meta_file)
                    }
                    if (method == 2)
                    {
                        private$meta_data = meta_data
                    }
                    if (method == 3)
                    {
                        stopifnot(.hasSlot(my_data, "meta.data"))
                        if(length(my_data@meta.data$cell)>0){ #DHmod: if barcode data is recorded in "cell" column of the meta data,  just use meta data as is.
                            private$meta_data=my_data@meta.data
                        }else{ #if not, make a "cell" column with barcode data in the meta data #DHmod
                            private$meta_data = my_data@meta.data %>% tibble::rownames_to_column(var = "cell") # assuming there is a meta.data slot!!!
                        }
                    }
                    
                    private$cluster_col = cluster_col
                    private$treatment_col = treatment_col
                    private$sampleId_col = sampleId_col
                    
                    if(!("cell" %in% names(private$meta_data))) stop("Please make sure there is a 'cell' column in the meta data!")
                    alignment = all(private$cells == private$meta_data[,"cell"])
                    stopifnot(alignment)
                    
                    cat("Cluster summary...\n")
                    print(ftable(private$meta_data[, private$sampleId_col], private$meta_data[, private$cluster_col]))
                    cat("\n")
                    
                    # extract MT rows
                    private$MT_rows = grep("^MT-", private$genes, ignore.case = TRUE)
                    cat("Mitochondrial genes found: \n")
                    print(private$genes[private$MT_rows])
                    cat("\n")
                    
                    # normalize by Dann's method
                    private$lib_sizes = Matrix::colSums(private$counts, na.rm = TRUE)
                    private$av_lib = mean(private$lib_sizes)
                    
                    private$norm_counts =  sweep_sparseMatrix(private$counts, 2, private$lib_sizes/private$av_lib, "/")
                    
                    # get pseudo-bulk counts
                    # note: using dense matrix for pseudobulk counts
                    if (is.null(mode) || mode == "group") {
                        private$mode = "group"
                        #print(private$sampleId_col)
                        #print(seq_along(private$meta_data[, private$sampleId_col]))
                        #print(length(private$meta_data[, private$sampleId_col]))
                        by_set <- split(seq_along(private$meta_data[, private$sampleId_col]), private$meta_data[, private$sampleId_col])
                        pseudoBulk_counts <- sapply(names(by_set), function(x) {
                            idx <- by_set[[x]]
                            pseudo.bulk.counts <- Matrix::rowSums(private$counts[, idx, drop = FALSE],  sum, na.rm = TRUE)
                        })
                    } else {
                        private$mode = "cluster"
                        by_set <- split(seq_along(private$meta_data[, private$sampleId_col]), private$meta_data[, private$sampleId_col])
                        
                        pseudoBulk_counts_ref <- sapply(names(by_set), function(x) {
                            idx <- intersect(by_set[[x]], which(private$meta_data$.GrouP == "ref"))
                            pseudo.bulk.counts <- Matrix::rowSums(private$counts[, idx, drop = FALSE],  sum, na.rm = TRUE)
                        })
                        
                        pseudoBulk_counts_alt <- sapply(names(by_set), function(x) {
                            idx <- intersect(by_set[[x]], which(private$meta_data$.GrouP == "alt"))
                            pseudo.bulk.counts <- Matrix::rowSums(private$counts[, idx, drop = FALSE],  sum, na.rm = TRUE)
                        })
                        
                        pseudoBulk_counts <- merge(pseudoBulk_counts_ref, pseudoBulk_counts_alt, by = 0, suffixes = c(".ref", ".alt")) %>%
                            column_to_rownames(var = "Row.names")
                        
                        pseudoBulk_meta <- data.frame(PseudoBulk.Id = colnames(pseudoBulk_counts), stringsAsFactors = FALSE) %>%
                            dplyr::mutate(id_temp = gsub(".alt$", "", gsub('.ref$', "", PseudoBulk.Id)),
                                          .GrouP = "ref")
                        pseudoBulk_meta$.GrouP[grep('alt', pseudoBulk_meta$PseudoBulk.Id)] <- 'alt'
                        
                        names(pseudoBulk_meta)[2] <- private$sampleId_col
                        
                        private$pseudoBulk_meta <- pseudoBulk_meta
                        
                    }
                    
                    private$pseudoBulk_counts <- pseudoBulk_counts
                    cat("Dimensions of pseudo-bulk counts data: \n")
                    print(dim(private$pseudoBulk_counts))
                },
                # -----------------------------------
                # QC_data()
                #
                # Makes QC plots for scRNAseq datafile
                # saves results to pdf out_file
                # -----------------------------------
                make_QCplots = function(out_file)
                {
                    pdf(out_file, width = 9, height = 7)
                    # QC plots
                    sce <- SingleCellExperiment(assays = list(counts=private$counts),
                                                colData = data.frame(cluster = private$meta_data[, private$cluster_col],
                                                                     sample = private$meta_data[, private$sampleId_col],
                                                                     treatment = private$meta_data[, private$treatment_col]))
                    
                    if (length(private$MT_rows)>0)
                    {
                        sce <- calculateQCMetrics(sce, feature_controls=list(Mt=private$MT_rows))
                        
                        sce_mt <- SingleCellExperiment(assays = list(counts=private$counts[-private$MT_rows,]),
                                                       colData = data.frame(cluster = private$meta_data[, private$cluster_col],
                                                                            sample = private$meta_data[, private$sampleId_col],
                                                                            treatment = private$meta_data[, private$treatment_col]))
                        sce_mt <- calculateQCMetrics(sce_mt)
                    } else {
                        sce <- calculateQCMetrics(sce)
                    }
                    
                    # Frequency histogram of total UMI counts per column
                    
                    if (length(private$MT_rows)>0)
                    {
                        par(mfrow=c(1,2))
                        hist(sce$total_counts/1e3, xlab="Library sizes (thousands)", main="Histogram of total UMI counts",
                             breaks=20, col="grey80", ylab="Number of cells")
                        
                        hist(sce_mt$total_counts/1e3, xlab="Library sizes (thousands)", main="Histogram of total UMI counts (exclude MT)",
                             breaks=20, col="grey80", ylab="Number of cells")
                    } else {
                        hist(sce$total_counts/1e3, xlab="Library sizes (thousands)", main="Histogram of total UMI counts",
                             breaks=20, col="grey80", ylab="Number of cells")
                    }
                    
                    # Frequency histogram of features in each cell
                    
                    if (length(private$MT_rows)>0)
                    {
                        par(mfrow=c(1,2))
                        hist(sce$total_features_by_counts, xlab="Number of expressed genes", main="Number of expressed genes for each cell",
                             breaks=20, col="grey80", ylab="Number of cells")
                        hist(sce_mt$total_features_by_counts, xlab="Number of expressed genes",
                             main="Number of expressed genes for each cell", sub="excluded MT",
                             breaks=20, col="grey80", ylab="Number of cells")
                    } else {
                        hist(sce$total_features_by_counts, xlab="Number of expressed genes", main="Number of expressed genes for each cell",
                             breaks=20, col="grey80", ylab="Number of cells")
                    }
                    
                    # percentage of counts due to MT
                    if (length(private$MT_rows)>0)
                    {
                        par(mfrow=c(1,2))
                        hist(sce$pct_counts_Mt, xlab="Mitochondrial proportion (%)",
                             ylab="Number of cells", breaks=20, main="Mitochondrial proportion for each cell", col="grey80")
                        plot(sce$total_features_by_counts, sce$pct_counts_Mt,
                             xlab="Number of expressed genes",
                             ylab="percentage of total count due to MT",
                             main = "MT percentage vs number of expressed genes") # % of total UMI count due to MT UMI vs total features for each cell
                    }
                    
                    # Histogram of log average UMI for each gene after normalization
                    #unique_clusters = sort(unique(private$cluster))
                    #graph_col = ceiling(length(unique_clusters)/4)
                    #par(mfrow=c(4,graph_col))
                    #
                    #for (acluster in unique_clusters)
                    #{
                    #  one_cluster = private$norm_counts[-private$MT_rows,private$cluster_col == acluster]
                    #  log_av_UMI = log(apply(one_cluster,1,mean))
                    #  hist(log_av_UMI,breaks=20,main=paste0("log avUMI for each gene after normalization for cluster ", acluster),sub="excluded MT genes",cex.main=0.75)
                    #}
                    #
                    ## Scatterplot of log[(%CV)^2] vs log average UMI for each gene after normalization? for each cluster?
                    #par(mfrow=c(2,graph_col))
                    #
                    #for (acluster in unique_clusters)
                    #{
                    #  one_cluster = private$norm_counts[-private$MT_rows,private$cluster_col == acluster]
                    #  log_av_UMI = log(apply(one_cluster,1,mean))
                    #  logCVsq = log(apply(one_cluster,1,my_CV)^2)
                    #  plot(log_av_UMI,logCVsq,main = paste0("log(%CV^2) vs log avUMI for cluster ", acluster),sub="excluded MT genes",cex.main=0.75)
                    #}
                    
                    # find % of total UMI count for top 50 genes for each cell
                    if (length(private$MT_rows)>0)
                    {
                        par(mfrow=c(1,2))
                        hist(sce$pct_counts_in_top_50_features,breaks=20,sub="included MT genes",main="Percentage of top 50 features")
                        hist(sce_mt$pct_counts_in_top_50_features,breaks=20,sub="excluded MT genes",main="Percentage of top 50 features")
                    } else {
                        hist(sce$pct_counts_in_top_50_features,breaks=20,main="Percentage of top 50 features")
                    }
                    
                    # How many cells/percentage of cells have UMI count > 0 for each gene
                    par(mfrow=c(1,2))
                    cellnum = apply_sparseMatrix(private$counts, 1, function(x) sum(x>0, na.rm = TRUE))
                    hist(cellnum, breaks = 20, main = "# of cells with UMI > 0 for each gene",
                         xlab="# of cells with UMI count > 0")
                    
                    cellperc = apply_sparseMatrix(private$counts, 1, function(x) sum(x>0, na.rm = TRUE))/length(private$cells)*100
                    hist(cellperc, breaks = 20, main = "% of cells with UMI > 0 for each gene",
                         xlab="percentage of cells with UMI count > 0")
                    
                    if (length(private$MT_rows)>0)
                    {
                        par(mfrow=c(1,2))
                        cellnum = apply_sparseMatrix(private$counts[-private$MT_rows,], 1, function(x) sum(x>0, na.rm = TRUE))
                        hist(cellnum, breaks = 20, main = "# of cells with UMI > 0 for each gene",
                             sub="excluded MT genes",xlab="# of cells with UMI count > 0")
                        
                        cellperc = apply_sparseMatrix(private$counts[-private$MT_rows,], 1, function(x) sum(x>0, na.rm = TRUE))/length(private$cells)*100
                        hist(cellperc, breaks = 20, main = "% of cells with UMI > 0 for each gene",
                             sub="excluded MT genes",xlab="percentage of cells with UMI count > 0")
                    }
                    
                    # How many cells per subject
                    # Library sizes of cells for each subject
                    par(mfrow=c(1,2))
                    par(mar=c(8, 8, 1, 1))
                    barplot(table(sce$sample), las = 3, main = "Number of cells per sample")
                    boxplot(total_counts ~ sample, data = colData(sce), las = 3, log = "y", main = "Library size (log scale) by sample")
                    
                    dev.off()
                },
                # ------------------
                # apply_filter(): Apply various QC filters to the data
                #
                # Apply a MT gene filter? T or F
                # Filter cells by library size: # cells with library size < lib_size_low or library size > lib_size_high will be removed
                # min cells per gene filter: at least 50 cells per gene
                # min genes per cell filter: at least 500 genes per cell
                # perc_filter T or F: if True, use a percentage of total cells for the threshold on cells per gene filter
                # if false, use min.cells per gene threshold for cells per gene filter
                # min.perc.per.gene = percentage threshold
                # --------------------------
                apply_filter = function(MTfilter = TRUE, lib_size_low = 200, lib_size_high = 20*10^6,
                                        min.cells.per.gene = 50, min.genes.per.cell = 500,
                                        min.perc.cells.per.gene = 0.01, perc_filter = TRUE)
                {
                    stopifnot(is.logical(MTfilter))
                    stopifnot(is.logical(perc_filter))
                    
                    stopifnot(is.numeric(lib_size_low))
                    stopifnot(lib_size_low >= 0)
                    stopifnot(length(lib_size_low) == 1)
                    
                    stopifnot(is.numeric(lib_size_high))
                    stopifnot(length(lib_size_high) == 1)
                    
                    stopifnot(is.numeric(min.cells.per.gene))
                    stopifnot(length(min.cells.per.gene) == 1)
                    stopifnot(min.cells.per.gene >= 0)
                    
                    stopifnot(is.numeric(min.genes.per.cell))
                    stopifnot(length(min.genes.per.cell) == 1)
                    stopifnot(min.genes.per.cell >= 0)
                    
                    stopifnot(is.numeric(min.perc.cells.per.gene))
                    stopifnot(length(min.perc.cells.per.gene) == 1)
                    stopifnot(min.perc.cells.per.gene >= 0)
                    stopifnot(min.perc.cells.per.gene <= 1)
                    
                    if (MTfilter & length(private$MT_rows) > 0)
                    {
                        private$counts = private$counts[-private$MT_rows, , drop = FALSE]
                        private$genes = private$genes[-private$MT_rows]
                        private$MT_rows = NULL
                        private$filter_info[["MT_gene"]]<-private$genes[private$MT_rows]
                    }
                    
                    # remove bad genes
                    cell_num = apply_sparseMatrix(private$counts, 1, function(x) sum(x>0, na.rm = TRUE))
                    
                    if (perc_filter)
                    {
                        min.cells.per.gene = ceiling(min.perc.cells.per.gene*length(private$cells))
                    }
                    
                    bad_genes = which(cell_num < min.cells.per.gene)
                    private$filter_info[["less_cell_num"]] <- private$genes[bad_genes]
                    if (length(bad_genes) > 0)
                    {
                        private$counts = private$counts[-bad_genes, , drop = FALSE]
                        private$genes = private$genes[-bad_genes]
                        private$MT_rows = grep("^MT-", private$genes, ignore.case = TRUE)
                    }
                    
                    private$lib_sizes = Matrix::colSums(private$counts, na.rm = TRUE)
                    private$av_lib = mean(private$lib_sizes)
                    
                    # remove bad cells
                    bad_cells = which((private$lib_sizes < lib_size_low) | (private$lib_sizes > lib_size_high))
                    ngenes_by_cell = apply_sparseMatrix(private$counts, 2, function(x) sum(x>0, na.rm = TRUE))
                    bad_cells_2 = which(ngenes_by_cell < min.genes.per.cell)
                    bad_cells_union = union(bad_cells, bad_cells_2)
                    
                    if (length(bad_cells_union) > 0)
                    {
                        private$counts = private$counts[,-bad_cells_union, drop = FALSE]
                        private$cells = private$cells[-bad_cells_union]
                        private$meta_data = private$meta_data[-bad_cells_union, ,drop = FALSE]
                        private$lib_sizes = private$lib_sizes[-bad_cells_union]
                        private$av_lib = mean(private$lib_sizes)
                    }
                    
                    private$norm_counts = sweep_sparseMatrix(private$counts, 2, private$lib_sizes/private$av_lib,"/")
                    
                    # get pseudo-bulk counts
                    if (private$mode == "group") {
                        by_set <- split(seq_along(private$meta_data[, private$sampleId_col]), private$meta_data[, private$sampleId_col])
                        pseudoBulk_counts <- sapply(names(by_set), function(x) {
                            idx <- by_set[[x]]
                            pseudo.bulk.counts <- Matrix::rowSums(private$counts[, idx, drop = FALSE],  sum, na.rm = TRUE)
                        })
                    } else {
                        by_set <- split(seq_along(private$meta_data[, private$sampleId_col]), private$meta_data[, private$sampleId_col])
                        
                        pseudoBulk_counts_ref <- sapply(names(by_set), function(x) {
                            idx <- intersect(by_set[[x]], which(private$meta_data$.GrouP == "ref"))
                            pseudo.bulk.counts <- Matrix::rowSums(private$counts[, idx, drop = FALSE],  sum, na.rm = TRUE)
                        })
                        
                        pseudoBulk_counts_alt <- sapply(names(by_set), function(x) {
                            idx <- intersect(by_set[[x]], which(private$meta_data$.GrouP == "alt"))
                            pseudo.bulk.counts <- Matrix::rowSums(private$counts[, idx, drop = FALSE],  sum, na.rm = TRUE)
                        })
                        
                        pseudoBulk_counts <- merge(pseudoBulk_counts_ref, pseudoBulk_counts_alt, by = 0, suffixes = c(".ref", ".alt")) %>%
                            column_to_rownames(var = "Row.names")
                        
                        pseudoBulk_meta <- data.frame(PseudoBulk.Id = colnames(pseudoBulk_counts), stringsAsFactors = FALSE) %>%
                            dplyr::mutate(id_temp = gsub(".alt$", "", gsub('.ref$', "", PseudoBulk.Id)))
                        
                        names(pseudoBulk_meta)[2] <- private$sampleId_col
                        
                        private$pseudoBulk_meta <- pseudoBulk_meta
                        
                    }
                    
                    private$pseudoBulk_counts <- pseudoBulk_counts
                    cat("Dimensions of scRNAseq counts data after filtering...\n")
                    print(dim(private$counts))
                    cat("\n")
                    cat("Cluster summary...\n")
                    print(ftable(private$meta_data[, private$sampleId_col], private$meta_data[, private$cluster_col]))
                    cat("\n")
                    cat("Mitochondrial genes found: \n")
                    print(private$genes[private$MT_rows])
                    cat("\n")
                    cat("Dimensions of pseudo-bulk counts data after filtering: \n")
                    print(dim(private$pseudoBulk_counts))
                },
                # same as above, but returns an R6 object
                # perc_filter T or F: if True, use a percentage of total cells for the threshold on cells per gene filter
                # if false, use min.cells per gene threshold for cells per gene filter
                # min.perc.per.gene = percentage threshold
                apply_filter_R6 = function(MTfilter = TRUE, lib_size_low = 200, lib_size_high = 20*10^6,
                                           min.cells.per.gene = 50, min.genes.per.cell = 500,
                                           min.perc.cells.per.gene = 0.01, perc_filter = TRUE)
                {
                    stopifnot(is.logical(MTfilter))
                    stopifnot(is.logical(perc_filter))
                    
                    stopifnot(is.numeric(lib_size_low))
                    stopifnot(lib_size_low >= 0)
                    stopifnot(length(lib_size_low) == 1)
                    
                    stopifnot(is.numeric(lib_size_high))
                    stopifnot(length(lib_size_high) == 1)
                    
                    stopifnot(is.numeric(min.cells.per.gene))
                    stopifnot(length(min.cells.per.gene) == 1)
                    stopifnot(min.cells.per.gene >= 0)
                    
                    stopifnot(is.numeric(min.genes.per.cell))
                    stopifnot(length(min.genes.per.cell) == 1)
                    stopifnot(min.genes.per.cell >= 0)
                    
                    stopifnot(is.numeric(min.perc.cells.per.gene))
                    stopifnot(length(min.perc.cells.per.gene) == 1)
                    stopifnot(min.perc.cells.per.gene >= 0)
                    stopifnot(min.perc.cells.per.gene <= 1)
                    
                    filter_counts = private$counts
                    filter_genes = private$genes
                    temp_cells = private$cells
                    temp_colData = private$meta_data
                    
                    if (MTfilter & length(private$MT_rows) > 0)
                    {
                        filter_counts = filter_counts[-private$MT_rows, ,drop = FALSE]
                        filter_genes = filter_genes[-private$MT_rows]
                    }
                    
                    # remove bad genes
                    cell_num = apply_sparseMatrix(private$counts, 1, function(x) sum(x>0, na.rm = TRUE))
                    
                    if (perc_filter)
                    {
                        min.cells.per.gene = ceiling(min.perc.cells.per.gene*length(private$cells))
                    }
                    
                    bad_genes = which(cell_num < min.cells.per.gene)
                    private$filter_info[["less_cell_num"]] <- private$genes[bad_genes]
                    if (length(bad_genes) > 0)
                    {
                        filter_counts = filter_counts[-bad_genes, ,drop = FALSE]
                        filter_genes = filter_genes[-bad_genes]
                    }
                    
                    temp_lib_sizes = Matrix::colSums(filter_counts, na.rm = TRUE)
                    
                    # remove bad cells
                    bad_cells = which((temp_lib_sizes < lib_size_low) | (temp_lib_sizes > lib_size_high))
                    ngenes_by_cell = apply_sparseMatrix(private$counts, 2, function(x) sum(x>0, na.rm = TRUE))
                    bad_cells_2 = which(ngenes_by_cell < min.genes.per.cell)
                    bad_cells_union = union(bad_cells, bad_cells_2)
                    
                    if (length(bad_cells_union) > 0)
                    {
                        filter_counts = filter_counts[,-bad_cells_union, drop = FALSE]
                        temp_colData = temp_colData[-bad_cells_union, ,drop = FALSE]
                        temp_cells = temp_cells[-bad_cells_union]
                    }
                    
                    cat("Summary after filtering...\n")
                    return(BiostatsSingleCell$new(count_data = filter_counts, meta_data = temp_colData, cluster_col = private$cluster_col, 
                                                  sampleId_col = private$sampleId_col, treatment_col = private$treatment_col))
                },
                # Return an R6 object with subsets of counts and meta-data which meet the provided condition 
                # condition must be logical! For example: subset(RIN >= 7) or subset(RIN >=7, Region %in% 'motor-premotor')
                subset = function(condition)
                {
                    try(if(missing(condition)) stop("Please provide condition!"))
                    
                    #temp_counts = private$counts
                    temp_colData = private$meta_data
                    
                    e <- substitute(condition)
                    r <- eval(e, temp_colData, parent.frame())
                    if (!is.logical(r)) 
                        stop("'condition' must be logical")
                    idx <- r & !is.na(r)
                    
                    subset_colData <- temp_colData[idx, ]
                    subset_cells <- subset_colData[['cell']]
                    
                    idx.cells <- match(subset_cells, colnames(private$counts))
                    idx.cells <- idx.cells[!is.na(idx.cells)]
                    
                    subset_counts <- private$counts[, idx.cells]
                    
                    cat("Summary after subsetting...\n")
                    return(BiostatsSingleCell$new(count_data = subset_counts, meta_data = subset_colData, cluster_col = private$cluster_col, 
                                                  sampleId_col = private$sampleId_col, treatment_col = private$treatment_col))
                },
                down_sample = function(number)
                {
                    n_cells = length(private$cells)
                    
                    stopifnot(is.numeric(number))
                    stopifnot(number >= 1)
                    stopifnot(number <= n_cells)
                    
                    my_frac = number/n_cells
                    
                    subset_colData = data.frame(private$meta_data %>% 
                                                    group_by(!!sym(private$sampleId_col)) %>% 
                                                    sample_frac(my_frac) %>% ungroup())
                    
                    subset_cells <- subset_colData[['cell']]
                    
                    idx.cells <- match(subset_cells, colnames(private$counts))
                    idx.cells <- idx.cells[!is.na(idx.cells)]
                    
                    subset_counts <- private$counts[, idx.cells]
                    
                    cat("Summary after downsampling...\n")
                    return(BiostatsSingleCell$new(count_data = subset_counts, meta_data = subset_colData, cluster_col = private$cluster_col, 
                                                  sampleId_col = private$sampleId_col, treatment_col = private$treatment_col))
                },
                # set the data analysis mode
                # mode 1: group mode (DEGs between groups within each cell type)
                # mode 2: cluster mode (DEGs between clusters within each group)
                set_group_mode = function(cluster_of_interest, ref_group, alt_group)
                {
                    stopifnot(is.character(cluster_of_interest))
                    stopifnot(is.character(ref_group))
                    stopifnot(is.character(alt_group))
                    stopifnot(length(cluster_of_interest) == 1)
                    stopifnot(cluster_of_interest %in% unique(private$meta_data[, private$cluster_col]))
                    stopifnot(ref_group %in% unique(private$meta_data[, private$treatment_col]))
                    stopifnot(alt_group %in% unique(private$meta_data[, private$treatment_col]))
                    
                    cluster_index <- private$meta_data[, private$cluster_col] == cluster_of_interest
                    private$meta_data <- private$meta_data %>%
                        dplyr::mutate(.FilteR = cluster_index,
                                      .GrouP = NA)
                    
                    ref_index = which(private$meta_data[, private$treatment_col] %in% ref_group)
                    alt_index = which(private$meta_data[, private$treatment_col] %in% alt_group)
                    stopifnot(length(intersect(ref_index, alt_index))==0) 
                    try(if(length(intersect(ref_index, which(cluster_index)))==0) stop("Error, 0 cell in the reference group!"))
                    # stopifnot(length(intersect(ref_index, which(cluster_index)))>0) 
                    try(if(length(intersect(alt_index, which(cluster_index)))==0) stop("Error, 0 cell in the alternative group!"))
                    # stopifnot(length(intersect(alt_index, which(cluster_index)))>0) 
                    
                    private$meta_data$.GrouP[ref_index] <- "ref"
                    private$meta_data$.GrouP[alt_index] <- "alt"
                    private$mode = "group"
                    #   # if one of group sample size id 0 then set .GrouP=NA                             ADD by Grace
                    # private$meta_data<-private$meta_data%>%dplyr::group_by(patientID,predicted.id)%>%
                    #   dplyr::mutate(.GrouP1= ifelse(sum(.GrouP=="ref")==0 | sum(.GrouP=="alt")==0 ,NA, 1))%>%
                    #   ungroup()%>%dplyr:::mutate(.GrouP=ifelse(.GrouP1==1,.GrouP,"EMPTY"))%>%dplyr::select(-.GrouP1)
                },
                set_cluster_mode = function(group_of_interest, ref_cluster, alt_cluster = "others")
                {
                    stopifnot(is.character(group_of_interest))
                    stopifnot(is.character(ref_cluster))
                    stopifnot(is.character(alt_cluster))
                    stopifnot(length(group_of_interest) == 1)
                    stopifnot(group_of_interest %in% unique(private$meta_data[, private$cluster_col]))
                    stopifnot(ref_cluster %in% unique(private$meta_data[, private$treatment_col]))
                    stopifnot(alt_cluster %in% c(unique(private$meta_data[, private$treatment_col]), "others"))
                    
                    group_index <- private$meta_data[, private$cluster_col] == group_of_interest
                    private$meta_data <- as.data.frame(private$meta_data) %>%
                        dplyr::mutate(.FilteR = group_index,
                                      .GrouP = NA)
                    
                    ref_index = which(private$meta_data[, private$treatment_col] %in% ref_cluster)
                    if (alt_cluster == "others") {
                        alt_index = which(!(private$meta_data[, private$treatment_col] %in% ref_cluster))
                    } else {
                        alt_index = which(private$meta_data[, private$treatment_col] %in% alt_cluster)
                    }
                    stopifnot(length(intersect(ref_index, alt_index))==0) 
                    stopifnot(length(intersect(ref_index, which(group_index)))>0) 
                    stopifnot(length(intersect(alt_index, which(group_index)))>0) 
                    
                    private$meta_data$.GrouP[ref_index] <- "ref"
                    private$meta_data$.GrouP[alt_index] <- "alt"
                    private$mode = "cluster"
                    #   # if one of group sample size id 0 then set .GrouP=NA                             ADD by Grace
                    private$meta_data<-private$meta_data%>%dplyr::group_by(patientID,predicted.id)%>%
                        dplyr::mutate(.GrouP1= ifelse(sum(.GrouP=="ref")==0 | sum(.GrouP=="alt")==0 ,NA, 1))%>%
                        ungroup()%>%dplyr:::mutate(.GrouP=ifelse(.GrouP1==1,.GrouP,"EMPTY"))%>%dplyr::select(-.GrouP1)
                    
                },
                # -------------------------
                # apply a filter to two groups (e.g. treated vs. placebo) for one cluster of interest
                # keep only genes with at least 50 cells (or at least 10% of cells (if larger than 50) from the smaller group) for both groups
                # and: both groups must have at least 50 cells  
                # or: either group (or both) must have at least 50 cells (not recommended)
                # filter out genes if both groups have 75% percentile at 0
                # filter out subjects if number of cells is less than 5
                # --------------
                apply_filter_contrasts_R6 = function(min.cells.per.gene = 50, min.perc.cells.per.gene = 0.10, perc.cells.filter = T,
                                                     min.cells.per.gene.type = "and", cells.per.gene.filter = T, 
                                                     perc.filter = T, perc.filter.type = "and", perc_threshold = 0.75, 
                                                     min.ave.pseudo.bulk.cpm = 1, pseudo.bulk.cpm.filter = T, min.cells.per.subj = 5)
                {
                    stopifnot(is.logical(cells.per.gene.filter))
                    stopifnot(is.logical(perc.filter))
                    stopifnot(is.logical(pseudo.bulk.cpm.filter))
                    stopifnot(is.logical(perc.cells.filter))
                    stopifnot(is.numeric(min.cells.per.gene))
                    stopifnot(length(min.cells.per.gene) == 1)
                    stopifnot(min.cells.per.gene >= 0)
                    stopifnot(is.numeric(min.ave.pseudo.bulk.cpm))
                    stopifnot(length(min.ave.pseudo.bulk.cpm) == 1)
                    stopifnot(min.ave.pseudo.bulk.cpm >= 0)
                    stopifnot(is.numeric(min.perc.cells.per.gene))
                    stopifnot(length(min.perc.cells.per.gene) == 1)
                    stopifnot(min.perc.cells.per.gene >= 0)
                    stopifnot(min.perc.cells.per.gene <= 1)
                    stopifnot(is.numeric(perc_threshold))
                    stopifnot(length(perc_threshold) == 1)
                    stopifnot(perc_threshold >= 0 & perc_threshold <= 1)
                    stopifnot(min.cells.per.gene.type == "or" | min.cells.per.gene.type == "and")
                    stopifnot(perc.filter.type == "or" | perc.filter.type == "and")
                    stopifnot(is.numeric(min.cells.per.subj))
                    stopifnot(min.cells.per.subj >= 0)
                    
                    filter_counts = private$counts
                    filter_genes = private$genes
                    filter_cells = private$cells
                    filter_meta_data = as.data.frame(private$meta_data)
                    mode = private$mode
                    
                    # only keep cells from the selected cluster
                    cluster_index <- which(filter_meta_data[, '.FilteR'] == TRUE)
                    
                    placebo_index = which(filter_meta_data[, '.GrouP'] %in% "ref")
                    treated_index = which(filter_meta_data[, '.GrouP'] %in% "alt")
                    stopifnot(length(intersect(placebo_index,treated_index))==0)
                    
                    all_index = c(intersect(cluster_index, placebo_index),
                                  intersect(cluster_index, treated_index)) # Just in case there are more than 2 groups
                    filter_meta_data <- filter_meta_data[all_index, ]
                    bad_genes = NULL
                    bad_genes_2 = NULL
                    
                    if (cells.per.gene.filter)
                    {
                        if (perc.cells.filter)
                        {
                            min.cells.per.gene = max(min.cells.per.gene, 
                                                     ceiling(min.perc.cells.per.gene*min(length(intersect(placebo_index,cluster_index)),length(intersect(treated_index,cluster_index)))))
                            message(paste0("Minimum cells per gene was set to ", min.cells.per.gene))
                        }
                        
                        placebo_count = apply_sparseMatrix(filter_counts[,intersect(cluster_index, placebo_index)], 1, function(x) sum(x > 0, na.rm = T))
                        treated_count = apply_sparseMatrix(filter_counts[,intersect(cluster_index, treated_index)], 1, function(x) sum(x > 0, na.rm = T))
                        if (min.cells.per.gene.type == "and") bad_genes = which(placebo_count < min.cells.per.gene | treated_count < min.cells.per.gene)
                        if (min.cells.per.gene.type == "or") bad_genes = which(placebo_count < min.cells.per.gene & treated_count < min.cells.per.gene)
                        private$filter_info[["less_per_cell_num_group"]]<-filter_genes[bad_genes]    # add by Grace 
                    }
                    
                    if (perc.filter)
                    {
                        placebo_perc = apply_sparseMatrix(private$norm_counts[,intersect(cluster_index, placebo_index)], 1, function(x) quantile(x, probs=perc_threshold, na.rm = T))
                        treated_perc = apply_sparseMatrix(private$norm_counts[,intersect(cluster_index, treated_index)], 1, function(x) quantile(x, probs=perc_threshold, na.rm = T))
                        if (perc.filter.type == "and") bad_genes_2 = which(placebo_perc == 0 & treated_perc == 0) # Both groups have 75% percentile at 0
                        if (perc.filter.type == "or") bad_genes_2 = which(placebo_perc == 0 | treated_perc == 0) # At least one group has 75% percentile at 0 
                        private$filter_info[["75%_percentile_zero_group"]]<-filter_genes[bad_genes_2]    # add by Grace
                    }
                    
                    bad_list = union(bad_genes, bad_genes_2)
                    if (length(bad_list) > 0)
                    {
                        filter_counts = filter_counts[-bad_list, ,drop = FALSE]
                        filter_genes = filter_genes[-bad_list]
                    }
                    
                    filter_counts = filter_counts[, all_index]
                    filter_cells = filter_cells[all_index]
                    
                    # get number of cells for each subject
                    ncells_subj <- table(filter_meta_data[, private$sampleId_col])
                    
                    # find subjects with number of cells < threshold (5 by default)
                    low_cells <- ncells_subj[which(ncells_subj < min.cells.per.subj)]
                    
                    # drop those subjects
                    if(length(low_cells) > 0) 
                    {
                        bad_subjs <- which(filter_meta_data[, private$sampleId_col] %in% names(low_cells))
                        filter_counts = filter_counts[, -bad_subjs, drop = FALSE]
                        filter_meta_data = filter_meta_data[-bad_subjs, ]
                        if (length(low_cells) == 1)
                        {
                            message(glue("Subject: {names(low_cells)} was excluded due to limited number of cells!"))
                        } else {
                            message(glue("{length(low_cells)} subjects ({paste(names(low_cells), collapse = ', ')}) were excluded due to limited number of cells!"))
                        }
                    }
                    
                    # get sample size (number of subjects) for each group
                    nsubjs_by_group <- table(unique(filter_meta_data[, c(private$sampleId_col, ".GrouP")])$.GrouP)
                    if (!(length(nsubjs_by_group) == 2 & min(nsubjs_by_group) >= 2)) stop("Error! Each group needs to have at least 2 subjects!")
                    
                    if (pseudo.bulk.cpm.filter)
                    {
                        # get pseudo-bulk counts
                        if (private$mode == "group") {
                            by_set <- split(seq_along(filter_meta_data[, private$sampleId_col]), filter_meta_data[, private$sampleId_col])
                            pseudoBulk_counts <- sapply(names(by_set), function(x) {
                                idx <- by_set[[x]]
                                pseudo.bulk.counts <- Matrix::rowSums(filter_counts[, idx, drop = FALSE],  sum, na.rm = TRUE)
                            })
                        } else {
                            by_set <- split(seq_along(filter_meta_data[, private$sampleId_col]), filter_meta_data[, private$sampleId_col])
                            # reference cluster
                            pseudoBulk_counts_ref <- sapply(names(by_set), function(x) {
                                idx <- intersect(by_set[[x]], which(filter_meta_data$.GrouP == "ref"))
                                pseudo.bulk.counts <- Matrix::rowSums(filter_counts[, idx, drop = FALSE],  sum, na.rm = TRUE)
                            })
                            # alternative cluster
                            pseudoBulk_counts_alt <- sapply(names(by_set), function(x) {
                                idx <- intersect(by_set[[x]], which(filter_meta_data$.GrouP == "alt"))
                                pseudo.bulk.counts <- Matrix::rowSums(filter_counts[, idx, drop = FALSE],  sum, na.rm = TRUE)
                            })
                            
                            pseudoBulk_counts <- merge(pseudoBulk_counts_ref, pseudoBulk_counts_alt, by = 0, suffixes = c(".ref", ".alt")) %>%
                                column_to_rownames(var = "Row.names")
                        }
                        
                        # get cpm and do the filtering
                        cpm.counts <- cpm(pseudoBulk_counts)
                        mean.cpm <- apply(cpm.counts, 1, mean)
                        low_cpm_genes = which(mean.cpm < min.ave.pseudo.bulk.cpm)
                        private$filter_info[["low_cpm"]]<-filter_genes[low_cpm_genes]    # add by Grace
                        cat("Summary after filtering...\n")
                        if (length(low_cpm_genes > 0))
                        {
                            print(paste0(length(low_cpm_genes), " were filtered out due to low average CPM"))
                            return(BiostatsSingleCell$new(count_data = filter_counts[-low_cpm_genes,], meta_data = filter_meta_data, cluster_col = private$cluster_col, 
                                                          sampleId_col = private$sampleId_col, treatment_col = private$treatment_col, mode = mode))
                        } else {
                            return(BiostatsSingleCell$new(count_data = filter_counts, meta_data = filter_meta_data, cluster_col = private$cluster_col, 
                                                          sampleId_col = private$sampleId_col, treatment_col = private$treatment_col, mode = mode))
                        }
                    } else {
                        print("Summary after filtering...")
                        return(BiostatsSingleCell$new(count_data = filter_counts, meta_data = filter_meta_data, cluster_col = private$cluster_col, 
                                                      sampleId_col = private$sampleId_col, treatment_col = private$treatment_col, mode = mode))
                    }
                },
                # run t-test pipeline
                # return sorted results
                t_test_pipeline = function()
                {
                    pseudoBulk <- private$pseudoBulk_counts
                    lib_sizes = apply(pseudoBulk, 2, sum)
                    av_lib = mean(lib_sizes)
                    norm_counts = sweep(pseudoBulk, 2, lib_sizes/av_lib, "/")
                    sample_ids <- colnames(norm_counts)
                    
                    if (private$mode == "group") {
                        meta_data <- unique(private$meta_data[, c(private$sampleId_col, ".GrouP")])
                        treated_index <- which(sample_ids %in% meta_data[[private$sampleId_col]][meta_data[[".GrouP"]] %in% 'alt'])
                        placebo_index <- which(sample_ids %in% meta_data[[private$sampleId_col]][meta_data[[".GrouP"]] %in% 'ref'])
                        t_pvalues = apply(norm_counts, 1, function(x) t.test(x[treated_index], x[placebo_index])$p.value)
                    } else {
                        meta_data <- private$pseudoBulk_meta
                        treated_index <- which(sample_ids %in% meta_data[["PseudoBulk.Id"]][meta_data[[".GrouP"]] %in% 'alt'])
                        placebo_index <- which(sample_ids %in% meta_data[["PseudoBulk.Id"]][meta_data[[".GrouP"]] %in% 'ref'])
                        t_pvalues = apply(norm_counts, 1, function(x) t.test(x[treated_index], x[placebo_index], paired = TRUE)$p.value) # paired t-test
                    }
                    
                    t_fdr = p.adjust(t_pvalues,method="fdr")
                    t_FC = my_FC(norm_counts, treated_index, placebo_index)
                    t_log2FC = log2(t_FC)
                    t_log2FC[is.infinite(t_log2FC)] = NA
                    t_log2FC[is.nan(t_log2FC)] = NA
                    
                    t_table = data.frame("ID" = private$genes, "log2FC" = t_log2FC, "Pvalue" = t_pvalues, "FDR" = t_fdr)
                    private$de_results[["t.test"]] = t_table
                    private$de_method = "t.test"
                    return(t_table %>% dplyr:::arrange(FDR))
                },
                # run ANCOVA pipeline
                # return sorted results
                ancova_pipeline = function(covs = NULL)
                {
                    stopifnot(is.null(covs) | is.character(covs))
                    
                    pseudoBulk <- private$pseudoBulk_counts
                    lib_sizes = apply(pseudoBulk, 2, sum)
                    av_lib = mean(lib_sizes)
                    norm_counts = sweep(pseudoBulk, 2, lib_sizes/av_lib, "/")
                    genes <- private$genes
                    meta_data <- private$meta_data
                    
                    if (!is.null(covs)) {
                        try(if(all(covs %in% names(meta_data)) == FALSE) stop("Error: can not find all covariates from the meta data!"))
                    }
                    
                    if (private$mode == "group") {
                        sample_ids <- data.frame(sampleID = colnames(norm_counts), stringsAsFactors = FALSE)
                        names(sample_ids) = private$sampleId_col
                        meta_data <- unique(private$meta_data[, c(private$sampleId_col, '.GrouP', covs)])
                    } else {
                        sample_ids <- private$pseudoBulk_meta
                        meta_data <- unique(private$meta_data[, c(private$sampleId_col, covs)])
                    }
                    meta_data <- sample_ids %>%
                        dplyr:::left_join(meta_data)
                    
                    other_group <- setdiff(unique(meta_data[, '.GrouP']), c('ref', 'alt'))
                    meta_data[, '.GrouP'] <- factor(meta_data[, '.GrouP'], 
                                                    levels = c('ref', 'alt', other_group))
                    meta_data[, '.GrouP'] <- relevel(meta_data[, '.GrouP'], ref = "ref")
                    
                    placebo_index <- which(meta_data[['.GrouP']] %in% 'ref')
                    treated_index <- which(meta_data[['.GrouP']] %in% 'alt')
                    col_index = c(placebo_index, treated_index)
                    group_label <- paste0('.GrouP', "alt")
                    ancova_table<-data.frame()
                    
                    for (gene in genes) {
                        
                        counts_frame = data.frame(counts = as.vector(t( norm_counts[gene, col_index])))
                        
                        glmer_frame =dplyr::: bind_cols(meta_data, counts_frame)
                        
                        #gene_id <- private$genes[gene]
                        
                        fml <- paste0("counts ~ ",paste(c( covs,'.GrouP'), collapse = " + "))
                        
                        ancova = lm(as.formula(fml),data  = glmer_frame)
                        ls.mod1 <- emmeans(ancova , pairwise ~ '.GrouP')
                        resw <- as.data.frame(summary(ls.mod1, type = "response", infer = TRUE))
                        a_p_val <- summary(ancova)$coefficients[group_label, 4]
                        a_FC <- resw$emmeans.emmean[[2]]/resw$emmeans.emmean[[1]]
                        a_log2FC = log2(a_FC)
                        a_log2FC[is.infinite(a_log2FC)] = NA
                        a_log2FC[is.nan(a_log2FC)] = NA
                        a_table = data.frame("ID" = gene , "log2FC" = a_log2FC, "Pvalue" = a_p_val)
                        ancova_table<-rbind(ancova_table,a_table)
                    }
                    
                    ancova_table$FDR = p.adjust(ancova_table$Pvalue, method="fdr")
                    
                    private$de_results[["ANCOVA"]] <- ancova_table
                    private$de_method = "ANCOVA"
                    
                    return(ancova_table %>% dplyr:::arrange(FDR))
                },
                # Run U-test pipeline
                # return sorted results
                u_test_pipeline = function()
                {
                    pseudoBulk <- private$pseudoBulk_counts
                    lib_sizes = apply(pseudoBulk, 2, sum)
                    av_lib = mean(lib_sizes)
                    norm_counts = sweep(pseudoBulk, 2, lib_sizes/av_lib, "/")
                    sample_ids <- colnames(norm_counts)
                    
                    if (private$mode == "group") {
                        meta_data <- unique(private$meta_data[, c(private$sampleId_col, ".GrouP")])
                        treated_index <- which(sample_ids %in% meta_data[[private$sampleId_col]][meta_data[[".GrouP"]] %in% 'alt'])
                        placebo_index <- which(sample_ids %in% meta_data[[private$sampleId_col]][meta_data[[".GrouP"]] %in% 'ref'])
                        u_pvalues = apply(norm_counts, 1, function(x) wilcox.test(x[treated_index],x[placebo_index])$p.value)
                    } else {
                        meta_data <- private$pseudoBulk_meta
                        treated_index <- which(sample_ids %in% meta_data[["PseudoBulk.Id"]][meta_data[[".GrouP"]] %in% 'alt'])
                        placebo_index <- which(sample_ids %in% meta_data[["PseudoBulk.Id"]][meta_data[[".GrouP"]] %in% 'ref'])
                        u_pvalues = apply(norm_counts, 1, function(x) wilcox.test(x[treated_index],x[placebo_index], paired = TRUE)$p.value) # paired u-test
                    }
                    
                    u_fdr = p.adjust(u_pvalues,method="fdr")
                    u_FC = my_FC(norm_counts, treated_index, placebo_index)
                    u_log2FC = log2(u_FC)
                    u_log2FC[is.infinite(u_log2FC)] = NA
                    u_log2FC[is.nan(u_log2FC)] = NA
                    
                    u_table = data.frame("ID" = private$genes, "log2FC" = u_log2FC, "Pvalue" = u_pvalues, "FDR" = u_fdr)
                    private$de_results[["u.test"]] = u_table
                    private$de_method = "u.test"
                    return(u_table %>% dplyr:::arrange(FDR))
                },
                # Run pooled U-test pipeline (cell level counts data ignoring sample ids)
                # return sorted results
                # NOT RECOMMENDED TO USE (might want to remove this function from the code)
                pooled_u_test_pipeline = function(placebo_group, treated_group)
                {
                    cell_ids <- colnames(private$norm_counts)
                    meta_data <- private$meta_data[, c("cell", '.GrouP')]
                    treated_index <- which(cell_ids %in% meta_data[["cell"]][meta_data[[".GrouP"]] %in% 'alt'])
                    placebo_index <- which(cell_ids %in% meta_data[["cell"]][meta_data[[".GrouP"]] %in% 'ref'])
                    u_pvalues = apply_sparseMatrix(private$norm_counts, 1, function(x) wilcox.test(x[treated_index],x[placebo_index])$p.value)
                    u_fdr = p.adjust(u_pvalues,method="fdr")
                    u_FC = my_FC(private$norm_counts, treated_index, placebo_index)
                    u_log2FC = log2(u_FC)
                    u_log2FC[is.infinite(u_log2FC)] = NA
                    u_log2FC[is.nan(u_log2FC)] = NA
                    
                    u_table = data.frame("ID" = row.names(private$norm_counts), "log2FC" = u_log2FC, "Pvalue" = u_pvalues, "FDR" = u_fdr)
                    private$de_results[["u.pool.test"]] = u_table
                    private$de_method = "u.pool.test"
                    return(u_table %>% dplyr:::arrange(FDR))
                },
                # run edge-R pipeline
                # return sorted results
                edgeR_pipeline = function(covs = NULL)
                {
                    stopifnot(is.null(covs) | is.character(covs))
                    
                    pseudoBulk <- private$pseudoBulk_counts
                    
                    if (!is.null(covs)) {
                        try(if(all(covs %in% names(private$meta_data)) == FALSE) stop("Error: can not find all covariates from the meta data!"))
                    }
                    
                    if (private$mode == "group") {
                        sample_ids <- data.frame(sampleID = colnames(pseudoBulk), stringsAsFactors = FALSE)
                        names(sample_ids) = private$sampleId_col
                        meta_data <- unique(private$meta_data[, c(private$sampleId_col, '.GrouP', covs)])
                    } else {
                        sample_ids <- private$pseudoBulk_meta
                        meta_data <- unique(private$meta_data[, c(private$sampleId_col, covs)])
                    }
                    
                    meta_data <- sample_ids %>%
                        dplyr:::left_join(meta_data)
                    
                    other_group <- setdiff(unique(meta_data[, '.GrouP']), c('ref', 'alt'))
                    meta_data[, '.GrouP'] <- factor(meta_data[, '.GrouP'], 
                                                    levels = c('ref', 'alt', other_group))
                    meta_data[, '.GrouP'] <- relevel(meta_data[, '.GrouP'], ref = "ref")
                    
                    cat("Forming DGEList...\n")
                    edgeR_object = DGEList(counts = pseudoBulk, group = meta_data[, '.GrouP'])
                    
                    cat("Calculating NormFactors...\n")
                    edgeR_object = calcNormFactors(edgeR_object)
                    
                    cat("Creating design matrix...\n")
                    if (private$mode == "group") {                                                 # add by Grace
                        fml <- paste0("~ ", paste(c('.GrouP', covs), collapse = " + "))
                    } else{
                        fml <- paste0("~ ", paste(c('.GrouP', private$sampleId_col), collapse = " + ")) 
                    }
                    design <- model.matrix(as.formula(fml), data = meta_data)
                    
                    cat("Estimating dispersions...\n")
                    edgeR_object <- estimateDisp(edgeR_object, design)
                    
                    cat("Fitting the model...\n")
                    fit <- glmQLFit(edgeR_object, design)
                    
                    cat("Running QLFTest...\n")
                    
                    idx <- grep(paste0('.GrouP', "alt"), colnames(design))
                    edgeR_results = glmQLFTest(fit, coef = idx)
                    edgeR_all = data.frame(topTags(edgeR_results, n = nrow(pseudoBulk))$table)
                    edgeR_all$ID = row.names(edgeR_all)
                    row.names(edgeR_all) = NULL
                    edgeR_all = edgeR_all[,c(6, 1:5)] %>%
                        dplyr::rename("log2FC" = "logFC",
                                      "Pvalue" = "PValue")
                    private$de_results[["edgeR"]] = edgeR_all
                    private$de_method = "edgeR"
                    
                    return(edgeR_all %>% dplyr:::arrange(FDR))
                },
                # limma-voom: TMM
                # return sorted results
                limma_pipeline = function(covs = NULL)
                {
                    stopifnot(is.null(covs) | is.character(covs))
                    
                    if (!is.null(covs)) {
                        try(if(all(covs %in% names(private$meta_data)) == FALSE) stop("Error: can not find all covariates from the meta data!"))
                    }
                    
                    pseudoBulk <- private$pseudoBulk_counts
                    
                    if (private$mode == "group") {
                        sample_ids <- data.frame(sampleID = colnames(pseudoBulk), stringsAsFactors = FALSE)
                        names(sample_ids) = private$sampleId_col
                        meta_data <- unique(private$meta_data[, c(private$sampleId_col, '.GrouP', covs)])
                    } else {
                        sample_ids <- private$pseudoBulk_meta
                        meta_data <- unique(private$meta_data[, c(private$sampleId_col, covs)])
                    }
                    
                    meta_data <- sample_ids %>%
                        dplyr:::left_join(meta_data)
                    
                    other_group <- setdiff(unique(meta_data[, '.GrouP']), c('ref', 'alt'))
                    meta_data[, '.GrouP'] <- factor(meta_data[, '.GrouP'], 
                                                    levels = c('ref', 'alt', other_group))
                    meta_data[, '.GrouP'] <- relevel(meta_data[, '.GrouP'], ref = "ref")
                    
                    cat("Forming DGEList...\n")
                    edgeR_object = DGEList(counts = pseudoBulk, group = meta_data[, '.GrouP'])
                    
                    cat("Calculating NormFactors...\n")
                    edgeR_object = calcNormFactors(edgeR_object)
                    
                    cat("Creating design matrix...\n")
                    
                    if (private$mode == "group") {                                                    # add by Grace
                        fml <- paste0("~ ", paste(c('.GrouP', covs), collapse = " + "))
                    } else{
                        fml <- paste0("~ ", paste(c('.GrouP', private$sampleId_col), collapse = " + ")) 
                    }
                    
                    design <- model.matrix(as.formula(fml), data = meta_data)
                    
                    cat("Running voom...\n")
                    v <- voom(edgeR_object, design)
                    #v <- voomWithQualityWeights(edgeR_object, design)
                    
                    cat("Linear Model Fit...\n")
                    fit <- lmFit(v, design)
                    
                    idx <- grep(paste0('.GrouP', "alt"), colnames(design))
                    
                    cat("Emperical Bayes...\n")
                    fit2 <- eBayes(fit)
                    
                    cat("Top Table...\n")
                    limma_table = topTable(fit2, coef = idx, number = nrow(pseudoBulk))
                    limma_table$ID = row.names(limma_table)
                    row.names(limma_table) = NULL
                    limma_table = limma_table[,c(7,1:6)] %>%
                        dplyr::rename("log2FC" = "logFC",
                                      "Pvalue" = "P.Value",
                                      "FDR" = "adj.P.Val")
                    private$de_results[["limma"]] <- limma_table
                    private$de_method = "limma"
                    
                    return(limma_table %>% dplyr:::arrange(FDR))
                },
                # limma-voom cell level: TMM
                # return sorted results
                limma_cell_level_pipeline = function(covs = NULL, covs_formula = NULL,remove_subj=TRUE,cell_level_covs=NULL,cell_level_covs_formula=NULL) 
                {
                    stopifnot(is.null(covs) | is.character(covs))
                    stopifnot(is.null(covs_formula) | is.character(covs_formula))
                    stopifnot(is.null(cell_level_covs) | is.character(cell_level_covs))
                    stopifnot(is.null(cell_level_covs_formula) | is.character(cell_level_covs_formula))
                    
                    genes <- private$genes
                    counts <- as.matrix(private$counts)
                    meta_data <- private$meta_data
                    
                    if (!is.null(covs)) {
                        try(if(all(covs %in% names(meta_data)) == FALSE) stop("Error: can not find all covariates from the meta data!"))
                    }
                    
                    if (!is.null(cell_level_covs)) {
                        try(if(all(cell_level_covs %in% names(meta_data)) == FALSE) stop("Error: can not find all covariates from the meta data!"))
                    }
                    
                    other_group <- setdiff(unique(meta_data[, '.GrouP']), c('ref', 'alt'))
                    meta_data[, '.GrouP'] <- factor(meta_data[, '.GrouP'], 
                                                    levels = c('ref', 'alt', other_group))
                    meta_data[, '.GrouP'] <- relevel(meta_data[, '.GrouP'], ref = "ref")
                    
                    cat("Forming DGEList...\n")
                    edgeR_object = DGEList(counts = counts, group = meta_data[, '.GrouP'])
                    
                    cat("Calculating NormFactors...\n")
                    edgeR_object = calcNormFactors(edgeR_object)
                    
                    cat("Creating design matrix...\n")
                    
                    if (is.null(covs_formula)) 
                    {
                        if(is.null(cell_level_covs_formula))
                        {
                            fml <- paste0("~ ", paste(c('.GrouP', covs,cell_level_covs), collapse = " + "))
                        } else {
                            fml <- paste0("~ ", paste(c('.GrouP', covs,cell_level_covs_formula), collapse = " + "))
                            
                        }
                    } else {
                        if(is.null(cell_level_covs_formula))
                        {
                            fml <- paste0("~ ", paste(c('.GrouP', covs_formula,cell_level_covs), collapse = " + "))
                        } else {
                            fml <- paste0("~ ", paste(c('.GrouP', covs_formula,cell_level_covs_formula), collapse = " + "))
                            
                        }
                    }
                    design <- model.matrix(as.formula(fml), data = meta_data)
                    
                    if (remove_subj)
                    {
                        if (is.null(cell_level_covs_formula)){
                            fml <- paste0("~ 0 + ", paste(c(private$sampleId_col,cell_level_covs),collapse= " + "))
                        } else {
                            fml <- paste0("~ 0 + ", paste(c(private$sampleId_col,cell_level_covs_formula),collapse= " + "))
                        }
                    } else {
                        # when remove subj = FALSE, fml will be same as calculated for design variable above so donât change value of fml for this case
                    }
                    
                    design_voom <- model.matrix(as.formula(fml), data = meta_data)
                    print(head(design_voom))
                    
                    cat("Running voom...\n")
                    v <- voom(edgeR_object, design_voom)
                    
                    individual <- meta_data[, private$sampleId_col]
                    dupcor <- duplicateCorrelation(v, design, block = individual)
                    
                    cat("Linear Model Fit...\n")
                    fit <- lmFit(v, design, block = individual, correlation = dupcor$consensus)
                    
                    idx <- grep(paste0('.GrouP', "alt"), colnames(design))
                    
                    cat("Emperical Bayes...\n")
                    fit2 <- eBayes(fit)
                    
                    cat("Top Table...\n")
                    limma_table = topTable(fit2, coef = idx, number = nrow(counts))
                    limma_table$ID = row.names(limma_table)
                    row.names(limma_table) = NULL
                    limma_table = limma_table[,c(7,1:6)] %>%
                        dplyr::rename("log2FC" = "logFC",
                                      "Pvalue" = "P.Value",
                                      "FDR" = "adj.P.Val")
                    private$de_results[["limma_cell"]] <- limma_table
                    private$de_method = "limma_cell"
                    
                    return(limma_table %>% dplyr:::arrange(FDR))
                },
                
                # limma_cell_level_pipeline = function(covs = NULL)
                # {
                #   stopifnot(is.null(covs) | is.character(covs))
                # 
                #   genes <- private$genes
                #   counts <- as.matrix(private$counts)
                #   meta_data <- private$meta_data
                # 
                #   if (!is.null(covs)) {
                #       try(if(all(covs %in% names(meta_data)) == FALSE) stop("Error: can not find all covariates from the meta data!"))
                #   }
                #   
                #   #if (private$mode == "group") {
                #   #  sample_ids <- data.frame(sampleID = colnames(pseudoBulk), stringsAsFactors = FALSE)
                #   #  names(sample_ids) = private$sampleId_col
                #   #  meta_data <- unique(private$meta_data[, c(private$sampleId_col, '.GrouP', covs)])
                #   #} else {
                #   #  sample_ids <- private$pseudoBulk_meta
                #   #  meta_data <- unique(private$meta_data[, c(private$sampleId_col, covs)])
                #   #}
                #   
                #   #meta_data <- sample_ids %>%
                #   #  dplyr:::left_join(meta_data)
                #   
                #   other_group <- setdiff(unique(meta_data[, '.GrouP']), c('ref', 'alt'))
                #   meta_data[, '.GrouP'] <- factor(meta_data[, '.GrouP'], 
                #                                   levels = c('ref', 'alt', other_group))
                #   meta_data[, '.GrouP'] <- relevel(meta_data[, '.GrouP'], ref = "ref")
                #   
                #   cat("Forming DGEList...\n")
                #   edgeR_object = DGEList(counts = counts, group = meta_data[, '.GrouP'])
                #   
                #   cat("Calculating NormFactors...\n")
                #   edgeR_object = calcNormFactors(edgeR_object)
                #   
                #   cat("Creating design matrix...\n")
                #   
                #   fml <- paste0("~ ", paste(c('.GrouP', covs), collapse = " + "))
                #   design <- model.matrix(as.formula(fml), data = meta_data)
                #   
                #   
                #   cat("Running voom...\n")
                #   v <- voom(edgeR_object, design)
                #   #v <- voomWithQualityWeights(edgeR_object, design)
                # 
                #   individual <- meta_data[, private$sampleId_col]
                #   dupcor <- duplicateCorrelation(v, design, block = individual)
                #   
                #   cat("Linear Model Fit...\n")
                #   fit <- lmFit(v, design, block = individual, correlation = dupcor$consensus)
                #   
                #   idx <- grep(paste0('.GrouP', "alt"), colnames(design))
                #   
                #   cat("Emperical Bayes...\n")
                #   fit2 <- eBayes(fit)
                #   
                #   cat("Top Table...\n")
                #   limma_table = topTable(fit2, coef = idx, number = nrow(counts))
                #   limma_table$ID = row.names(limma_table)
                #   row.names(limma_table) = NULL
                #   limma_table = limma_table[,c(7,1:6)] %>%
                #     dplyr::rename("log2FC" = "logFC",
                #                   "Pvalue" = "P.Value",
                #                   "FDR" = "adj.P.Val")
                #   private$de_results[["limma_cell"]] <- limma_table
                #   private$de_method = "limma_cell"
                #   
                #   return(limma_table %>% dplyr:::arrange(FDR))
                # },
                # run DESeq2 pipeline
                # return sorted results
                DESeq2_pipeline = function(covs = NULL, shrink = TRUE, parallel = FALSE)
                {
                    stopifnot(is.logical(shrink))
                    stopifnot(is.logical(parallel))
                    stopifnot(is.null(covs) | is.character(covs))
                    
                    if (!is.null(covs)) {
                        try(if(all(covs %in% names(private$meta_data)) == FALSE) stop("Error: can not find all covariates from the meta data!"))
                    }
                    
                    start_time <- proc.time()
                    pseudoBulk <- private$pseudoBulk_counts
                    
                    if (private$mode == "group") {
                        sample_ids <- data.frame(sampleID = colnames(pseudoBulk), stringsAsFactors = FALSE)
                        names(sample_ids) = private$sampleId_col
                        meta_data <- unique(private$meta_data[, c(private$sampleId_col, '.GrouP', covs)])
                    } else {
                        sample_ids <- private$pseudoBulk_meta
                        meta_data <- unique(private$meta_data[, c(private$sampleId_col, covs)])
                    }
                    
                    meta_data <- sample_ids %>%
                        dplyr:::left_join(meta_data)
                    
                    other_group <- setdiff(unique(meta_data[, '.GrouP']), c('ref', 'alt'))
                    meta_data[, '.GrouP'] <- factor(meta_data[, '.GrouP'], 
                                                    levels = c('ref', 'alt', other_group))
                    meta_data[, '.GrouP'] <- relevel(meta_data[, '.GrouP'], ref = "ref")
                    
                    if (private$mode == "group") {                                                 # add by Grace
                        fml <- paste0("~ ", paste(c('.GrouP', covs), collapse = " + "))
                    } else{
                        fml <- paste0("~ ", paste(c('.GrouP', private$sampleId_col), collapse = " + ")) 
                    }
                    design = as.formula(fml)
                    dds = DESeqDataSetFromMatrix(countData = pseudoBulk, colData = meta_data, design = design)
                    dds = DESeq(dds, parallel = parallel)
                    
                    if (shrink)
                    {
                        res = lfcShrink(dds, coef = paste('.GrouP', "alt", "vs", "ref", sep = "_"),
                                        type = "apeglm", parallel = parallel)
                    } else {
                        res = results(dds, contrast = c('.GrouP', "alt", "ref"), parallel=parallel)
                    }
                    
                    res_frame = data.frame(res)
                    res_frame$ID = row.names(res_frame)
                    row.names(res_frame) = NULL
                    
                    res_frame = res_frame %>%
                        dplyr::rename("log2FC" = "log2FoldChange",
                                      "Pvalue" = "pvalue",
                                      "FDR" = "padj")
                    
                    if(shrink)
                    {
                        res_frame = res_frame[,c(6,1:5)]
                        private$de_results[["DESeq2.shrink"]] <- res_frame
                        private$de_method = "DESeq2.shrink"
                    } else {
                        res_frame = res_frame[,c(7,1:6)]
                        private$de_results[["DESeq2"]] <- res_frame
                        private$de_method = "DESeq2"
                    }
                    
                    print(proc.time() - start_time)
                    return(res_frame %>% dplyr:::arrange(FDR))
                },
                # run nebula pipeline 
                nebula_pipeline = function(covs = NULL, method = "HL") { ##LP
                    stopifnot(is.null(covs) | is.character(covs))
                    stopifnot(is.character(method) & method %in% c("LN", "HL"))
                    
                    libsizes <- private$lib_sizes
                    meta_data <- private$meta_data
                    counts <- private$counts
                    
                    # Nebula assumes that "the cells of the same subject have to be grouped!!!"
                    order_ids <- order(meta_data[, private$sampleId_col])
                    libsizes <- libsizes[order_ids]
                    meta_data <- meta_data[order_ids, ]
                    counts <- counts[ ,order_ids]
                    
                    if (!is.null(covs)) {
                        try(if(all(covs %in% names(meta_data)) == FALSE) stop("Error: can not find all covariates from the meta data!"))
                    }
                    
                    placebo_index <- which(meta_data[['.GrouP']] %in% 'ref')
                    treated_index <- which(meta_data[['.GrouP']] %in% 'alt')
                    col_index = c(placebo_index, treated_index)
                    
                    other_group <- setdiff(unique(meta_data[, '.GrouP']), c('ref', 'alt'))
                    meta_data[, '.GrouP'] <- factor(meta_data[, '.GrouP'],
                                                    levels = c('ref', 'alt', other_group))
                    meta_data[, '.GrouP'] <- relevel(meta_data[, '.GrouP'], ref = "ref")
                    
                    meta_data <- meta_data[col_index, c(private$sampleId_col, '.GrouP', covs)]
                    
                    fml = paste0(" ~ ", paste(c('.GrouP',covs),collapse = " + "))  ##LP
                    df = model.matrix(as.formula(fml),data=meta_data)
                    
                    data_g = group_cell(count=counts[, col_index],id=meta_data[, private$sampleId_col],pred=df)
                    if(is.null(data_g)){
                        re = nebula(counts[, col_index], meta_data[, private$sampleId_col], pred = df, offset = libsizes[col_index], method = method) ##LP
                    }else{
                        re = nebula(data_g$count,data_g$id,pred=data_g$pred, offset = libsizes[colnames(data_g$count)], method = method)  ##LP
                    }
                    final_table = data.frame("ID" = re$summary[,"gene"],re$summary[,grep("GrouP",colnames(re$summary))])
                    colnames(final_table) = c("ID","logFC","se","Pvalue")
                    
                    if(method == "LN")
                    {
                        final_table = cbind(final_table, algorithm = re$algorithm)
                    } else { 
                        final_table = cbind(final_table, algorithm = "HL") 
                    }
                    
                    final_table$log2FC = log2(exp(final_table$logFC))
                    final_table$FDR = p.adjust(final_table[,"Pvalue"],method="fdr")
                    
                    #LP: adding reciprocals of cell-level overdispersions 
                    cell.over <- data.frame("ID"=re$summary$gene, "cell"=re$overdispersion$Cell, "convergence"=re$convergence)
                    final_table <- final_table %>% left_join(cell.over, by="ID")
                    
                    #LP: calculating CPS
                    cps <- mean(table(meta_data[, private$sampleId_col]))
                    final_table$metric <- cps/final_table$cell
                    
                    private$de_results[["NEBULA"]] <- final_table
                    private$de_method = "NEBULA"
                    return(list(res.tab=final_table[,c(1,6,4,7,5,9,10)] %>% dplyr::arrange(FDR), res.ls=re))  ##LP
                },
                # run glmmTMB pipeline
                # return sorted results
                # by default 4 cores will be used
                glmmTMB_pipeline = function(covs = NULL, family = "nbinom2", cores = 4, detection_rate = FALSE, input_lib_sizes = NULL) {
                    stopifnot(family %in% c("nbinom2", "nbinom1", "poisson", "nbinom2zi", "nbinom1zi"))
                    stopifnot(is.null(covs) | is.character(covs))
                    stopifnot(cores%%1 == 0 & cores >= 1)
                    stopifnot(is.logical(detection_rate))
                    
                    # Number of cores to use 
                    registerDoMC(cores = cores)
                    
                    genes <- private$genes
                    libsizes <- private$lib_sizes
                    meta_data <- private$meta_data
                    
                    stopifnot(is.null(input_lib_sizes) | (length(input_lib_sizes) == nrow(meta_data)))
                    
                    if (!is.null(covs)) {
                        try(if(all(covs %in% names(meta_data)) == FALSE) stop("Error: can not find all covariates from the meta data!"))
                    }
                    
                    placebo_index <- which(meta_data[['.GrouP']] %in% 'ref')
                    treated_index <- which(meta_data[['.GrouP']] %in% 'alt')
                    col_index = c(placebo_index, treated_index)
                    ngenes <- length(genes)
                    
                    other_group <- setdiff(unique(meta_data[, '.GrouP']), c('ref', 'alt'))
                    meta_data[, '.GrouP'] <- factor(meta_data[, '.GrouP'], 
                                                    levels = c('ref', 'alt', other_group))
                    meta_data[, '.GrouP'] <- relevel(meta_data[, '.GrouP'], ref = "ref")
                    
                    meta_data <- meta_data[col_index, c(private$sampleId_col, '.GrouP', covs)]
                    
                    if(detection_rate) {
                        cdr <- scale(colSums(private$counts>0)) %>% as.data.frame %>% rename(cdr = V1)
                    }
                    
                    if(detection_rate) {
                        fml <- paste0("counts ~ ",
                                      paste(c('.GrouP', covs), collapse = " + "),
                                      glue(" + (1|{private$sampleId_col}) + offset(loglib) + cdr"))
                    } else {
                        fml <- paste0("counts ~ ",
                                      paste(c('.GrouP', covs), collapse = " + "),
                                      glue(" + (1|{private$sampleId_col}) + offset(loglib)"))
                    }
                    group_label <- paste0('.GrouP', "alt")
                    
                    analyze_gene <- function(gene)
                    {
                        counts_frame = data.frame(counts = as.vector(t(private$counts[gene, col_index])),
                                                  loglib = log(private$lib_sizes[col_index]))
                        if(detection_rate) {
                            glmer_frame = dplyr:::bind_cols(meta_data, cdr, counts_frame)
                        } else {
                            glmer_frame = dplyr:::bind_cols(meta_data, counts_frame)
                        }
                        
                        gene_id <- private$genes[gene]
                        p_val <- NA
                        log2FC <- NA
                        fail2converge <- FALSE 
                        warning <- NA
                        my_AIC <- NA
                        my_BIC <- NA
                        
                        if (family %in% c("nbinom2", "nbinom1", "poisson")) 
                        {
                            tryCatch({
                                glmm_NegBin = glmmTMB(as.formula(fml),
                                                      data  = glmer_frame,
                                                      family = family,
                                                      REML = FALSE,
                                                      control = glmmTMBControl(optCtrl = list(iter.max=1e5,eval.max=1e5)))
                                
                                if(glmm_NegBin$fit$convergence == 0) 
                                {
                                    p_val <- summary(glmm_NegBin)$coefficients$cond[group_label, 4]
                                    log2FC <- log2(exp(summary(glmm_NegBin)$coefficients$cond[group_label, 1]))
                                    my_AIC <- AIC(glmm_NegBin)
                                    my_BIC <- BIC(glmm_NegBin)
                                } else {
                                    p_val <- NA
                                    log2FC <- NA
                                    fail2converge <- TRUE
                                    
                                }
                            },warning = function(war){
                                fail2converge <<- FALSE 
                                warning <<- conditionMessage(war)
                                print(war)
                            },error = function(err){
                                fail2converge <<- FALSE 
                                print(err)
                            },finally = {
                                if (exists("glmm_NegBin")) rm(glmm_NegBin)
                            })
                        }
                        
                        if (family %in% "nbinom2zi")
                        {
                            tryCatch({
                                zifml <- as.formula(paste0("~ ", '.GrouP'))
                                glmm_NegBin = glmmTMB(as.formula(fml),
                                                      data  = glmer_frame,
                                                      family = nbinom2,
                                                      REML = FALSE,
                                                      zi = zifml, 
                                                      control = glmmTMBControl(optCtrl =list(iter.max=1e5,eval.max=1e5)))
                                
                                
                                if(glmm_NegBin$fit$convergence == 0) 
                                {
                                    p_val <- summary(glmm_NegBin)$coefficients$cond[group_label, 4]
                                    log2FC <- log2(exp(summary(glmm_NegBin)$coefficients$cond[group_label, 1]))
                                    my_AIC <- AIC(glmm_NegBin)
                                    my_BIC <- BIC(glmm_NegBin)
                                } else {
                                    p_val <- NA
                                    log2FC <- NA
                                    fail2converge <- TRUE
                                }
                                
                            },warning = function(war){
                                fail2converge <<- FALSE 
                                warning <<- conditionMessage(war)
                                print(war)
                            },error = function(err){
                                fail2converge <<- FALSE 
                            },finally = {
                                if (exists("glmm_NegBin")) rm(glmm_NegBin)
                            })
                        }
                        
                        if (family %in% "nbinom1zi")
                        {
                            tryCatch({
                                zifml <- as.formula(paste0("~ ", '.GrouP'))
                                glmm_NegBin = glmmTMB(as.formula(fml),
                                                      data  = glmer_frame,
                                                      family = nbinom1,
                                                      REML = FALSE,
                                                      zi = zifml,
                                                      control = glmmTMBControl(optCtrl =list(iter.max=1e5,eval.max=1e5)))
                                
                                warning <- names(summary(warnings(glmm_NegBin)))[1]
                                
                                if(glmm_NegBin$fit$convergence == 0) 
                                {
                                    p_val <- summary(glmm_NegBin)$coefficients$cond[group_label, 4]
                                    log2FC <- log2(exp(summary(glmm_NegBin)$coefficients$cond[group_label, 1]))
                                } else{
                                    p_val <- NA
                                    log2FC <- NA
                                    fail2converge <- TRUE
                                }
                                
                            },warning = function(war){
                                fail2converge <<- FALSE 
                                warning <<- conditionMessage(war)
                                print(war)
                            },error = function(err){
                                fail2converge <<- FALSE 
                            },finally = {
                                if (exists("glmm_NegBin")) rm(glmm_NegBin)
                            })
                        }
                        message("Gene #", gene, " is done.")
                        data.frame(ID = gene_id, log2FC = log2FC, Pvalue = p_val, Error = fail2converge, Warning = warning,AICvalue = my_AIC, BICvalue = my_BIC)
                    }
                    
                    glmm_table = foreach(gene=seq(ngenes),.combine=rbind) %dopar% {analyze_gene(gene)}
                    glmm_table$FDR = p.adjust(glmm_table$Pvalue, method="fdr")
                    
                    private$de_results[["glmmTMB"]] <- glmm_table
                    private$de_method = "glmmTMB"
                    nTotal <- nrow(glmm_table)
                    nError <- sum(glmm_table$Error, na.rm = T)
                    percError <- paste0(round(nError/nTotal*100, 2), "%")
                    message(paste0("Performed DE analysis for ",  nTotal, " genes."))
                    message(paste0("Model failed to converge for ",  nError, " (", percError, ") genes!"))
                    return(glmm_table %>% dplyr:::arrange(FDR))
                },
                # run MAST pipeline
                # return sorted results
                MAST_pipeline = function(covs = NULL, method = "glmer", ebayes = FALSE, detection_rate = TRUE) {
                    stopifnot(method %in% c("glm", "glmer", "bayesglm"))
                    stopifnot(is.null(covs) | is.character(covs))
                    stopifnot(is.logical(ebayes))
                    stopifnot(is.logical(detection_rate))
                    
                    genes <- private$genes
                    meta_data <- private$meta_data
                    
                    if (!is.null(covs)) {
                        try(if(all(covs %in% names(meta_data)) == FALSE) stop("Error: can not find all covariates from the meta data!"))
                    }
                    
                    sce <- SingleCellExperiment(assays = list(counts=as.matrix(private$counts)),
                                                colData = meta_data)
                    norm.sce <- normalize(sce)
                    sca <- SceToSingleCellAssay(norm.sce)
                    
                    if(detection_rate) {
                        cdr <- as.data.frame(colSums(as.matrix(private$counts)>0)/nrow(private$counts)*100) 
                        colnames(cdr)<-"cdr"
                        colData(sca)$cdr <- cdr$cdr
                    }
                    
                    if(detection_rate) {
                        fml <- paste0(" ~ ",
                                      paste(c('.GrouP', covs), collapse = " + "),
                                      glue(" + (1|{private$sampleId_col}) + cdr"))
                    } else {
                        fml <- paste0(" ~ ",
                                      paste(c('.GrouP', covs), collapse = " + "),
                                      glue(" + (1|{private$sampleId_col})"))
                    }
                    grp <- factor(colData(sca)$.GrouP)
                    grp <- relevel(grp, "ref")
                    colData(sca)$.GrouP <- grp
                    
                    zlm.fit <- zlm(as.formula(fml), sca = sca, method = method, ebayes = ebayes)
                    summaryGrp <- summary(zlm.fit, doLRT = '.GrouPalt')
                    summaryDt <- summaryGrp$datatable
                    fcHurdle <- merge(summaryDt[contrast=='.GrouPalt' & component=='H',.(ID = primerid, Pvalue = `Pr(>Chisq)`)], #hurdle P values
                                      summaryDt[contrast=='.GrouPalt' & component=='logFC', .(ID = primerid, log2FC = coef, ci.hi, ci.lo)], by='ID') #logFC coefficients
                    
                    fcHurdle[, FDR:=p.adjust(Pvalue, 'fdr')]
                    setorder(fcHurdle, FDR)
                    private$de_results[["MAST"]] <- as.data.frame(fcHurdle)
                    private$de_method = "MAST"
                    return(as.data.frame(fcHurdle))
                },
                # function for the user to get violin plot
                # This function can takes either private R6 objects or user provided counts and meta data
                violinPlot = function(countsdata = NULL, meta = NULL, cellinfo = NULL, xgroup = NULL, colorgroup = NULL, gene.name = "", de.method = "glmmTMB", batch = NULL)
                {
                    edb <- EnsDb.Hsapiens.v86 
                    # use private R6 objects
                    if(is.null(countsdata) & is.null(meta) & is.null(xgroup) & is.null(colorgroup))
                    {
                        stopifnot(is.character(gene.name) & gene.name != "")
                        stopifnot(is.null(batch) | is.character(batch))
                        stopifnot(is.null(de.method) | is.character(de.method))
                        stopifnot(is.character(cellinfo))
                        
                        if(is.null(de.method)) {
                            de.method = private$de_method
                        } else {
                            if(!(de.method %in% names(private$de_results))) stop(paste0("DE results not availabe for the '", de.method, "' method. Please make sure to run the DE analysis using that method before generating the violin plot"))
                        }
                        
                        treatment = private$treatment_col
                        meta_data <- private$meta_data
                        genes <- private$genes
                        
                        stopifnot(gene.name %in% genes)
                        
                        if(is.data.frame(private$norm_counts)){
                            plotdat <- as.data.frame(t(private$norm_counts[rownames(private$norm_counts)==gene.name, ]))%>%rownames_to_column(cellinfo)
                        } else {
                            plotdat <- as.data.frame(private$norm_counts[rownames(private$norm_counts)==gene.name, ])%>%rownames_to_column(cellinfo)
                        }
                        
                        colnames(plotdat)[2] <- "Normalized_Count"
                        plotdat <- plotdat %>%
                            dplyr::mutate(logcount = log2(Normalized_Count + 1)) %>%
                            dplyr:::left_join(meta_data)
                        
                        # get biotype for gene
                        goids <- genes(edb, columns = c("gene_name","gene_biotype"), filter = GeneNameFilter(gene.name))
                        
                        # get FDR log2FC and pvalue
                        fcdat <- private$de_results[[de.method]] %>% dplyr::filter(ID == gene.name)
                        fcdat$log2FC <- round(fcdat$log2FC, digits=2)
                        fcdat$FDR <- formatC(fcdat$FDR, format = "e", digits = 2)
                        fcdat$Pvalue <- formatC(fcdat$Pvalue, format = "e", digits = 2)
                        if(is.null(batch)) {
                            p <- ggplot(plotdat, aes_string(private$sampleId_col, "logcount")) +
                                geom_violin(aes_string(fill = treatment), alpha = 0.2) +
                                geom_jitter(position = position_jitter(0.2),size=0.1) +
                                stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                             geom = "crossbar", 
                                             width = 0.5,
                                             position = position_dodge(width = .25))+
                                theme_bw() +
                                labs(y = "log2(Normalized Counts + 1)", title = paste0(gene.name, ", ", goids$gene_biotype,", ",de.method, ", log2FC = ", fcdat$log2FC,", FDR = ", fcdat$FDR, ", Pvalue = ", fcdat$Pvalue))+
                                theme(plot.title = element_text(size = 10))
                        } else {
                            stopifnot(batch %in% colnames(plotdat))
                            p <- ggplot(plotdat, aes_string(private$sampleId_col, "logcount")) +
                                geom_violin(aes_string(fill = treatment), alpha = 0.2) +
                                geom_jitter(aes_string(colour = batch), position = position_jitter(0.2),size=0.1) +
                                stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                             geom = "crossbar", 
                                             width = 0.5,
                                             position = position_dodge(width = .25))+
                                theme_bw() +
                                labs(y = "log2(Normalized Counts + 1)", title = paste0(gene.name, ", ",goids$gene_biotype,", ", de.method, ", log2FC = ", fcdat$log2FC,", FDR = ", fcdat$FDR, ", Pvalue = ", fcdat$Pvalue))
                        }
                        return(p)
                    } else {
                        if(is.null(countsdata)) stop ("Please provide counts data")
                        if(is.null(meta)) stop ("Please provide meta data")
                        
                        stopifnot(is.character(cellinfo))
                        stopifnot(is.character(xgroup))
                        stopifnot(is.character(colorgroup))
                        stopifnot(is.character(gene.name) & gene.name != "")
                        stopifnot(is.null(batch) | is.character(batch))
                        stopifnot(cellinfo %in% colnames(meta))
                        stopifnot(xgroup %in% colnames(meta))
                        stopifnot(colorgroup %in% colnames(meta))
                        stopifnot(identical(colnames(countsdata),meta[[cellinfo]])) 
                        
                        colorby <-factor(meta[[colorgroup]])
                        group<-factor(meta[[xgroup]])
                        genes <- rownames(countsdata)
                        stopifnot(gene.name %in% genes)
                        
                        if(is.data.frame(countsdata)){
                            plotdat <- as.data.frame(t(countsdata[rownames(countsdata)==gene.name, ]))%>%rownames_to_column(cellinfo)
                        } else {
                            plotdat <- as.data.frame(countsdata[rownames(countsdata)==gene.name, ])%>%rownames_to_column(cellinfo)
                        }
                        
                        colnames(plotdat)[2] <- "Normalized_Count"
                        plotdat <- plotdat %>%
                            dplyr::mutate(logcount= log2(Normalized_Count+1)) %>%
                            dplyr:::left_join(meta)
                        # get biotype for gene  
                        # ensembl <-useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
                        # goids = getBM(attributes = c('ensembl_gene_id', 'gene_biotype','external_gene_name'), 
                        #               filters = 'external_gene_name', 
                        #               values = gene.name, 
                        #               mart = ensembl)
                        # if (!requireNamespace("BiocManager", quietly = TRUE))
                        #   install.packages("BiocManager")
                        
                        goids<-genes(edb,columns=c("gene_name","gene_biotype"),filter=GeneNameFilter(gene.name))
                        if(is.null(batch)) {
                            p <- ggplot(plotdat, aes_string(xgroup, "logcount")) +
                                geom_violin(aes_string(fill = colorby), alpha = 0.2) +
                                geom_jitter(position = position_jitter(0.2),size=0.1) +
                                stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                             geom = "crossbar", 
                                             width = 0.5,
                                             position = position_dodge(width = .25))+
                                theme_bw() +
                                labs(y = "log2(Normalized Counts + 1)", title = paste0(gene.name,", ",goids$gene_biotype))
                        } else {
                            
                            stopifnot(batch %in% colnames(meta))
                            p <- ggplot(plotdat, aes_string(xgroup, "logcount")) +
                                geom_violin(aes_string(fill = colorby), alpha = 0.2) +
                                geom_jitter(aes_string(colour = batch), position = position_jitter(0.2),size=0.1) +
                                stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                             geom = "crossbar", 
                                             width = 0.5,
                                             position = position_dodge(width = .25))+
                                theme_bw() +
                                labs(y = "log2(Normalized Counts + 1)", title = paste0(gene.name,", ",goids$gene_biotype))
                        }
                        return(p)
                    }
                },
                # volcano plot
                volcanoPlot = function(FDR_threshold = 0.05, FC_threshold = 1, de.method = NULL, title = "Volcano Plot")
                {
                    stopifnot(is.numeric(FDR_threshold) & FDR_threshold >= 0 & FDR_threshold <= 1)
                    stopifnot(is.numeric(FC_threshold) & FC_threshold >= 1)
                    stopifnot(is.null(title) | is.character(title))
                    stopifnot(is.null(de.method) | is.character(de.method))
                    
                    if(is.null(de.method)) {
                        de.method = private$de_method
                    } else {
                        if(!(de.method %in% names(private$de_results))) stop(paste0("DE results not availabe for the '", de.method, "' method. Please make sure to run the DE analysis using that method before generate the violin plot"))
                    }
                    
                    if (FC_threshold == 1) {
                        subtitle = paste0("Significant: FDR < ", FDR_threshold, "\n",
                                          "Not Significant: FDR >= ", FDR_threshold, "\n")
                        res <- private$de_results[[de.method]] %>%
                            dplyr::mutate(note =
                                              case_when(
                                                  FDR < FDR_threshold ~ "Significant",
                                                  TRUE ~ "Not significant"
                                              )) %>%
                            dplyr::mutate(note = factor(note, levels = c("Significant", "Not significant")))
                        cols <- c("Significant" = "#F8766D", "Not significant" = "#619CFF")
                        p <- ggplot(res, aes(log2FC, -log10(Pvalue))) +
                            geom_point(aes(color = note)) +
                            labs(title = title, subtitle = subtitle) +
                            scale_colour_manual(values = cols) +
                            theme_bw() +
                            theme(legend.title = element_blank())
                    } else {
                        log2fc_threshold = log2(FC_threshold)
                        subtitle = paste0("Up-regulated: FDR < ", FDR_threshold, " & log2FC > ", round(log2fc_threshold, 3), "\n",
                                          "Down-regulated: FDR < ", FDR_threshold, " & log2FC < ", -round(log2fc_threshold, 3), "\n")
                        res <- private$de_results[[de.method]] %>%
                            dplyr::mutate(note =
                                              case_when(
                                                  FDR < FDR_threshold & log2FC > log2fc_threshold ~ "Up-regulated",
                                                  FDR < FDR_threshold & log2FC < -log2fc_threshold ~ "Down-regulated",
                                                  TRUE ~ "Not significant"
                                              )) %>%
                            dplyr::mutate(note = factor(note, levels = c("Up-regulated", "Down-regulated", "Not significant")))
                        cols <- c("Up-regulated" = "#F8766D", "Down-regulated" = "#00BA38", "Not significant" = "#619CFF")
                        p <- ggplot(res, aes(log2FC, -log10(Pvalue))) +
                            geom_point(aes(color = note)) +
                            geom_vline(xintercept = log2fc_threshold, linetype = "dotted") +
                            geom_vline(xintercept = -log2fc_threshold, linetype = "dotted") +
                            labs(title = title, subtitle = subtitle) +
                            scale_colour_manual(values = cols) +
                            theme_bw() +
                            theme(legend.title = element_blank())
                    }
                    return(p)
                },
                # getter functions
                assayData = function()
                {
                    return(private$counts)
                },
                pseudoBulkData = function()
                {
                    return(private$pseudoBulk_counts)
                },
                pData = function()
                {
                    return(private$meta_data)
                },
                get_cells = function()
                {
                    return(private$cells)
                },
                get_genes = function()
                {
                    return(private$genes)
                },
                get_MT = function()
                {
                    return(private$MT_rows)
                },
                get_lib_sizes = function()
                {
                    return(private$lib_sizes)
                },
                get_av_lib_size = function()
                {
                    return(private$av_lib)
                },
                get_norm_counts = function()
                {
                    return(private$norm_counts)
                },
                get_DE_results = function()
                {
                    return(private$de_results)
                },
                get_filter_info = function()          
                {
                    return(private$filter_info)
                }
            )
    )

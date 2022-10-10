#!/usr/bin/env Rscript
#args = commandArgs(trailingOnly = TRUE)

suppressMessages(suppressWarnings(require(ggplot2)))
suppressMessages(suppressWarnings(require(readr)))
suppressMessages(suppressWarnings(require(SummarizedExperiment)))
suppressMessages(suppressWarnings(require(SingleCellExperiment)))
suppressMessages(suppressWarnings(require(Seurat)))
suppressMessages(suppressWarnings(require(tidyr)))

#Big Block of Functions

PlotSmoothers <- function(models,  
                                      gene,
                                      Xcolnames, 
                                      nPoints = 100, 
                                      lwd = 2,
                                      size = 2/3,
                                      xlab = "Pseudotime",
                                      ylab = "Log(expression + 1)",
                                      border = FALSE,
                                      sample = 1,
                                      alpha = 2/3,
                                      pointCol = NULL,
                                      curvesCols = NULL,
                                      plotLineages = TRUE)
{

  counts(models) <- assay(models, "X")
  counts = counts(models)
  
  if (is.null(names(models))) {
    rownames(models) <- rownames(counts) <- seq_len(nrow(models))
    #print("in the is.null loop for some reason")
    message(paste0(
      "The sce object has no rownames. Assuming that the counts and the sce ",
      "objects are ordered in the same way"))
  }
  # input is singleCellExperiment object.
  if(length(gene) > 1) stop("Only provide a single gene's ID with the ",
                            "gene argument.")
  # check if all gene IDs provided are present in the models object.
  if (is(gene, "character")) {

    if (!all(gene %in% names(models))) {
      stop("The gene ID is not present in the models object.")
    }
   
    id = gene
  } else id <- gene
  
  #Create dm, X slingshotColData

  df = colData(models)

  dm = df[ , grepl( "dm." , names( df ) ) ]
  dm= as.data.frame(dm)
  colnames(dm) = gsub(pattern = "dm.", replacement = "", x = colnames(dm))

  X = df[ , grepl( "X." , names( df ) ) ]
  X = as.data.frame(X)

  colnames(X) = Xcols
              
  toMatch <- c("pseudotime.", "cellWeights.")
  
  slingshotColData = df[ , grepl((paste(toMatch,collapse="|")), names( df )) ]
  
  slingshotColData$pseudotime_1 = NULL
  slingshotColData$pseudotime_2 = NULL
  
  conditions = df$conditions

   #Construct Beta
  
  df2 = rowData(models)
  
  all_beta = df2[ , grepl( "beta." , names( df2 ) ) ]
  beta = all_beta[id,]
  beta = as.data.frame(beta)

  #dm <- colData(models)$tradeSeq$dm # design matrix
  #y <- unname(counts[names(models),][id,]) 
  y <- unname(counts[id,])
  y = round(y)

  #X <- colData(models)$tradeSeq$X # linear predictor
  #slingshotColData <- colData(models)$slingshot
  pseudotime <- slingshotColData[,grep(x = colnames(slingshotColData),
                                       pattern = "pseudotime")]
  nLineages <- length(grep(x = colnames(dm), pattern = "t[1-9]"))
  #beta <- rowData(models)$tradeSeq[id, ]$beta[[1]]
  if(any(is.na(beta))){
    stop("Some coefficients for this gene are NA. Cannot plot this gene.")
  }
  #conditions <- colData(models)$tradeSeq$conditions
  nConditions <- nlevels(conditions)
  
  # Construct time variable based on cell assignments.
  lcol <- timeAll <- rep(0, nrow(dm))
  # Loop over all curves, i.e. all lineqges and conditions
  for (jj in seq_len(nLineages)) {
    for (kk in seq_len(nConditions)){
      for (ii in seq_len(nrow(dm))) {
        if (dm[ii, paste0("l", jj, "_", kk)] == 1) {
          timeAll[ii] <- dm[ii, paste0("t", jj)]
          lcol[ii] <- paste0("Lineage ", jj, "_", levels(conditions)[kk])
        } else {
          next
        }
      }
    }
  }
  
  if (!is.null(pointCol)) {
    if (length(pointCol) == 1) {
      col <- colData(models)[, pointCol]
    } else if (length(pointCol) == ncol(models)) {
      col <- pointCol
    } else {
      col <- lcol
      message(paste("pointCol should have length of either 1 or the number of cells,",
                    "reverting to default color scheme, by lineages and conditions."))
    }
  } else {
    col <- lcol
  }
  

  # plot raw data
  df <- data.frame("time" = timeAll,
                   "gene_count" = y,
                   "pCol" = as.character(col),
                   "lineage" = as.character(lcol))
  
  
  
  # Reorder curves according to the levels of conditions
  combs <- paste0("Lineage ", seq_len(nLineages), "_")
  combs <- lapply(combs, function(lin) paste0(lin, levels(conditions)))
  combs <- do.call('c', combs)
  df$lineage <- factor(df$lineage, levels = combs)
  rows <- sample(seq_len(nrow(df)), nrow(df) * sample, replace = FALSE)
  df <- df[rows, ]

  # Create Basic Plot

  p <- ggplot(df, aes(x = time, y = log1p(gene_count))) +
    labs(x = xlab, y = ylab) +
    theme_classic()
  if(is.null(pointCol)){
    p <- p +
      geom_point(size = size, aes(col = lineage)) +
      scale_color_viridis_d(alpha = alpha)
  } else {
    p <- p +
      geom_point(size = size, alpha = alpha, aes(col = pCol)) +
      scale_color_discrete() +
      labs(col = "Cell labels")
  }
  

  # Predict and plot smoothers across the range

  if (plotLineages) {
    if (!is.null(curvesCols)) {
      if (length(curvesCols) != nLineages * nConditions) {
        curvesCols <- viridis::viridis(nLineages * nConditions)
        message("Incorrect number of lineage colors. Default to viridis")
      }
    } else {
      curvesCols <- viridis::viridis(nLineages * nConditions)
    }
    
    for (jj in seq_len(nLineages)) {
      for(kk in seq_len(nConditions)){
        df <- .getPredictRangeDf(dm, lineageId = jj, conditionId = kk,
                                 nPoints = nPoints)
        Xdf <- predictGAM(lpmatrix = X,
                          df = df,
                          pseudotime = pseudotime,
                          conditions = conditions)
        yhat <-  c(exp(t(Xdf %*% t(beta)) + df$offset))
        if (border) {
          p <- p +
            geom_line(data = data.frame("time" = df[, paste0("t", jj)],
                                        "gene_count" = yhat,
                                        "lineage" = as.character(paste0(jj, "_", kk))),
                      lwd = lwd + 1, colour = "white")
          
        }
        p <- p +
          geom_line(data = data.frame("time" = df[, paste0("t", jj)],
                                      "gene_count" = yhat,
                                      "lineage" = as.character(paste0(jj, "_", kk))),
                    lwd = lwd,
                    col = curvesCols[jj * nConditions - (nConditions - kk)])
      }
    }
  }

  return(p)
}


## get predictor matrix for a range of pseudotimes of a smoother.
.getPredictRangeDf <- function(dm, lineageId, conditionId = NULL, nPoints = 100){
  vars <- dm[1, ]
  if ("y" %in% colnames(vars)) {
    vars <- vars[!colnames(vars) %in% "y"]
    off <- 1
  } else {
    off <- 0
  }
  offsetId <- grep(x = colnames(vars), pattern = "offset")
  offsetName <- colnames(vars)[offsetId]
  offsetName <- substr(offsetName, start = 8, stop = nchar(offsetName) - 1)
  names(vars)[offsetId] <- offsetName
  # set all times on 0
  vars[, grep(colnames(vars), pattern = "t[1-9]")] <- 0
  # set all lineages on 0
  vars[, grep(colnames(vars), pattern = "l[1-9]")] <- 0
  # duplicate to nPoints
  vars <- rbind(vars, vars[rep(1, nPoints - 1), ])
  # set range of pseudotime for lineage of interest
  if (is.null(conditionId)) {
    lineageIds <- grep(colnames(vars), pattern = paste0("l", lineageId, "($|_)"))
  } else {
    lineageIds <- grep(colnames(vars), pattern = paste0("l", lineageId,
                                                        "_", conditionId, "$"))
  }
  if (length(lineageIds) == 1){
    lineageData <- dm[dm[, lineageIds + off] == 1,
                      paste0("t", lineageId)]
  } else {
    lineageData <- dm[rowSums(dm[, lineageIds + off]) == 1,
                      paste0("t", lineageId)]
  }
  # make sure lineage starts at zero
  if(min(lineageData) / max(lineageData) < .01) {
    lineageData[which.min(lineageData)] <- 0
  }
  vars[, lineageIds] <- 1 / length(lineageIds)
  # set lineage
  vars[, paste0("t", lineageId)] <- seq(min(lineageData),
                                        max(lineageData),
                                        length = nPoints)
  # set offset
  vars[, offsetName] <- mean(dm[, grep(x = colnames(dm),
                                       pattern = "offset")])
  return(vars)
}


# Predicting fits ----
# lpmatrix given X and design
predictGAM <- function(lpmatrix, df, pseudotime, conditions = NULL){
  # this function is an alternative of predict.gam(model, newdata = df, type = "lpmatrix")
  # INPUT:
  # lpmatrix is the linear predictor matrix of the GAM model
  # df is a data frame of values for which we want the lpmatrix
  # pseudotime is the n x l matrix of pseudotimes
  # conditions is the vector of conditions, if present.
  
  # if pseudotime is vector, make it a matrix.
  if(is.null(dim(pseudotime))) pseudotime <- matrix(pseudotime,ncol=1)
  
  condPresent <- !is.null(conditions)
  if(condPresent) nConditions <- nlevels(conditions)
  
  # for each curve, specify basis function IDs for lpmatrix
  allBs <- grep(x = colnames(lpmatrix), pattern = "[0-9]):l[1-9]")
  
  if(!condPresent){
    lineages <- sub(pattern = "s\\(", replacement = "",
                    x = colnames(lpmatrix[,allBs]))
    lineages <- sub(pattern = "\\):.*", replacement = "",
                    x = lineages)
    nCurves <- length(unique(lineages))
    for (ii in seq_len(nCurves)) {
      assign(paste0("id",ii), allBs[which(lineages == paste0("t", ii))])
    }
  } else if(condPresent){
    lineages <- sub(pattern = "s\\(t", replacement = "",
                    x = colnames(lpmatrix[,allBs]))
    lineages <- sub(pattern = "\\):.*", replacement = "",
                    x = lineages)
    nLineages <- length(unique(lineages))
    curves <- sub(pattern = ".*:l", replacement = "",
                  x = colnames(lpmatrix[,allBs]))
    curves <- sub(pattern = "\\..*", replacement = "",
                  x = curves)
    nCurves <- length(unique(curves))
    for (ii in seq_len(nLineages)) {
      for(kk in seq_len(nConditions))
        assign(paste0("id", ii, "_", kk), allBs[which(curves == paste0(ii, "_", kk))])
    }
  }
  
  
  # specify lineage assignment for each cell (i.e., row of lpmatrix)
  if(!condPresent){
    lineageID <- apply(lpmatrix, 1, function(x){
      for (ii in seq_len(nCurves)) {
        if (!all(x[get(paste0("id", ii))] == 0)) {
          return(ii)
        }
      }
    })
  } else if(condPresent){
    # first number is lineage, second number is condition.
    lineageID <- apply(lpmatrix, 1, function(x){
      for (ii in seq_len(nLineages)) {
        # loop over lineages
        for(kk in seq_len(nConditions)){
          # loop over conditions
          if (!all(x[get(paste0("id", ii, "_", kk))] == 0)) {
            return(as.numeric(paste0(ii, kk)))
          }
        }
      }
    })
  }
  
  
  # fit splinefun for each basis function based on assigned cells
  if(!condPresent) {
    for (ii in seq_len(nCurves)) { # loop over curves
      for (jj in seq_len(length(allBs) / nCurves)) { #within curve, loop over basis functions
        assign(paste0("l",ii,".",jj),
               stats::splinefun(x = pseudotime[lineageID == ii, ii],
                                y = lpmatrix[lineageID == ii, #only cells for lineage
                                             get(paste0("id", ii))[jj]],
                                ties = mean)) #basis function
      }
    }
  } else if(condPresent) {
    for (ii in  seq_len(nLineages)) {
      # loop over curves
      for(kk in seq_len(nConditions)){
        for (jj in seq_len(length(allBs) / (nLineages * nConditions))) {
          #within curve, loop over basis functions
          assign(paste0("l",ii, "_", kk,".",jj),
                 stats::splinefun(
                   x = pseudotime[lineageID == as.numeric(paste0(ii, kk)), ii],
                   y = lpmatrix[lineageID == as.numeric(paste0(ii, kk)), #only cells for lineage
                                get(paste0("id", ii, "_", kk))[jj]],
                   ties = mean)) #basis function
        }
      }
    }
  }
  
  
  # use input to estimate X for each basis function
  Xout <- matrix(0, nrow = nrow(df), ncol = ncol(lpmatrix))
  if(!condPresent){
    for (ii in seq_len(nCurves)) { # loop over curves
      if (all(df[, paste0("l", ii)] == 1)) { # only predict if weight = 1
        for (jj in seq_len(length(allBs) / nCurves)) { # within curve, loop over basis functions
          f <- get(paste0("l", ii, ".", jj))
          Xout[, get(paste0("id", ii))[jj]] <- f(df[, paste0("t", ii)])
        }
      }
    }
  } else if(condPresent){
    # for (ii in (seq_len(nCurves)[seq(2, nCurves, by=2)])/2) {
    for (ii in seq_len(nLineages)) {
      # loop over curves
      for(kk in seq_len(nConditions)){
        # loop over conditions
        if (all(df[, paste0("l", ii, "_", kk)] != 0)) { # only predict if weight = 1
          for (jj in seq_len(length(allBs) / (nLineages * nConditions))) { 
            # within curve, loop over basis functions
            f <- get(paste0("l", ii, "_", kk, ".", jj))
            Xout[, get(paste0("id", ii, "_", kk))[jj]] <- f(df[, paste0("t", ii)])
          }
        }
      }
    }
  }
  
  
  # add fixed covariates as in df
  dfSmoothID <- grep(x = colnames(df), pattern = "[t|l][1-9]")
  dfOffsetID <- grep(x = colnames(df), pattern = "offset")
  Xout[, -allBs] <- df[, -c(dfSmoothID, dfOffsetID)]
  
  # return
  colnames(Xout) <- colnames(lpmatrix)
  return(Xout)
}

#Getting Smoother lines for each condition
predictSmoother <- function(models,gene,nPoints = 100,tidy = TRUE){

# check if all gene IDs provided are present in the models object.
if (is(gene, "character")) {
  if (!all(gene %in% rownames(models))) {
    stop("Not all gene IDs are present in the models object.")
    }
  id <- match(gene, rownames(models))
    } else id <- gene

#get dm, X, and slingshotColData
df = colData(models)

dm = df[ , grepl( "dm." , names( df ) ) ]
dm = as.data.frame(dm)
colnames(dm) = gsub(pattern = "dm.", replacement = "", x = colnames(dm))

X = df[ , grepl( "X." , names( df ) ) ]
X = as.data.frame(X)

Xcols = c("U","s(t1):l1_1.1", "s(t1):l1_1.2", "s(t1):l1_1.3", "s(t1):l1_1.4", "s(t1):l1_1.5", "s(t1):l1_1.6",
"s(t1):l1_1.7", "s(t1):l1_1.8", "s(t1):l1_1.9", "s(t1):l1_2.1", "s(t1):l1_2.2", "s(t1):l1_2.3", "s(t1):l1_2.4",
"s(t1):l1_2.5", "s(t1):l1_2.6", "s(t1):l1_2.7", "s(t1):l1_2.8" ,"s(t1):l1_2.9", "s(t2):l2_1.1", "s(t2):l2_1.2",
"s(t2):l2_1.3", "s(t2):l2_1.4", "s(t2):l2_1.5", "s(t2):l2_1.6", "s(t2):l2_1.7", "s(t2):l2_1.8", "s(t2):l2_1.9",
"s(t2):l2_2.1", "s(t2):l2_2.2", "s(t2):l2_2.3", "s(t2):l2_2.4", "s(t2):l2_2.5", "s(t2):l2_2.6", "s(t2):l2_2.7",
"s(t2):l2_2.8","s(t2):l2_2.9")

colnames(X) = Xcols
              
toMatch <- c("pseudotime.", "cellWeights.")
  
slingshotColData = df[ , grepl((paste(toMatch,collapse="|")), names( df )) ]
  
slingshotColData$pseudotime_1 = NULL
slingshotColData$pseudotime_2 = NULL
  
#Construct Beta
  
df2 = rowData(models)
  
all_beta = df2[ , grepl( "beta." , names( df2 ) ) ]
beta = all_beta[id,]
beta = as.data.frame(beta)

# get tradeSeq info
#dm <- colData(models)$tradeSeq$dm # design matrix
#X <- colData(models)$tradeSeq$X # linear predictor

#slingshotColData <- colData(models)$slingshot
pseudotime <- slingshotColData[,grep(x = colnames(slingshotColData),
                                                 pattern = "pseudotime")]
                                      
if (is.null(dim(pseudotime))) pseudotime <- matrix(pseudotime, ncol = 1)
  #betaMat <- rowData(models)$tradeSeq$beta[[1]]
  #beta <- as.matrix(betaMat[id,])
  #rownames(beta) <- gene
  df2 = rowData(models)
  
  all_beta = df2[ , grepl( "beta." , names( df2 ) ) ]
  beta = all_beta[id,]
  beta = as.data.frame(beta)

  #condPresent <- suppressWarnings({
  #!is.null(SummarizedExperiment::colData(models)$tradeSeq$conditions)
  #})

  condPresent = 1 

if(!condPresent){
  yhatMat <- .predictSmooth(dm = dm,
    X = X,
    beta = beta,
    pseudotime = pseudotime,
    gene = gene,
    nPoints = nPoints,
    tidy = tidy)
} else if(condPresent){
  #conditions <- SummarizedExperiment::colData(models)$tradeSeq$conditions
  conditions = df$conditions
  yhatMat <- .predictSmooth_conditions(dm = dm,
                                      X = X,
                                      beta = beta,
                                      pseudotime = pseudotime,
                                      gene = gene,
                                      nPoints = nPoints,
                                      conditions = conditions,
                                      tidy = tidy)
}
return(yhatMat)
  
} # end of function


.predictSmooth_conditions <- function(dm, X, beta, pseudotime, gene, nPoints, conditions, tidy){
  
  nCurves <- length(grep(x = colnames(dm), pattern = "t[1-9]"))
  nConditions <- nlevels(conditions)

  # get predictor matrix
  if (tidy) out <- list()
  for (jj in seq_len(nCurves)) {
    if (tidy) out_cond <- list()
    for(kk in seq_len(nConditions)){
      df <- .getPredictRangeDf(dm, lineageId = jj, conditionId = kk,
                               nPoints = nPoints)
      Xdf <- predictGAM(lpmatrix = X,
                        df = df,
                        pseudotime = pseudotime,
                        conditions = conditions)
      if(kk == 1) XallCond <- Xdf
      if(kk > 1) XallCond <- rbind(XallCond, Xdf)
      if (tidy) {
        out_cond[[kk]] <- data.frame(lineage = jj, time = df[, paste0("t",jj)],
                                     condition = levels(conditions)[kk])
      }
    }
    if (jj == 1) Xall <- XallCond
    if (jj > 1) Xall <- rbind(Xall, XallCond)
    if (tidy) out[[jj]] <- do.call(rbind, out_cond)
  }
  if (tidy) outAll <- do.call(rbind, out)

  # loop over all genes
  yhatMat <- matrix(NA, nrow = length(gene), ncol = nCurves * nConditions * nPoints)
  rownames(yhatMat) <- gene
  pointNames <- expand.grid(1:nCurves, 1:nConditions)
  baseNames <- paste0("lineage", pointNames[,1], "_condition",
                      levels(conditions)[pointNames[,2]])
  colnames(yhatMat) <- c(sapply(baseNames, paste0, "_point",1:nPoints))
  for (jj in 1:length(gene)) {
    yhat <- c(exp(t(Xall %*% t(beta[as.character(gene[jj]), ,
                                    drop = FALSE])) +
                    df$offset[1]))
    yhatMat[jj, ] <- yhat
  }
  ## return output
  if (!tidy) {
    return(yhatMat)
  } else {
    outList <- list()
    for (gg in seq_len(length(gene))){
      curOut <- outAll
      curOut$gene <- gene[gg]
      curOut$yhat <- yhatMat[gg,]
      outList[[gg]] <- curOut
    }
    return(do.call(rbind, outList))
  }
}

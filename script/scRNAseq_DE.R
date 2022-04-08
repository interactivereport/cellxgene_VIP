## wrapper function for pipeline_class.R
## to have a better user experience, easier usage

args <- commandArgs(trailingOnly=TRUE)
if(length(args)==0) scRNAseq_DE_path <- normalizePath(dirname(sys.frame(1)$ofile))

scRNAseq_DE <- function(
    countRDS,
    metaRDS,
    output,
    method,
    column_sample,
    column_cluster,
    column_group=NULL,
    grp_ref=NULL,
    grp_alt=NULL,
    method_model=NULL,
    column_covars=NULL,
    
    min.cells.per.gene = 3,
    min.genes.per.cell = 250,
    min.perc.cells.per.gene = 0.00,
    perc_filter = TRUE,

    R6_min.cells.per.gene = 3,
    R6_min.perc.cells.per.gene = 0.1,
    R6_min.cells.per.gene.type = "or",
    R6_cells.per.gene.filter = TRUE,
    R6_perc.cells.filter = TRUE,
    R6_perc.filter = FALSE,
    R6_perc.filter.type = "and",
    R6_perc_threshold = 0.90,
    R6_min.ave.pseudo.bulk.cpm = 1,
    R6_pseudo.bulk.cpm.filter = FALSE,
    R6_min.cells.per.subj = 3,
    
    core=1,
    qsub=FALSE,
    addSRC=NULL
){
    env <- as.list(environment())
    checkInput(env)
    meta <- readRDS(metaRDS)
    if(qsub){
        jID <- paste0("scDE",sample(100,1),"_")
        saveRDS(env,file=paste0(output,"/",jID,"env.rds"))
        sh <- c("#!/bin/bash",
                "#$ -N jID",
                "#$ -wd wkPath",
                "#$ -pe node core",
                "#$ -l h_rt=50:00:00",
                "#$ -o jID.log",
                "#$ -e jID.log",
                "#- End UGE embedded arguments",
                ": > $SGE_STDOUT_PATH",
                "cat $PE_HOSTFILE",
                "Rscript cmdPath/scRNAseq_DE.R cmdPath grpInterest")
        
        for(one in unique(meta[,column_cluster])){
            oneSh <- gsub("jID",paste0(jID,one),sh)
            oneSh <- gsub("wkPath",output,oneSh)
            oneSh <- gsub("core",core,oneSh)
            oneSh <- gsub("cmdPath",scRNAseq_DE_path,oneSh)
            oneSh <- gsub("grpInterest",paste(jID,one,addSRC),oneSh)
            strSH <- paste0(output,"/",jID,one,".sh")
            cat(paste(oneSh,collapse="\n"),"\n",sep="",file=strSH)
            system(paste("qsub",strSH))
        }

    }else{
        if(!is.null(addSRC)) for(one in unlist(strsplit(addSRC,";"))) source(one)
        source(paste0(scRNAseq_DE_path,"/pipeline_class.R"))
        for(one in unique(meta[,column_cluster])){
            print(system.time(
                scRNAseq_DE_one(env,one)
                ))
        }
    }
    
}

checkInput <- function(env){
    if(!file.exists(env$countRDS) || !file.exists(env$metaRDS)){
        stop("Either count RDS file or meta RDS file is missing!")
    }
    meta <- readRDS(env$metaRDS)
    for(one in c(env$column_sample,env$column_cluster,env$column_group,env$column_covars)){
        if(!one%in%colnames(meta))
            stop(paste(one,"is not in the sample meta table!"))
    }
    if(!is.null(env$column_group)){
        if(is.null(env$grp_ref) || is.null(env$grp_alt))stop(paste("grp_ref and grp_alt are required for",env$column_group))
        for(one in c(env$grp_ref,env$grp_alt)){
            if(!one%in%unique(meta[,env$column_group]))
                stop(paste(one,"is not in the column",env$column_group,"from sample meta table!"))
        }
    }
    allMethods <- c("t_test", "u_test","edgeR","limma","DESeq2","MAST","limma_cell_level","glmmTMB","nebula","ancova")
    if(!env$method%in%allMethods)
        stop(paste0("method (",env$method,")is not supported!\nSelect from ",
                    paste(allMethods,collapse=", ")))
    if(env$method=="nebula"){
        if(!env$method_model%in%c("LN", "HL")) stop("method_model has to be LN or HL for nebula method!")
        if(env$method_model!="HL") warning("method_model is recommended to be HL for nebula method!")
    }
    else if(env$method=="glmmTMB"){
        if(!env$method_model%in% c("nbinom2", "nbinom1", "poisson", "nbinom2zi", "nbinom1zi"))
            stop("method_model has to be nbinom2, nbinom1, poisson, nbinom2zi or nbinom1zi for glmmTMB method!")
        if(env$method_model!="nbinom2") warning("method_model is recommended to be nbinom2 for nebula method!")
    }
    system(paste("mkdir -p",env$output))
}

scRNAseq_DE_one <- function(
    env,
    cluster_interest,
    strSrc=NULL
){
    if(!is.null(strSrc)) source(paste0(strSrc,"/pipeline_class.R"))
    message("===== read counts and meta information =====")
    counts <- readRDS(env$countRDS)
    allMeta <- readRDS(env$metaRDS)[colnames(counts),]
    allMeta$cell <- rownames(allMeta)
    message("===== ",env$method,":",cluster_interest," =====")
    strOut <- paste0(env$output,"/",env$method,"_",env$column_cluster,"/")
    system(paste("mkdir -p",strOut))
    if(!is.null(env$column_group)){
        strF <- file.path(strOut,paste0(env$grp_alt,".vs.",env$grp_ref,"_",gsub("_",".",cluster_interest),".QC.pdf"))
        sce <- BiostatsSingleCell$new(count_data = counts,
                                      meta_data = allMeta,
                                      sampleId_col = env$column_sample,
                                      cluster_col = env$column_cluster,
                                      treatment_col = env$column_group)
        sce$set_group_mode(cluster_of_interest = cluster_interest, ref_group = env$grp_ref, alt_group =env$grp_alt)

    }else{
        strF <- file.path(strOut,paste0(cluster_interest,".vs.Rest","_",gsub("_",".",env$column_cluster),".QC.pdf"))
        intrestGrp <- as.character(allMeta[,env$column_cluster])
        intrestGrp[intrestGrp!=cluster_interest] <- "others"
        meta <- cbind(allMeta,all="all",intrestGrp=intrestGrp)
        
        sce <- BiostatsSingleCell$new(count_data = counts,
                                      meta_data = meta,
                                      sampleId_col = env$column_sample,
                                      cluster_col = "all",
                                      treatment_col = "intrestGrp")
        sce$set_group_mode(cluster_of_interest = "all", ref_group = "others", alt_group = cluster_interest)
        
    }
    sce$make_QCplots(strF)
    sce$apply_filter(min.cells.per.gene = env$min.cells.per.gene, min.genes.per.cell = env$min.genes.per.cell,
                     min.perc.cells.per.gene = env$min.perc.cells.per.gene,perc_filter = env$perc_filter) # 0% expression requirement
    strF <- gsub("QC.pdf","csv",strF)
    
    if(tryCatch({
        sce_qc <- sce$apply_filter_contrasts_R6(min.cells.per.gene = env$R6_min.cells.per.gene,
                                                min.perc.cells.per.gene = env$R6_min.perc.cells.per.gene,
                                                perc.cells.filter = env$R6_perc.cells.filter,
                                                min.cells.per.gene.type = env$R6_min.cells.per.gene.type,
                                                cells.per.gene.filter = env$R6_cells.per.gene.filter,
                                                perc.filter = env$R6_perc.filter,
                                                perc.filter.type = env$R6_perc.filter.type,
                                                perc_threshold = env$R6_perc_threshold,
                                                min.ave.pseudo.bulk.cpm = env$R6_min.ave.pseudo.bulk.cpm,
                                                pseudo.bulk.cpm.filter = env$R6_pseudo.bulk.cpm.filter,
                                                min.cells.per.subj = env$R6_min.cells.per.subj)
        T
    },error=function(err){
        print(err)
        F
    })){
        covars <- env$column_covars
        system.time(de <- switch(env$method,
                                  't_test' = sce_qc$t_test_pipeline(),
                                  'u_test' = sce_qc$u_test_pipeline(),
                                  'edgeR'= sce_qc$edgeR_pipeline(covs = covars),
                                  'limma'= sce_qc$limma_pipeline(covs = covars),
                                  'DESeq2'= sce_qc$DESeq2_pipeline(covs = covars),
                                  'glmmTMB'= sce_qc$glmmTMB_pipeline(covs = covars, family = env$method_model, cores=env$core,detection_rate = FALSE),
                                  'MAST'= sce_qc$MAST_pipeline(covs = covars, detection_rate = TRUE),
                                  'limma_cell_level'= sce_qc$limma_cell_level_pipeline(covs = covars),
                                  'ancova'= sce_qc$ancova_pipeline(covs = covars),
                                  'nebula'= sce_qc$nebula_pipeline(covs = covars,method=env$method_model)))
        write.csv(de, file=strF, row.names = FALSE)
        p <- sce_qc$volcanoPlot(FDR_threshold = 0.05, FC_threshold = 2, title = env$method)
        ggsave(gsub("csv","png",strF))
    }
    
}

if(length(args)>2){
    if(length(args)>3){
        for(one in unlist(strsplit(args[4],";"))) source(one)
    }
    print(system.time(
        scRNAseq_DE_one(readRDS(paste0(args[2],"env.rds")),
                                args[3],
                                args[1])
        ))
    
}


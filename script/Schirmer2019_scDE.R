rm(list=ls())
closeAllConnections()
graphics.off()
message("Loading ...")
source("scRNAseq_DE.R",chdir=T)
strCount <- "MS_Nature_2019_Schirmer.counts.rds"
strMeta <- "MS_Nature_2019_Schirmer.meta.rds"

## run DE analysis -----------
message("Run DE analysis ...")

scRNAseq_DE(strCount,strMeta,
    output="nebula_HL",
    method="nebula",
    column_sample="sample",
    column_cluster="cell_type",
    column_group="diagnosis",
    grp_ref="Control",
    grp_alt="MS",
    method_model= "HL",
    column_covars=c("sex","mt_frac","Capbatch"),

    core=1,
    qsub=F # if qsub from sge is avaible, you can set it T
)

scRNAseq_DE(strCount,strMeta,
            output="nebula_LN",
            method="nebula",
            column_sample="sample",
            column_cluster="cell_type",
            column_group="diagnosis",
            grp_ref="MS",
            grp_alt="Control",
            method_model= "LN",
            column_covars=c("sex","mt_frac","Capbatch"),
            
            core=1,
            qsub=F # if qsub from sge is avaible, you can set it T
)

scRNAseq_DE(strCount,strMeta,
            output="glmmTMB_nbinom2",
            method="glmmTMB",
            column_sample="sample",
            column_cluster="cell_type",
            column_group="diagnosis",
            grp_ref="MS",
            grp_alt="Control",
            method_model= "nbinom2",
            column_covars=c("sex","mt_frac","Capbatch"),
            
            core=1,
            qsub=F # if qsub from sge is avaible, you can set it T
)

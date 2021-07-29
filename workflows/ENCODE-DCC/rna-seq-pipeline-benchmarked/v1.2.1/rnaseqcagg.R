#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
##USAGE: Rscript rnaseqcagg.R job_id workflow_instance_identifier workflow_identifier rep test qcfile outputfile


#args = commandArgs(trailingOnly=TRUE)
args = commandArgs(TRUE)
if (length(args)<6)
  stop("Expecting job inst work [files]")

job_id<-args[1]
workflow_inst<-args[2]
workflow_id<-args[3]
rep<-args[4]
test<-args[5]
qcfile<-args[6]
outputfile <- args[7]

read.table(qcfile,sep="\t",header=FALSE,quote="") %>% dplyr::rename("variable"=V1,"value"=V2) %>% dplyr::mutate(WorkflowId=workflow_id,
                                            WorkflowInstanceID=workflow_inst,
                                            JobRunID=job_id,
                                            rep=rep,
                                            test=test) %>% write.table(row.names=FALSE,sep="\t",file=outputfile,quote=FALSE)

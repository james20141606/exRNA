#! /usr/bin/env Rscript 
wholeanno <-read.table('wholeannotation.csv',sep=',',header=T)[seq(1,64),]
wholeanno$Class <- "Ctrl"
wholeanno[which(wholeanno$Disease_status=="HCC_before_surgery"),"Class"] <- "HCC"
suppressMessages(library("scImpute"))
print ('package loaded')
reads.qc <- readRDS("05.matrix/hcc_lulab.sequentialMap.homer.merged.clean.binnednotrim61sample.rds")
sampleLables <- c()
for(i in colnames(reads.qc)){tmp <- as.character(wholeanno[which(wholeanno$Sample_ID==i),"Class"]); sampleLables <- c(sampleLables,tmp)}
print ('start impuate')
scimpute(count_path = "05.matrix/hcc_lulab.sequentialMap.homer.merged.clean.binnednotrim61sample.rds.csv", infile = "csv", 
         outfile = "txt", out_dir = "05.matrix/imputation/binned", Kcluster = 2, ncores = 2, labels = as.vector(sampleLables), labeled = TRUE)

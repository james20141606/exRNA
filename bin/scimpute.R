#! /usr/bin/env Rscript 
suppressMessages(library("scImpute"))
reads.qc <- readRDS("05.matrix/hcc_lulab.sequentialMap.homer.merged.clean.binnednotrim61sample.rds")
sampleLables <- c()
for(i in colnames(reads.qc)){tmp <- as.character(wholeanno[which(wholeanno$Sample_ID==i),"Class"]); sampleLables <- c(sampleLables,tmp)}

scimpute(count_path = "05.matrix/hcc_lulab.sequentialMap.homer.merged.clean.binnednotrim61sample.rds.csv", infile = "csv", 
         outfile = "txt", out_dir = "05.matrix/imputation/binned", Kcluster = 2, ncores = 2, labels = as.vector(sampleLables), labeled = TRUE)
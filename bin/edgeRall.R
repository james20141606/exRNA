setwd('61samplesmatrix')
library("edgeR")
library('limma')
commonpath = '04.counts/hcc_lulab.sequentialMap.homer.'
savepath = '05.diffexp/hcc_lulab.sequentialMap.homer.'

afterind = c(1,  2,  3,  6,  8, 11, 12, 16, 18, 19, 31, 32, 33, 34, 35, 36, 38,
             40, 42) +1
beforeind=c(4,  5,  7,  9, 10, 13, 14, 15, 17, 20, 21, 22, 23, 24, 25, 26, 27,
            28, 29, 30, 37, 39, 41, 43, 44, 45, 46, 58, 59, 60, 61)+1
normalind = c(47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57)+1

dgListGroups <-c('After', 'After', 'After', 'Before', 'Before', 'After', 'Before',
                 'After', 'Before', 'Before', 'After', 'After', 'Before', 'Before',
                 'Before', 'After', 'Before', 'After', 'After', 'Before', 'Before',
                 'Before', 'Before', 'Before', 'Before', 'Before', 'Before',
                 'Before', 'Before', 'Before', 'After', 'After', 'After', 'After',
                 'After', 'After', 'Before', 'After', 'Before', 'After', 'Before',
                 'After', 'Before', 'Before', 'Before', 'Before', 'normal',
                 'normal', 'normal', 'normal', 'normal', 'normal', 'normal',
                 'normal', 'normal', 'normal', 'normal', 'Before', 'Before',
                 'Before', 'Before')
for (rnaname in c('miRNA','lncRNA','piRNA','snRNA','srpRNA','snoRNA','tRNA','Y_RNA','mRNA')){
####################################################################################
####################################################################################
raw_ <-read.table(paste(commonpath,rnaname,".mx",sep=""),header=T,sep='\t')

####################################################################################
raw <-raw_[c(afterind,beforeind)]
dgListGroups <-c(afterind,beforeind)
dgList <- DGEList(counts=raw, genes=raw_[1],group=factor(dgListGroups))
countsPerMillion <- cpm(dgList)
countCheck <- countsPerMillion >2
keep <- which(rowSums(countCheck)>=2)
dgList <- dgList[keep,]
dgList <- calcNormFactors(dgList,method="TMM")
design.mat <- model.matrix(~0 + dgList$samples$group)
colnames(design.mat) <- levels(dgList$samples$group)
d2 <- estimateGLMCommonDisp(dgList,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat)
d2 <- estimateGLMTagwiseDisp(d2,design.mat)
fit <- glmFit(d2,design.mat)
lrt <- glmLRT(fit,contrast=c(1,-1))
edgeR_result <- topTags(lrt,n=nrow(dgList))
if(!file.exists(paste(savepath,rnaname,".AftervsBefore.edgeR.tsv",sep=""))){
  write.table(edgeR_result$table,file =paste(savepath,rnaname,".AftervsBefore.edgeR.tsv",sep=""),col.names=T,row.names = F, 
              sep="\t", quote=F)}
####################################################################################
raw <-raw_[c(afterind,normalind)]
dgListGroups <-c(afterind,normalind)
dgList <- DGEList(counts=raw, genes=raw_[1],group=factor(dgListGroups))
countsPerMillion <- cpm(dgList)
countCheck <- countsPerMillion >2
keep <- which(rowSums(countCheck)>=2)
dgList <- dgList[keep,]
dgList <- calcNormFactors(dgList,method="TMM")
design.mat <- model.matrix(~0 + dgList$samples$group)
colnames(design.mat) <- levels(dgList$samples$group)
d2 <- estimateGLMCommonDisp(dgList,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat)
d2 <- estimateGLMTagwiseDisp(d2,design.mat)
fit <- glmFit(d2,design.mat)
lrt <- glmLRT(fit,contrast=c(1,-1))
edgeR_result <- topTags(lrt,n=nrow(dgList))
if(!file.exists(paste(savepath,rnaname,".AftervsNormal.edgeR.tsv",sep=""))){
  write.table(edgeR_result$table,file =paste(savepath,rnaname,".AftervsNormal.edgeR.tsv",sep=""),col.names=T,row.names = F, 
              sep="\t", quote=F)}
####################################################################################
raw <-raw_[c(beforeind,normalind)]
dgListGroups <-c(beforeind,normalind)
dgList <- DGEList(counts=raw, genes=raw_[1],group=factor(dgListGroups))
countsPerMillion <- cpm(dgList)
countCheck <- countsPerMillion >2
keep <- which(rowSums(countCheck)>=2)
dgList <- dgList[keep,]
dgList <- calcNormFactors(dgList,method="TMM")
design.mat <- model.matrix(~0 + dgList$samples$group)
colnames(design.mat) <- levels(dgList$samples$group)
d2 <- estimateGLMCommonDisp(dgList,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat)
d2 <- estimateGLMTagwiseDisp(d2,design.mat)
fit <- glmFit(d2,design.mat)
lrt <- glmLRT(fit,contrast=c(1,-1))
edgeR_result <- topTags(lrt,n=nrow(dgList))
if(!file.exists(paste(savepath,rnaname,".BeforevsNormal.edgeR.tsv",sep=""))){
  write.table(edgeR_result$table,file =paste(savepath,rnaname,".BeforevsNormal.edgeR.tsv",sep=""),col.names=T,row.names = F, 
              sep="\t", quote=F)}
}
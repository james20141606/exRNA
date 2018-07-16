setwd('61samplesmatrix')
#rnanames = 'miRNA'
rnaname = 'miRNA'
commonpath = '04.counts/hcc_lulab.sequentialMap.homer.'
savepath = '05.diffexp/hcc_lulab.sequentialMap.homer.'
mx <- read.table(paste(commonpath,rnaname,".mx",sep=""),sep="\t",header=T)
design <- read.table("05.diffexp/design.tsv",sep="\t",header=T)
filter_genes <- apply(
  mx[,2:ncol(mx)],
  1,
  function(x) length(x[x > 2]) >= 2
)
mx_filterGenes <- mx[filter_genes,]
####################################################################################
####################################################################################
library(DESeq2)
countData <- mx_filterGenes
colData <- design
dds <- DESeqDataSetFromMatrix(countData, colData, design=~Treatment, tidy=TRUE)
norm <- rlog(dds,blind=FALSE)
norm_matrix <- assay(norm)
norm_df <- data.frame(Gene=rownames(norm_matrix), norm_matrix)
if(!file.exists(paste(savepath,rnaname,".DESeq2.rlog.mx",sep=""))){
  write.table(norm_df, paste(savepath,rnaname,".DESeq2.rlog.mx",sep=""), row.names = FALSE,sep="\t")}
deg <- DESeq(dds)

res <- results(deg,contrast=c("Treatment","normal","Before"),tidy=TRUE)
merged_res <- merge(norm_df,res,by.x="Gene",by.y="row")
if(!file.exists(paste(savepath,rnaname,".NormalvsBefore.DESeq2.tsv",sep=""))){
  write.table(merged_res,file=paste(savepath,rnaname,".NormalvsBefore.DESeq2.tsv",sep=""),sep="\t",row.names=FALSE)
}

res1 <- results(deg,contrast=c("Treatment","normal","After"),tidy=TRUE)
merged_res1 <- merge(norm_df,res1,by.x="Gene",by.y="row")
if(!file.exists(paste(savepath,rnaname,".NormalvsAfter.DESeq2.tsv",sep=""))){
  write.table(merged_res1,file=paste(savepath,rnaname,".NormalvsAfter.DESeq2.tsv",sep=""),sep="\t",row.names=FALSE)
}

res2 <- results(deg,contrast=c("Treatment","After","Before"),tidy=TRUE)
merged_res2 <- merge(norm_df,res2,by.x="Gene",by.y="row")
if(!file.exists(paste(savepath,rnaname,".AftervsBefore.DESeq2.tsv",sep=""))){
  write.table(merged_res2,file=paste(savepath,rnaname,".AftervsBefore.DESeq2.tsv",sep=""),sep="\t",row.names=FALSE)
}

####################################################################################
####################################################################################
library(edgeR)
#Read Data in
countData <- mx_filterGenes[,-1]
rownames(countData) <- mx_filterGenes[,1]
design <- read.table("05.diffexp/design.tsv",sep="\t",header=T)
colData <- design
y <- DGEList(countData, samples=colData, group=colData$Treatment)
y <- calcNormFactors(y)
design <-model.matrix(~Treatment, data=colData)
y <- estimateDisp(y, design)
et <- exactTest(y)
res <- topTags(et,Inf)
tidyResult <- data.frame(Gene=rownames(res$table), res$table)
if(!file.exists(paste(savepath,rnaname,".NCvsHCC.edgeR.classic.tsv",sep=""))){
  write.table(tidyResult,file=paste(savepath,rnaname,".NCvsHCC.edgeR.classic.tsv",sep=""),sep="\t",row.names=FALSE)}
fit <- glmFit(y,design)
lrt <- glmLRT(fit,contrast = c(1,-1))
FDR <- p.adjust(lrt$table$PValue, method="BH")
padj_lrt <- cbind(lrt$table,FDR)
fit_df <- lrt$fitted.values
if(!file.exists(paste(savepath,rnaname,".homer.edgeR.TMM.mx",sep=""))){
  write.table(fit_df,file = paste(savepath,rnaname,".homer.edgeR.TMM.mx",sep=""),row.names = T, sep="\t", quote=F)}
merged_lrt <- merge(fit_df,padj_lrt,by="row.names")
colnames(merged_lrt)[1] <- "Genes"
if(!file.exists(paste(savepath,rnaname,".NCvsHCC.edgeR.tsv",sep=""))){
  write.table(merged_lrt,file =paste(savepath,rnaname,".NCvsHCC.edgeR.tsv",sep=""),row.names = F, sep="\t", quote=F)}
####################################################################################
####################################################################################

cpmMx <- read.table(paste(commonpath,rnaname,".homer.rpm.m",sep=""),sep="\t",header=T)
filter_cpm <- apply(
  mx[,2:ncol(cpmMx)],
  1,
  function(x) length(x[x > 0]) >= 2
)
mx_filterCPM <- cpmMx[filter_cpm,]
myFun <- function(x){
  x = as.numeric(x)
  v1 = x[2:4]
  v2 = x[5:10]
  out <- wilcox.test(v1,v2,exact = F) 
  out <- out$p.value
}
p_value <- apply(mx_filterCPM,1,myFun)
p_value[is.nan(p_value)] <- 1
FDR <- p.adjust(p_value,method = "BH")
mx_filterCPM$avgNC <- apply(mx_filterCPM[,2:4],1,mean)
mx_filterCPM$avgHCC <- apply(mx_filterCPM[,5:10],1,mean)
mx_filterCPM$log2fc <- log2((mx_filterCPM$avgNC+0.25)/(mx_filterCPM$avgHCC+0.25))
results <- cbind(mx_filterCPM,p_value,FDR)
if(!file.exists(paste(savepath,rnaname,".NCvsHCC.wilcox.tsv",sep=""))){
  write.table(results,file = paste(savepath,rnaname,".NCvsHCC.wilcox.tsv",sep=""),row.names = F, sep="\t", quote=F)}
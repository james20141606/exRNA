setwd('61samplesmatrix')
#rnanames = 'miRNA'
commonpath = '04.counts/hcc_lulab.sequentialMap.homer.'
savepath = '05.diffexp/hcc_lulab.sequentialMap.homer.'

for (rnaname in c('miRNA','lncRNA','piRNA','snRNA','srpRNA','snoRNA','tRNA','Y_RNA','mRNA')){

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
}
setwd('61samplesmatrix')
#rnanames = 'miRNA'
commonpath = '04.counts/hcc_lulab.sequentialMap.homer.'
savepath = '05.diffexp/hcc_lulab.sequentialMap.homer.'
afterind = c(1,  2,  3,  6,  8, 11, 12, 16, 18, 19, 31, 32, 33, 34, 35, 36, 38,
             40, 42) +1
beforeind=c(4,  5,  7,  9, 10, 13, 14, 15, 17, 20, 21, 22, 23, 24, 25, 26, 27,
            28, 29, 30, 37, 39, 41, 43, 44, 45, 46, 58, 59, 60, 61)+1
normalind = c(47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57)+1
for (rnaname in c('miRNA','lncRNA','piRNA','snRNA','srpRNA','snoRNA','tRNA','Y_RNA','mRNA')){
  
  cpmMx <- read.table(paste('04.counts/',rnaname,".homer.rpkm.mx",sep=""),sep="\t",header=T)
  filter_cpm <- apply(
    cpmMx[,2:ncol(cpmMx)],
    1,
    function(x) length(x[x > 0]) >= 2
  )
  mx_filterCPM <- cpmMx[filter_cpm,]
  ####################################################################################
  myFun <- function(x){
    x = as.numeric(x)
    v1 = x[normalind]
    v2 = x[beforeind]
    out <- wilcox.test(v1,v2,exact = F) 
    out <- out$p.value
  }
  p_value <- apply(mx_filterCPM,1,myFun)
  p_value[is.nan(p_value)] <- 1
  FDR <- p.adjust(p_value,method = "BH")
  mx_filterCPM$avgNC <- apply(mx_filterCPM[,normalind],1,mean)
  mx_filterCPM$avgHCC <- apply(mx_filterCPM[,beforeind],1,mean)
  mx_filterCPM$log2fc <- log2((mx_filterCPM$avgNC+0.25)/(mx_filterCPM$avgHCC+0.25))
  results <- cbind(mx_filterCPM,p_value,FDR)
  if(!file.exists(paste(savepath,rnaname,".NormalvsBefore.wilcox.tsv",sep=""))){
    write.table(results,file = paste(savepath,rnaname,".NormalvsBefore.wilcox.tsv",sep=""),row.names = F, sep="\t", quote=F)}
  ####################################################################################
  myFun <- function(x){
    x = as.numeric(x)
    v1 = x[normalind]
    v2 = x[afterind ]
    out <- wilcox.test(v1,v2,exact = F) 
    out <- out$p.value
  }
  p_value <- apply(mx_filterCPM,1,myFun)
  p_value[is.nan(p_value)] <- 1
  FDR <- p.adjust(p_value,method = "BH")
  mx_filterCPM$avgNC <- apply(mx_filterCPM[,normalind],1,mean)
  mx_filterCPM$avgHCC <- apply(mx_filterCPM[,afterind],1,mean)
  mx_filterCPM$log2fc <- log2((mx_filterCPM$avgNC+0.25)/(mx_filterCPM$avgHCC+0.25))
  results <- cbind(mx_filterCPM,p_value,FDR)
  if(!file.exists(paste(savepath,rnaname,".NormalvsAfter.wilcox.tsv",sep=""))){
    write.table(results,file = paste(savepath,rnaname,".NormalvsAfter.wilcox.tsv",sep=""),row.names = F, sep="\t", quote=F)}
  ####################################################################################
  myFun <- function(x){
    x = as.numeric(x)
    v1 = x[afterind ]
    v2 = x[beforeind]
    out <- wilcox.test(v1,v2,exact = F) 
    out <- out$p.value
  }
  p_value <- apply(mx_filterCPM,1,myFun)
  p_value[is.nan(p_value)] <- 1
  FDR <- p.adjust(p_value,method = "BH")
  mx_filterCPM$avgNC <- apply(mx_filterCPM[,afterind ],1,mean)
  mx_filterCPM$avgHCC <- apply(mx_filterCPM[,beforeind],1,mean)
  mx_filterCPM$log2fc <- log2((mx_filterCPM$avgNC+0.25)/(mx_filterCPM$avgHCC+0.25))
  results <- cbind(mx_filterCPM,p_value,FDR)
  if(!file.exists(paste(savepath,rnaname,".AftervsBefore.wilcox.tsv",sep=""))){
    write.table(results,file = paste(savepath,rnaname,".AftervsBefore.wilcox.tsv",sep=""),row.names = F, sep="\t", quote=F)}
  
}
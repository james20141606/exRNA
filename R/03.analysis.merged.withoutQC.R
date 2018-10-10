library(SingleCellExperiment)
library(scater)
library(scRNA.seq.funcs)
options(stringsAsFactors = FALSE)

setwd("./hcc_lulab/")


raw_mx <- read.table("hcc_lulab.sequentialMap.homer.merged.mx", sep = "\t")
anno <- read.table("hcc_lulab.sample.info.txt", sep = "\t", header=T)
anno$Class <- "Ctrl"
anno[which(anno$Disease_status=="HCC_before_surgery"),"Class"] <- "HCC"

mx <- raw_mx[, (names(raw_mx) %in% anno$Sample_ID)]

##---------------------------------------
## construct singleCellExperiment object
reads <- SingleCellExperiment(
    assays = list(counts = as.matrix(mx)),
    colData = anno)

keep_feature <- rowSums(counts(reads) > 0) > 0
reads <- reads[keep_feature, ]

# stableRNA <- isSpike(reads, "stableRNA") <- rownames(reads) %in% c("mature_miRNA:hsa-miR-99a-5p", "mature_miRNA:hsa-miR-30a-5p", "mature_miRNA:hsa-miR-221-3p")

# reads <- calculateQCMetrics(
#    reads,
#    feature_controls = list(
#        stableRNA = isSpike(reads, "stableRNA")
#    )
# )

reads <- calculateQCMetrics(reads)

##----------------------------
## filter samples and genes

hist(reads$total_counts,breaks = 100)
abline(v=10000000, col="red")
filter_by_total_counts <- (reads$total_counts > 10000000)
table(filter_by_total_counts)

hist(reads$total_features,breaks = 100)
abline(v=15000, col="red")
filter_by_expr_features <- (reads$total_features > 15000)
table(filter_by_expr_features)


# sample filtering
reads$use <- (
    # sufficient features (genes)
    filter_by_expr_features &
    # sufficient molecules counted
    filter_by_total_counts
)
table(reads$use)


## gene QC
# filter genes with too low expression
plotQC(reads, type = "highest-expression")
# top 20 
top20 <- c("Y_RNA_ENSG00000201778.1", "miRNA_ENSG00000284565.1", "miRNA_ENSG00000283364.1", "lncRNA_ENSG00000264066.6", "miRNA_ENSG00000284190.1", "miRNA_ENSG00000199179.3", "miRNA_ENSG00000199075.1", "miRNA_ENSG00000207789.1", "miRNA_ENSG00000284440.1", "miRNA_ENSG00000207778.3", "lncRNA_ENSG00000267391.4", "miRNA_ENSG00000199150.3", "miRNA_ENSG00000199161.1", "miRNA_ENSG00000199085.3", "miRNA_ENSG00000274705.2", "miRNA_ENSG00000283450.1", "miRNA_ENSG00000284520.1", "lncRNA_ENSG00000234741.7", "miRNA_ENSG00000199153.1", "miRNA_ENSG00000283733.1")

# number of cells with non-zero expression
plotQC(reads, type = "exprs-freq-vs-mean")
plotScater(reads, block1 = "Class",
     colour_by = "Batch", nfeatures = 500, exprs_values = "counts")

filter_genes <- apply(counts(reads[, colData(reads)$use]), 1, function(x) length(x[x >= 2]) >= 10)
table(filter_genes)
rowData(reads)$use <- filter_genes
reducedDim(reads) <- NULL
dim(reads[rowData(reads)$use, colData(reads)$use])


assay(reads, "logcounts_raw") <- log2(counts(reads) + 1)
reads.qc <- reads[rowData(reads)$use, colData(reads)$use]
endog_genes <- !rowData(reads.qc)$is_feature_control
# save the data
saveRDS(reads.qc, file = "hcc_lulab.sequentialMap.homer.merged.clean.rds")

##-------------------------
## normalization
# log2 raw counts
plotExpression(reads.qc, top20,
    colour_by = "Class", shape_by = "Disease_status",
    x = "Class", exprs_values="logcounts_raw")
plotPCA(
    reads.qc[endog_genes, ],
    exprs_values = "logcounts_raw",
    colour_by = "Class",
    shape_by = "Disease_status",
    size_by = "total_features"
)

plotTSNE(
    reads.qc[endog_genes, ],
    exprs_values = "logcounts_raw",
    perplexity = 8,
    colour_by = "Class",
    sharp_by = "Disease_status",
    size_by = "total_features",
    rand_seed = 123456,
    ntop = 300
)


# CPM
logcounts(reads.qc) <- log2(calculateCPM(reads.qc, use.size.factors = FALSE) + 1)

pdf("reads.qc.pca.pdf")
for(n in assayNames(reads.qc)) {
    print(
        plotPCA(
            reads.qc[endog_genes, ],
            colour_by = "Class",
            size_by = "total_features",
            shape_by = "Disease_status",
            exprs_values = n
        ) +
        ggtitle(n)
    )
}
dev.off()


# Effectiveness 1
pdf("reads.qc.pca_Batch.pdf")
for(n in assayNames(reads.qc)) {
    print(
        plotPCA(
            reads.qc[endog_genes, ],
            colour_by = "Batch",
            size_by = "total_features",
            shape_by = "Class",
            exprs_values = n
        ) +
        ggtitle(n)
    )
}
dev.off()

# Effectiveness 2
pdf("reads.qc.plotRLE.pdf")
res <- list()
for(n in assayNames(reads.qc)) {
    res[[n]] <- suppressWarnings(calc_cell_RLE(assay(reads.qc, n)))
}
par(mar=c(6,4,1,1))
boxplot(res, las=2)
dev.off()

# Effectiveness 3
pdf("reads.qc.plotQC.pdf")
for(n in assayNames(reads.qc)) {
    print(
        plotQC(
            reads.qc[endog_genes, ],
            type = "expl",
            exprs_values = n,
            variables = c(
                "total_features",
                "total_counts",
                "Batch",
                "Class",
                "Disease_status"
            )
        ) +
        ggtitle(n)
    )
}
dev.off()

saveRDS(reads.qc, file = "hcc_lulab.sequentialMap.homer.merged.withoutQC.rds")

##---------------------
## Clustering
library(SingleCellExperiment)
library(scater)
library(pcaMethods)
library(pcaReduce)
library(SC3)
library(pheatmap)
library(mclust)
set.seed(1234567)
setwd("./hcc_lulab/")

reads.qc <- readRDS("hcc_lulab.sequentialMap.homer.merged.clean.impute.rds")

# using SC3
reads.qc <- sc3_estimate_k(reads.qc)
metadata(reads.qc)$sc3$k_estimation
plotPCA(reads.qc, exprs_values = "ruvs20", colour_by = "Class", size_by = "total_features", shape_by = "Disease_status")
# run sc3
reads.qc <- sc3(reads.qc, ks = 2, biology = TRUE)
# consensus matrix
sc3_plot_consensus(reads.qc, k = 2, show_pdata = "Class")
# Silhouette plot
sc3_plot_silhouette(reads.qc, k = 2)
# heatmap
sc3_plot_expression(reads.qc, k = 2, show_pdata = "Class")
# identify marker genes
sc3_plot_markers(reads.qc, k = 2, show_pdata = "Class")



##-------------------------------
## diff exp 
library(scRNA.seq.funcs)
library(edgeR)
library(monocle)
library(MAST)
library(ROCR)
set.seed(1)

cts <- counts(reads.qc)
group <- as.data.frame(colData(reads.qc)[,1:4])

group <- read.table("hcc_lulab.sample.info.txt", sep = "\t", header=T)
group$Class <- "Ctrl"
group[which(group$Disease_status=="HCC_before_surgery"),"Class"] <- "HCC"
group <- group[group$Sample_ID %in% colnames(reads.qc),]

# using DESeq2

library(DESeq2)
#Read Data in
colData <- group
countData <- data.frame(geneID = row.names(cts), ceiling(cts))
rownames(countData) <- NULL

# generate Dataset
dds <- DESeqDataSetFromMatrix(countData, colData, design=~Class, tidy=TRUE)
deg <- DESeq(dds)
res <- results(deg,tidy=TRUE)
res_ordered <- res[order(res$padj),]

deseq2.diffexp.fdr05 <- res_ordered[which(res_ordered$padj<=0.05),]
deseq2.diffexp.fdr05.list <- deseq2.diffexp.fdr05$row
deseq2.diffexp.fdr001 <- res_ordered[which(res_ordered$padj<=0.001),]
deseq2.diffexp.fdr001.list <- deseq2.diffexp.fdr001$row

deseq2.diffexp.top.200 <- res_ordered[1:200,]
deseq2.diffexp.top.200.list <- deseq2.diffexp.top.200$row

write.table(res_ordered,file="hcc_lulab.sequentialMap.homer.merged.CtrlvsHCC.DESeq2.tsv",sep="\t", quote = F, row.names=FALSE)


# using edgeR
library(edgeR)
countData <- counts(reads.qc)
design <- read.table("hcc_lulab.sample.info.txt", sep = "\t", header=T)
design$Class <- "Ctrl"
design[which(design$Disease_status=="HCC_before_surgery"),"Class"] <- "HCC"
design <- design[design$Sample_ID %in% colnames(reads.qc),]
colData <- design
design <-model.matrix(~Class, data=colData)

y <- DGEList(countData, samples=colData, group=colData$Class)
y <- calcNormFactors(y)
y <- estimateDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
top <- topTags(lrt, n=nrow(countData))$table
empirical <- rownames(countData)[which(!(rownames(countData) %in% rownames(top)[1:5000]))]


top_ordered <- top[order(top$FDR),]

edgeR.diffexp.fdr05 <- top_ordered[which(top_ordered$FDR<=0.05),]
edgeR.diffexp.fdr05.list <- rownames(edgeR.diffexp.fdr05)
edgeR.diffexp.fdr001 <- top_ordered[which(top_ordered$FDR<=0.001),]
edgeR.diffexp.fdr001.list <- rownames(edgeR.diffexp.fdr001)

edgeR.diffexp.top.200 <- top_ordered[1:200,]
edgeR.diffexp.top.200.list <- rownames(edgeR.diffexp.top.200)


write.table(top_ordered,file="hcc_lulab.sequentialMap.homer.merged.CtrlvsHCC.edgeR.tsv",sep="\t",quote = F, row.names=T)


##--------------------------
## feature selection
library(scRNA.seq.funcs)
library(matrixStats)
library(scran)
library(M3Drop)
library(RColorBrewer)
library(SingleCellExperiment)
set.seed(1)
setwd("./hcc_lulab")

reads.qc <- readRDS("hcc_lulab.sequentialMap.homer.merged.withoutQC.rds")
anno <- read.table("hcc_lulab.sample.info.txt", sep = "\t", header=T)
anno$Class <- "Ctrl"
anno[which(anno$Disease_status=="HCC_before_surgery"),"Class"] <- "HCC"

sampleLables <- c()
for(i in colnames(reads.qc)){tmp <- as.character(anno[which(anno$Sample_ID==i),"Class"]); sampleLables <- c(sampleLables,tmp)}

# highly variable genes without spike-ins
Brennecke_HVG <- BrenneckeGetVariableGenes(
    counts(reads.qc),
    fdr = 0.05,
    minBiolDisp = 0.5
)
HVG_genes <- Brennecke_HVG$Gene

reads.qc.HVG_genes.top.200 <- Brennecke_HVG[1:200, "Gene"]
reads.qc.HVG_genes.fdr05 <- Brennecke_HVG[which(Brennecke_HVG$q.value<=0.05),"Gene"]
reads.qc.HVG_genes.fdr001 <- Brennecke_HVG[which(Brennecke_HVG$q.value<=0.001),"Gene"]


norm.counts <- assay(reads.qc, i="logcounts")
reads.qc.fit <- trendVar(norm.counts)
reads.qc.decomp <- decomposeVar(norm.counts, reads.qc.fit)
reads.qc.decomp.ordered <- reads.qc.decomp[order(reads.qc.decomp$FDR),]

reads.qc.top_hvgs.200 <- rownames(reads.qc.decomp.ordered[1:200,])
reads.qc.top_hvgs.fdr05 <- rownames(reads.qc.decomp.ordered[which(reads.qc.decomp.ordered$FDR<=0.05),])
reads.qc.top_hvgs.fdr001 <- rownames(reads.qc.decomp.ordered[which(reads.qc.decomp.ordered$FDR<=0.001),])


# highly Dropout genes 
M3Drop_genes <- M3DropFeatureSelection(
    counts(reads.qc),
    mt_method = "fdr",
    mt_threshold = 0.05
)
M3Drop_genes_ordered <- M3Drop_genes[order(M3Drop_genes$q.value),]
reads.qc.top_m3drop.top.200 <- rownames(M3Drop_genes_ordered[1:200,])
reads.qc.top_m3drop.fdr05 <- rownames(M3Drop_genes_ordered[which(M3Drop_genes_ordered$q.value<=0.05),])
reads.qc.top_m3drop.fdr001 <- rownames(M3Drop_genes_ordered[which(M3Drop_genes_ordered$q.value<=0.001),])


## correlation genes
norm.counts <- assay(reads.qc, i="logcounts")
cor_mat <- cor(t(norm.counts), method = "spearman") # Gene-gene correlations
diag(cor_mat) <- rep(0, times = nrow(norm.counts))
score <- apply(cor_mat, 1, function(x) {max(abs(x))}) #Correlation of highest magnitude
names(score) <- rownames(norm.counts);
score <- score[order(-score)]

Cor_genes.top.200 <- names(score[1:200])
Cor_genes.top.500 <- names(score[1:500])
Cor_genes.top.1000 <- names(score[1:1000])


# PCA loading
pca <- prcomp(assay(reads.qc, i="logcounts"))

# plot projection
sample_colors <- brewer.pal(max(3,length(unique(sampleLables))), "Set2")
plot(
    pca$rotation[,1], 
    pca$rotation[,2], 
    pch = 16, 
    col = sample_colors[as.factor(sampleLables)]
) 

score <- rowSums(abs(pca$x[,c(1,2)])) 
names(score) <- rownames(assay(reads.qc, i="counts"))
score <- score[order(-score)]

PCA_genes.top.200 <- names(score[1:200])
PCA_genes.top.500 <- names(score[1:500])
PCA_genes.top.1000 <- names(score[1:1000])

## visulization
library(pheatmap)
norm.counts <- assay(reads.qc, i="logcounts")
anno <- read.table("hcc_lulab.sample.info.txt", sep = "\t", header=T)
anno$Class <- "Ctrl"
anno[which(anno$Disease_status=="HCC_before_surgery"),"Class"] <- "HCC"
group <- anno[match(colnames(norm.counts), anno$Sample_ID),-1]
rownames(group) <- anno[match(colnames(norm.counts), anno$Sample_ID),1]


# edgeR
pheatmap(logcounts(reads.qc)[rownames(logcounts(reads.qc)) %in% edgeR.diffexp.top.200.list,], scale = "row", cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=group)
pheatmap(norm.ruvs20[rownames(norm.counts) %in% edgeR.diffexp.top.200.list,], scale = "row", cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=group)
# deseq2
pheatmap(logcounts(reads.qc)[rownames(logcounts(reads.qc)) %in% deseq2.diffexp.top.200.list,], scale = "row", cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=group)
pheatmap(norm.ruvs20[rownames(norm.counts) %in% deseq2.diffexp.top.200.list,], scale = "row", cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=group)

M3DropExpressionHeatmap(HVG_genes, counts(reads.qc), cell_labels = sampleLables)


# 7 miRNA
# miR-122, miR-192, miR-21, miR-223, miR-26a, miR-27a and miR-801
miRNA7 <- c("miRNA_ENSG00000284440.1", "miRNA_ENSG00000283926.1", "miRNA_ENSG00000284190.1", "miRNA_ENSG00000284567.1", "miRNA_ENSG00000199075.1", "miRNA_ENSG00000207808.1")
pheatmap(logcounts(reads.qc)[rownames(logcounts(reads.qc)) %in% miRNA7,], scale = "row", cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=group)
M3DropExpressionHeatmap(miRNA7, counts(reads.qc), cell_labels = sampleLables)


##------------------
## plot heatmap

genelist <- list(miRNA7, edgeR.diffexp.top.200.list, edgeR.diffexp.fdr05.list, edgeR.diffexp.fdr001.list, deseq2.diffexp.top.200.list, deseq2.diffexp.fdr05.list, deseq2.diffexp.fdr001.list, reads.qc.HVG_genes.top.200, reads.qc.HVG_genes.fdr05, reads.qc.HVG_genes.fdr001, reads.qc.top_hvgs.200, reads.qc.top_hvgs.fdr05, reads.qc.top_hvgs.fdr001, reads.qc.top_m3drop.top.200, reads.qc.top_m3drop.fdr05, reads.qc.top_m3drop.fdr001, Cor_genes.top.200, Cor_genes.top.500, Cor_genes.top.1000, PCA_genes.top.200, PCA_genes.top.500, PCA_genes.top.1000)
names(genelist) <- c("miRNA7", "edgeR.diffexp.top.200.list", "edgeR.diffexp.fdr05.list", "edgeR.diffexp.fdr001.list", "deseq2.diffexp.top.200.list", "deseq2.diffexp.fdr05.list", "deseq2.diffexp.fdr001.list", "reads.qc.HVG_genes.top.200", "reads.qc.HVG_genes.fdr05", "reads.qc.HVG_genes.fdr001", "reads.qc.top_hvgs.200", "reads.qc.top_hvgs.fdr05", "reads.qc.top_hvgs.fdr001", "reads.qc.top_m3drop.top.200", "reads.qc.top_m3drop.fdr05", "reads.qc.top_m3drop.fdr001", "Cor_genes.top.200", "Cor_genes.top.500", "Cor_genes.top.1000", "PCA_genes.top.200", "PCA_genes.top.500", "PCA_genes.top.1000")

save(genelist,file="genelist.withoutQC.Rdata")

pdf("reads.qc.selectgenes.heatmap.pdf")
for (i in 1:length(genelist)){
title <- names(genelist[i])
pheatmap(logcounts(reads.qc)[rownames(logcounts(reads.qc)) %in% genelist[[i]],], scale = "row", cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=group, main=title)
M3DropExpressionHeatmap(genelist[[i]], counts(reads.qc), cell_labels = sampleLables)
}
dev.off()


##---------------------
## machine learning

library(SingleCellExperiment)
library(scater)
options(stringsAsFactors = FALSE)

reads.qc <- readRDS("hcc_lulab.sequentialMap.homer.merged.withoutQC.rds")
load("genelist.withoutQC.Rdata")

# logistic regression in R
library(caret)

anno <- read.table("hcc_lulab.sample.info.txt", sep = "\t", header=T)
anno$Class <- "Ctrl"
anno[which(anno$Disease_status=="HCC_before_surgery"),"Class"] <- "HCC"
sampleLables <- c()
for(i in colnames(reads.qc)){tmp <- as.character(anno[which(anno$Sample_ID==i),"Class"]); sampleLables <- c(sampleLables,tmp)}

# train the model 
inputF <- t(assay(reads.qc, i="logcounts"))
inputF <- apply(inputF, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
inputF <- data.frame(sampleID = rownames(inputF), lables = sampleLables, inputF)
rownames(inputF) <- NULL

library(ROCR)

outlist <- list() 
pdf("reads.qc.selectgenes.glm.roc.pdf")
for (i in 1:length(genelist)){
title <- names(genelist[i])

model<- train(inputF[,colnames(inputF) %in% genelist[[i]]], inputF[,2], method="glm", family=binomial())

pred_raw = predict(model, newdata=inputF[,colnames(inputF) %in% genelist[[i]]], type="raw")
pred_prob = predict(model, newdata=inputF[,colnames(inputF) %in% genelist[[i]]], type="prob")

stat <- confusionMatrix(data=pred_raw, as.factor(inputF$lables))
acc <- paste("acc:", stat$overall[1], sep=" ")
sen <- paste("acc:", stat$byClass[1], sep=" ")
spec <- paste("spec:", stat$byClass[2], sep=" ")

pred_raw.1 <- as.vector(pred_raw)
pred_raw.1[pred_raw.1 == "HCC"] <- 1
pred_raw.1[pred_raw.1 == "Ctrl"] <- 0
pred_raw.1 <- as.numeric(pred_raw.1)

pred_prob <- pred_prob$HCC

goldlable <- inputF$lables
goldlable[goldlable == "HCC"] <- 1
goldlable[goldlable == "Ctrl"] <- 0
goldlable <- as.numeric(goldlable)

pred <- prediction(pred_prob, goldlable)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
auc <- performance(pred, measure = "auc")
auc_val <- auc@y.values[[1]]
addtext <- paste("AUC:", auc_val, sep=" ")

out <- c(title,stat$overall[1],stat$byClass[1], stat$byClass[2], auc_val)
outlist[[i]] <- out

plot(perf, col=rainbow(10), main=title)
mtext(addtext, side=1, line = 0)
mtext(acc, side=1, line = 1)
mtext(sen, side=1, line = 2)
mtext(spec, side=1, line = 3)
}
dev.off()
outmx <- do.call(rbind, outlist)
outmx <- as.data.frame(outmx)
colnames(outmx) <- c("genesets","Accuracy", "Sensitivity", "Specificity", "AUC")
write.table(outmx, "./hcc_lulab.sequentialMap.homer.merged.withoutQC.glm.stat.txt", sep="\t",quote=F, row.names=F, col.names=T)


# output genelist


 



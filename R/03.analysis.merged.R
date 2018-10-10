library(SingleCellExperiment)
library(scater)
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
abline(v=19000000, col="red")
filter_by_total_counts <- (reads$total_counts > 19000000)
table(filter_by_total_counts)

hist(reads$total_features,breaks = 100)
abline(v=15000, col="red")
filter_by_expr_features <- (reads$total_features > 15000)
table(filter_by_expr_features)

# automitic filter
reads <- plotPCA(
    reads,
    size_by = "total_features",
    shape_by = "Disease_status",
    pca_data_input = "pdata",
    detect_outliers = TRUE,
    return_SCE = TRUE
)

filter_by_outlier <- !reads$outlier

# sample filtering
reads$use <- (
    # sufficient features (genes)
    filter_by_expr_features &
    # sufficient molecules counted
    filter_by_total_counts &
    # automatic filtering
    filter_by_outlier
)
table(reads$use)


## gene QC
# filter genes with too low expression
plotQC(reads, type = "highest-expression")
# top 20 
top20 <- c("Y_RNA_ENSG00000201778.1", "miRNA_ENSG00000284565.1", "miRNA_ENSG00000283364.1", "lncRNA_ENSG00000264066.6", "miRNA_ENSG00000284190.1", "miRNA_ENSG00000199179.3", "miRNA_ENSG00000199075.1", "miRNA_ENSG00000207789.1", "miRNA_ENSG00000284440.1", "miRNA_ENSG00000207778.3", "lncRNA_ENSG00000267391.4", "miRNA_ENSG00000199150.3", "miRNA_ENSG00000199161.1", "miRNA_ENSG00000199085.3", "miRNA_ENSG00000274705.2", "miRNA_ENSG00000283450.1", "miRNA_ENSG00000284520.1", "lncRNA_ENSG00000234741.7", "miRNA_ENSG00000199153.1", "miRNA_ENSG00000283733.1")

top3 <- c("Y_RNA_ENSG00000201778.1", "miRNA_ENSG00000284565.1", "miRNA_ENSG00000284440.1")
plotExpression(reads.qc, top3,
    exprs_values = "logcounts",
    colour_by = "Class", shape_by = "Disease_status")



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

##---------------------
## Imputation
library("scImpute")
library(scran)
write.csv(counts(reads.qc), "hcc_lulab.sequentialMap.homer.merged.clean.rds.csv")

reads.qc <- readRDS("hcc_lulab.sequentialMap.homer.merged.clean.rds")
sampleLables <- c()
for(i in colnames(reads.qc)){tmp <- as.character(anno[which(anno$Sample_ID==i),"Class"]); sampleLables <- c(sampleLables,tmp)}

scimpute(count_path = "hcc_lulab.sequentialMap.homer.merged.clean.rds.csv", infile = "csv", outfile = "txt", out_dir = "./", Kcluster = 2, ncores = 2, labels = as.vector(sampleLables), labeled = TRUE)

res.qc <- read.table("scimpute_count.txt")
reads.qc.impute <- SingleCellExperiment(assays = list(counts = as.matrix(res.qc)), colData = colData(reads.qc))
reads.qc.impute <- calculateQCMetrics(reads.qc.impute)
assay(reads.qc.impute, "logcounts_raw") <- log2(counts(reads.qc.impute) + 1)


# without pre-cluster
scimpute(count_path = "hcc_lulab.sequentialMap.homer.merged.clean.rds.csv", infile = "csv", outfile = "txt", out_dir = "./withoutCluster.", ncores = 2, Kcluster = 2)
res.qc <- read.table("./withoutCluster.scimpute_count.txt")
reads.qc.impute <- SingleCellExperiment(assays = list(counts = as.matrix(res.qc)), colData = colData(reads.qc))
reads.qc.impute <- calculateQCMetrics(reads.qc.impute)
assay(reads.qc.impute, "logcounts_raw") <- log2(counts(reads.qc.impute) + 1)


##-------------------------
## normalization
# log2 raw counts
plotExpression(reads.qc.impute, top20,
    colour_by = "Class", shape_by = "Disease_status",
    x = "Class", exprs_values="logcounts_raw")
plotPCA(
    reads.qc.impute[endog_genes, ],
    exprs_values = "logcounts_raw",
    colour_by = "Class",
    shape_by = "Disease_status",
    size_by = "total_features"
)

plotTSNE(
    reads.qc.impute[endog_genes, ],
    exprs_values = "logcounts_raw",
    perplexity = 8,
    colour_by = "Class",
    sharp_by = "Disease_status",
    size_by = "total_features",
    rand_seed = 123456,
    ntop = 300
)


# CPM
logcounts(reads.qc.impute) <- log2(calculateCPM(reads.qc.impute, use.size.factors = FALSE) + 1)

# scran (CPM)
library(scran)
# define cluster for each sample
sampleLables <- c()
for(i in colnames(reads.qc.impute)){tmp <- as.character(anno[which(anno$Sample_ID==i),"Class"]); sampleLables <- c(sampleLables,tmp)}
sampleLables <- replace(sampleLables, which(sampleLables=="HCC"),1)
sampleLables <- replace(sampleLables, which(sampleLables=="Ctrl"),2)
sampleLables <- as.numeric(sampleLables)

# sampleLables <- quickCluster(reads.qc.impute, min.size = 10)
reads.qc.impute <- computeSumFactors(reads.qc.impute, sizes = 5, clusters = sampleLables)
reads.qc.impute <- normalize(reads.qc.impute)


plotPCA(
    reads.qc.impute,
    exprs_values = "logcounts",
    colour_by = "Class",
    shape_by = "Disease_status",
    size_by = "total_features",
    ntop = 300
)

plotTSNE(
    reads.qc.impute,
    exprs_values = "logcounts",
    perplexity = 8,
    colour_by = "Class",
    shape_by = "Disease_status",
    size_by = "total_features",
    rand_seed = 123456,
    ntop = 300
)

plotRLE(
    reads.qc.impute,
    exprs_mats = list(Raw = "counts", scran = "logcounts"),
    exprs_logged = c(TRUE, TRUE),
    colour_by = "Disease_status",
    legend = "row"
)

##------------------------------
## Dealing with confounders
## check confounders
plotPCA(
    reads.qc.impute[endog_genes, ],
    exprs_values = "logcounts_raw",
    colour_by = "Batch",
    size_by = "total_features"
)

plotQC(
    reads.qc.impute[endog_genes, ],
    type = "find-pcs",
    exprs_values = "logcounts_raw",
    variable = "total_features"
)

plotQC(
    reads.qc.impute[endog_genes, ],
    type = "expl",
    exprs_values = "logcounts_raw",
    variables = c(
        "total_features",
        "total_counts",
        "Batch",
        "Disease_status",
        "Class"
    )
)

## remove batch effect using enfogenous genes
library(RUVSeq)
library(sva)
library(scRNA.seq.funcs)

# using RUVseq (RUVs)
# RUVs uses centered (technical) replicate/negative control samples for which the covariates of interest are constant
scIdx <- matrix(-1, ncol = max(table(reads.qc.impute$Class)), nrow = 2)
tmp <- which(reads.qc.impute$Class == "HCC")
scIdx[1, 1:length(tmp)] <- tmp
tmp <- which(reads.qc.impute$Class == "Ctrl")
scIdx[2, 1:length(tmp)] <- tmp
cIdx <- rownames(reads.qc.impute)
ruvs <- RUVs(logcounts(reads.qc.impute), cIdx, k = 1, scIdx = scIdx, isLog = TRUE)
assay(reads.qc.impute, "ruvs1") <- ruvs$normalizedCounts
ruvs <- RUVs(logcounts(reads.qc.impute), cIdx, k = 5, scIdx = scIdx, isLog = TRUE)
assay(reads.qc.impute, "ruvs5") <- ruvs$normalizedCounts
ruvs <- RUVs(logcounts(reads.qc.impute), cIdx, k = 10, scIdx = scIdx, isLog = TRUE)
assay(reads.qc.impute, "ruvs10") <- ruvs$normalizedCounts
ruvs <- RUVs(logcounts(reads.qc.impute), cIdx, k = 20, scIdx = scIdx, isLog = TRUE)
assay(reads.qc.impute, "ruvs20") <- ruvs$normalizedCounts


library(edgeR)
countData <- counts(reads.qc.impute)
design <- read.table("hcc_lulab.sample.info.txt", sep = "\t", header=T)
design$Class <- "Ctrl"
design[which(design$Disease_status=="HCC_before_surgery"),"Class"] <- "HCC"
design <- design[design$Sample_ID %in% colnames(reads.qc.impute),]
colData <- design
design <-model.matrix(~Class, data=colData)

y <- DGEList(countData, samples=colData, group=colData$Class)
y <- calcNormFactors(y)
y <- estimateDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
top <- topTags(lrt, n=nrow(countData))$table
empirical <- rownames(countData)[which(!(rownames(countData) %in% rownames(top)[1:5000]))]


ruvg <- RUVg(logcounts(reads.qc.impute), empirical, k = 10, isLog = TRUE)
assay(reads.qc.impute, "ruvg10") <- ruvg$normalizedCounts
ruvg <- RUVg(logcounts(reads.qc.impute), empirical, k = 20, isLog = TRUE)
assay(reads.qc.impute, "ruvg20") <- ruvg$normalizedCounts


# using Combat
combat_data <- logcounts(reads.qc.impute)
mod_data <- as.data.frame(t(combat_data))
# Basic batch removal
mod0 = model.matrix(~ 1, data = mod_data) 
# Preserve biological variability
mod1 = model.matrix(~ reads.qc.impute$Class, data = mod_data) 
# adjust for total genes detected
mod2 = model.matrix(~ reads.qc.impute$total_features, data = mod_data)
assay(reads.qc.impute, "combat") <- ComBat(
    dat = t(mod_data), 
    batch = factor(reads.qc.impute$Batch), 
    mod = mod1,
    par.prior = TRUE,
    prior.plots = FALSE
)


# using GLM
glm_fun <- function(g, batch, class) {
  model <- glm(g ~ batch + class)
  model$coef[1] <- 0 # replace intercept with 0 to preserve reference batch.
  return(model$coef)
}
effects <- apply(
    logcounts(reads.qc.impute), 
    1, 
    glm_fun, 
    batch = reads.qc.impute$Batch, 
    class = reads.qc.impute$Class
)
corrected <- logcounts(reads.qc.impute) - t(effects[as.numeric(factor(reads.qc.impute$Batch)), ])
assay(reads.qc.impute, "glm") <- corrected


# check PCA
pdf("reads.qc.impute.remove_confounders.pca.pdf")
for(n in assayNames(reads.qc.impute)) {
    print(
        plotPCA(
            reads.qc.impute[endog_genes, ],
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
pdf("reads.qc.impute.remove_confounders.pca_Batch.pdf")
for(n in assayNames(reads.qc.impute)) {
    print(
        plotPCA(
            reads.qc.impute[endog_genes, ],
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
pdf("reads.qc.impute.remove_confounders.plotRLE.pdf")
res <- list()
for(n in assayNames(reads.qc.impute)) {
    res[[n]] <- suppressWarnings(calc_cell_RLE(assay(reads.qc.impute, n)))
}
par(mar=c(6,4,1,1))
boxplot(res, las=2)
dev.off()

# Effectiveness 3
pdf("reads.qc.impute.remove_confounders.plotQC.pdf")
for(n in assayNames(reads.qc.impute)) {
    print(
        plotQC(
            reads.qc.impute[endog_genes, ],
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

saveRDS(reads.qc.impute, file = "hcc_lulab.sequentialMap.homer.merged.clean.impute.rds")

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

reads.qc.impute <- readRDS("hcc_lulab.sequentialMap.homer.merged.clean.impute.rds")

# using SC3
reads.qc.impute <- sc3_estimate_k(reads.qc.impute)
metadata(reads.qc.impute)$sc3$k_estimation
plotPCA(reads.qc.impute, exprs_values = "ruvs20", colour_by = "Class", size_by = "total_features", shape_by = "Disease_status")
# run sc3
reads.qc.impute <- sc3(reads.qc.impute, ks = 2, biology = TRUE)
# consensus matrix
sc3_plot_consensus(reads.qc.impute, k = 2, show_pdata = "Class")
# Silhouette plot
sc3_plot_silhouette(reads.qc.impute, k = 2)
# heatmap
sc3_plot_expression(reads.qc.impute, k = 2, show_pdata = "Class")
# identify marker genes
sc3_plot_markers(reads.qc.impute, k = 2, show_pdata = "Class")



##-------------------------------
## diff exp 
library(scRNA.seq.funcs)
library(edgeR)
library(monocle)
library(MAST)
library(ROCR)
set.seed(1)

cts <- counts(reads.qc.impute)
group <- as.data.frame(colData(reads.qc.impute)[,1:4])

group <- read.table("hcc_lulab.sample.info.txt", sep = "\t", header=T)
group$Class <- "Ctrl"
group[which(group$Disease_status=="HCC_before_surgery"),"Class"] <- "HCC"
group <- group[group$Sample_ID %in% colnames(reads.qc.impute),]

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

write.table(res_ordered,file="hcc_lulab.sequentialMap.homer.merged.clean.impute.CtrlvsHCC.DESeq2.tsv",sep="\t", quote = F, row.names=FALSE)


# using edgeR
library(edgeR)
countData <- counts(reads.qc.impute)
design <- read.table("hcc_lulab.sample.info.txt", sep = "\t", header=T)
design$Class <- "Ctrl"
design[which(design$Disease_status=="HCC_before_surgery"),"Class"] <- "HCC"
design <- design[design$Sample_ID %in% colnames(reads.qc.impute),]
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


write.table(top_ordered,file="hcc_lulab.sequentialMap.homer.merged.clean.impute.CtrlvsHCC.edgeR.tsv",sep="\t",quote = F, row.names=T)


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

reads.qc.impute <- readRDS("hcc_lulab.sequentialMap.homer.merged.clean.impute.rds")
anno <- read.table("hcc_lulab.sample.info.txt", sep = "\t", header=T)
anno$Class <- "Ctrl"
anno[which(anno$Disease_status=="HCC_before_surgery"),"Class"] <- "HCC"

sampleLables <- c()
for(i in colnames(reads.qc.impute)){tmp <- as.character(anno[which(anno$Sample_ID==i),"Class"]); sampleLables <- c(sampleLables,tmp)}

# highly variable genes without spike-ins
Brennecke_HVG <- BrenneckeGetVariableGenes(
    counts(reads.qc.impute),
    fdr = 0.05,
    minBiolDisp = 0.5
)
HVG_genes <- Brennecke_HVG$Gene

reads.qc.impute.HVG_genes.top.200 <- Brennecke_HVG[1:200, "Gene"]
reads.qc.impute.HVG_genes.fdr05 <- Brennecke_HVG[which(Brennecke_HVG$q.value<=0.05),"Gene"]
reads.qc.impute.HVG_genes.fdr001 <- Brennecke_HVG[which(Brennecke_HVG$q.value<=0.001),"Gene"]


norm.ruvs20 <- assay(reads.qc.impute, i="ruvs20")
reads.qc.impute.fit <- trendVar(norm.ruvs20)
reads.qc.impute.decomp <- decomposeVar(norm.ruvs20, reads.qc.impute.fit)
reads.qc.impute.decomp.ordered <- reads.qc.impute.decomp[order(reads.qc.impute.decomp$FDR),]

reads.qc.impute.top_hvgs.200 <- rownames(reads.qc.impute.decomp.ordered[1:200,])
reads.qc.impute.top_hvgs.fdr05 <- rownames(reads.qc.impute.decomp.ordered[which(reads.qc.impute.decomp.ordered$FDR<=0.05),])
reads.qc.impute.top_hvgs.fdr001 <- rownames(reads.qc.impute.decomp.ordered[which(reads.qc.impute.decomp.ordered$FDR<=0.001),])


# highly Dropout genes 
M3Drop_genes <- M3DropFeatureSelection(
    counts(reads.qc.impute),
    mt_method = "fdr",
    mt_threshold = 0.05
)
M3Drop_genes_ordered <- M3Drop_genes[order(M3Drop_genes$q.value),]
reads.qc.impute.top_m3drop.top.200 <- rownames(M3Drop_genes_ordered[1:200,])
reads.qc.impute.top_m3drop.fdr05 <- rownames(M3Drop_genes_ordered[which(M3Drop_genes_ordered$q.value<=0.05),])
reads.qc.impute.top_m3drop.fdr001 <- rownames(M3Drop_genes_ordered[which(M3Drop_genes_ordered$q.value<=0.001),])


# comparing methods
M3DropExpressionHeatmap(HVG_genes, reads.qc.impute.list.expr_matrix, cell_labels = reads.qc.impute.list.lables)
M3DropExpressionHeatmap(reads.qc.impute.top_hvgs.100, reads.qc.impute.list.expr_matrix, cell_labels = reads.qc.impute.list.lables)


## correlation genes
norm.ruvs20 <- assay(reads.qc.impute, i="ruvs20")
cor_mat <- cor(t(norm.ruvs20), method = "spearman") # Gene-gene correlations
diag(cor_mat) <- rep(0, times = nrow(norm.ruvs20))
score <- apply(cor_mat, 1, function(x) {max(abs(x))}) #Correlation of highest magnitude
names(score) <- rownames(norm.ruvs20);
score <- score[order(-score)]

Cor_genes.top.200 <- names(score[1:200])
Cor_genes.top.500 <- names(score[1:500])
Cor_genes.top.1000 <- names(score[1:1000])


# PCA loading
pca <- prcomp(assay(reads.qc.impute, i="ruvs20"))

# plot projection
sample_colors <- brewer.pal(max(3,length(unique(sampleLables))), "Set2")
plot(
    pca$rotation[,1], 
    pca$rotation[,2], 
    pch = 16, 
    col = sample_colors[as.factor(sampleLables)]
) 

score <- rowSums(abs(pca$x[,c(1,2)])) 
names(score) <- rownames(assay(reads.qc.impute, i="ruvs20"))
score <- score[order(-score)]

PCA_genes.top.200 <- names(score[1:200])
PCA_genes.top.500 <- names(score[1:500])
PCA_genes.top.1000 <- names(score[1:1000])

## visulization
library(pheatmap)
norm.ruvs20 <- assay(reads.qc.impute, i="ruvs20")
anno <- read.table("hcc_lulab.sample.info.txt", sep = "\t", header=T)
anno$Class <- "Ctrl"
anno[which(anno$Disease_status=="HCC_before_surgery"),"Class"] <- "HCC"
group <- anno[match(colnames(norm.ruvs20), anno$Sample_ID),-1]
rownames(group) <- anno[match(colnames(norm.ruvs20), anno$Sample_ID),1]


# edgeR
pheatmap(logcounts(reads.qc.impute)[rownames(logcounts(reads.qc.impute)) %in% edgeR.diffexp.top.200.list,], scale = "row", cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=group)
pheatmap(norm.ruvs20[rownames(norm.ruvs20) %in% edgeR.diffexp.top.200.list,], scale = "row", cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=group)
# deseq2
pheatmap(logcounts(reads.qc.impute)[rownames(logcounts(reads.qc.impute)) %in% deseq2.diffexp.top.200.list,], scale = "row", cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=group)
pheatmap(norm.ruvs20[rownames(norm.ruvs20) %in% deseq2.diffexp.top.200.list,], scale = "row", cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=group)

M3DropExpressionHeatmap(HVG_genes, counts(reads.qc.impute), cell_labels = sampleLables)


# 7 miRNA
# miR-122, miR-192, miR-21, miR-223, miR-26a, miR-27a and miR-801
miRNA7 <- c("miRNA_ENSG00000284440.1", "miRNA_ENSG00000283926.1", "miRNA_ENSG00000284190.1", "miRNA_ENSG00000284567.1", "miRNA_ENSG00000199075.1", "miRNA_ENSG00000207808.1")
pheatmap(logcounts(reads.qc.impute)[rownames(logcounts(reads.qc.impute)) %in% miRNA7,], scale = "row", cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=group)
M3DropExpressionHeatmap(miRNA7, counts(reads.qc.impute), cell_labels = sampleLables)


##------------------
## plot heatmap

genelist <- list(miRNA7, edgeR.diffexp.top.200.list, edgeR.diffexp.fdr05.list, edgeR.diffexp.fdr001.list, deseq2.diffexp.top.200.list, deseq2.diffexp.fdr05.list, deseq2.diffexp.fdr001.list, reads.qc.impute.HVG_genes.top.200, reads.qc.impute.HVG_genes.fdr05, reads.qc.impute.HVG_genes.fdr001, reads.qc.impute.top_hvgs.200, reads.qc.impute.top_hvgs.fdr05, reads.qc.impute.top_hvgs.fdr001, reads.qc.impute.top_m3drop.top.200, reads.qc.impute.top_m3drop.fdr05, reads.qc.impute.top_m3drop.fdr001, Cor_genes.top.200, Cor_genes.top.500, Cor_genes.top.1000, PCA_genes.top.200, PCA_genes.top.500, PCA_genes.top.1000)
names(genelist) <- c("miRNA7", "edgeR.diffexp.top.200.list", "edgeR.diffexp.fdr05.list", "edgeR.diffexp.fdr001.list", "deseq2.diffexp.top.200.list", "deseq2.diffexp.fdr05.list", "deseq2.diffexp.fdr001.list", "reads.qc.impute.HVG_genes.top.200", "reads.qc.impute.HVG_genes.fdr05", "reads.qc.impute.HVG_genes.fdr001", "reads.qc.impute.top_hvgs.200", "reads.qc.impute.top_hvgs.fdr05", "reads.qc.impute.top_hvgs.fdr001", "reads.qc.impute.top_m3drop.top.200", "reads.qc.impute.top_m3drop.fdr05", "reads.qc.impute.top_m3drop.fdr001", "Cor_genes.top.200", "Cor_genes.top.500", "Cor_genes.top.1000", "PCA_genes.top.200", "PCA_genes.top.500", "PCA_genes.top.1000")

save(genelist,file="genelist.Rdata")

pdf("reads.qc.impute.remove_confounders.selectgenes.heatmap.pdf")
for (i in 1:length(genelist)){
title <- names(genelist[i])
pheatmap(logcounts(reads.qc.impute)[rownames(logcounts(reads.qc.impute)) %in% genelist[[i]],], scale = "row", cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=group, main=title)
pheatmap(norm.ruvs20[rownames(norm.ruvs20) %in% genelist[[i]],], scale = "row", cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=group, main=title)
M3DropExpressionHeatmap(genelist[[i]], counts(reads.qc.impute), cell_labels = sampleLables)
}
dev.off()


## feature selection

inputF1 <- t(assay(reads.qc.impute, i="ruvs20"))
inputF1 <- apply(inputF1, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
sampleLables <- c()
for(i in colnames(reads.qc.impute)){tmp <- as.character(anno[which(anno$Sample_ID==i),"Class"]); sampleLables <- c(sampleLables,tmp)}
inputF1 <- data.frame(sampleID = rownames(inputF1), lables = sampleLables, inputF1)
rownames(inputF1) <- NULL


inputF2 <- deseq2.diffexp.fdr001.list
inputF2 <- rownames(reads.qc.impute) 

# lasso
########        lasso algorithm
#mx      =       read.table(inputF1,sep="\t",header=T);#the matrix in the fixed format
mx <- inputF1
mx      =       mx[,-1];
dataall =       as.matrix(mx[,-1]);
classesall      =       as.matrix(mx[,1]);
FI <- inputF2
#FI      =       as.vector(as.matrix(read.table(inputF2,head=F,sep="\t")));
FN      =       colnames(dataall);
ID      =       which(is.element(FN,FI));

dataall =       dataall[,ID];

library("glmnet")
cvob1   =       cv.glmnet(dataall,as.factor(classesall),family="multinomial");
coefs   <-      coef(cvob1,s="lambda.1se");
classesID       =       unique(classesall)[,1];

for(i in 1:length(classesID)){
        coefs_i =       as.matrix(coefs[[i]]);
        if(i==1){
                out     =       coefs_i;
        }else{
                out     =       cbind(out,coefs_i);
        }
}
colnames(out)   =       classesID;
out     =       data.frame(out);
coef    =       rownames(out);
coef    =       data.frame(coef);
out2    =       cbind(coef,out);

write.table(out2,file=outputF,quote=F,row.names=F,col.names=T,sep="\t");

##---------------------
## machine learning

library(SingleCellExperiment)
library(scater)
options(stringsAsFactors = FALSE)

reads.qc.impute <- readRDS("hcc_lulab.sequentialMap.homer.merged.clean.impute.rds")
load("genelist.Rdata")

# logistic regression in R
library(caret)

anno <- read.table("hcc_lulab.sample.info.txt", sep = "\t", header=T)
anno$Class <- "Ctrl"
anno[which(anno$Disease_status=="HCC_before_surgery"),"Class"] <- "HCC"
sampleLables <- c()
for(i in colnames(reads.qc.impute)){tmp <- as.character(anno[which(anno$Sample_ID==i),"Class"]); sampleLables <- c(sampleLables,tmp)}


# define training control
train_control<- trainControl(method="cv", number=10)
# train the model 
inputF <- t(assay(reads.qc.impute, i="ruvs20"))
inputF <- apply(inputF, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
inputF <- data.frame(sampleID = rownames(inputF), lables = sampleLables, inputF)
rownames(inputF) <- NULL

library(ROCR)


outlist <- list() 
pdf("reads.qc.impute.remove_confounders.selectgenes.glm.roc.pdf")
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
write.table(outmx, "./hcc_lulab.sequentialMap.homer.merged.glm.stat.txt", sep="\t",quote=F, row.names=F, col.names=T)


outlist <- list()
pdf("reads.qc.impute.remove_confounders.selectgenes.rf.roc.pdf")
for (i in 1:length(genelist)){
title <- names(genelist[i])

model<- train(inputF[,colnames(inputF) %in% genelist[[i]]], inputF[,2], method="rf", family=binomial())

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
write.table(outmx, "./hcc_lulab.sequentialMap.homer.merged.rf.stat.txt", sep="\t",quote=F, row.names=F, col.names=T)

# output genelist
for (i in 1:length(genelist)){
title <- names(genelist[i])
glist <- genelist[[i]] 
filename <- paste(title,".list.txt",sep="")
write.table(glist,filename,sep="\t",quote=F,row.names=F,col.names=F)
}

# bash
for i in `ls | grep 'list.txt' | grep -v "remove_confounders"`
do j=${i%.*};
awk 'BEGIN{FS=OFS="\t"}{split($1,a,"RNA_"); print $1,a[1]"RNA",a[2]}' ./$j.txt > ./$j.foo
mv ./$j.foo ./$j.txt;
cut -f 2 ./$j.txt | sort | uniq -c | awk 'BEGIN{FS=" "; OFS="\t"}{print $2,$1}' | sed -e "s/^/${j}\t/g" > ./$j.stat
done;






# print cv scores
summary(model)

pred = predict(model, newdata=inputF[,colnames(inputF) %in% deseq2.diffexp.fdr001.list])
confusionMatrix(data=pred, inputF$lables)

# roc 
library(ROCR)
# Compute AUC for predicting Class with the model
prob = predict(model, newdata=inputF[,colnames(inputF) %in% deseq2.diffexp.fdr001.list], type="prob")
pred = predict(model, newdata=inputF[,colnames(inputF) %in% deseq2.diffexp.fdr001.list])
perf <- performance(prediction(prob$HCC,pred), measure = "tpr", x.measure = "fpr")
plot(perf)
auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]
auc

pred <- prediction(pred, inputF$labels)
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
plot(perf, col=rainbow(10))


# random forest
##      train a model,
Rscript bin/train_model.rf.R input/GSE71008.reads_NvsS1.train.mx input/GSE71008.reads_NvsS1.FS      output/GSE71008.reads_NvsS1.model.rf
##      predict the matrix by the trained model
Rscript bin/predict_model.rf.R  output/GSE71008.reads_NvsS1.model.rf  input/GSE71008.reads_NvsS1.test.mx  input/GSE71008.reads_NvsS1.FS  output/GSE71008.reads_NvsS1.probs.rf
##      if the true labels are known, we can test the model by comparing the predicted label with the true label
cut -f 2 input/GSE71008.reads_NvsS1.test.mx > output/GSE71008.reads_NvsS1.true.labels
Rscript bin/test_model.R     output/GSE71008.reads_NvsS1.true.labels     output/GSE71008.reads_NvsS1.probs.rf output/GSE71008.reads_NvsS1.predicted.labels.rf      output/GSE71008.reads_NvsS1.performance.rf
##      get the ROC and PPV curve, not to draw it yet. using CDS as positive, others are negative
Rscript bin/get_ROC_PPV.R    output/GSE71008.reads_NvsS1.true.labels      output/GSE71008.reads_NvsS1.probs.rf "S1"   output/GSE71008.reads_NvsS1..S1.auc_value.rf output/GSE71008.reads_NvsS1.S1.roc_curve.rf output/GSE71008.reads_NvsS1.S1.ppv_curve.rf


 



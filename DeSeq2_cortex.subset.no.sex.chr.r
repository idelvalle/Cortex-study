
rm(list=ls())
setwd("~/Desktop/adrenal sex dimorphism study/RNAseq/deseq2_all/Cortex/")

library("DESeq2")
library(RColorBrewer)
library(gplots)
library(ggplot2)
library("ggbeeswarm")
library(reshape2)
library(plyr)
library("DESeq2")
library(RColorBrewer)
library(gplots)
library(ggplot2)
library("ggbeeswarm")
library(reshape2)
library(plyr)
library("dplyr")
library("magrittr")
library("pheatmap")
library("RColorBrewer")
library("genefilter")
library("clusterProfiler")
library(org.Hs.eg.db)

countdata <- read.table(file="counts.cortex.csv", header=TRUE, sep=",", row.names=1)

head(countdata)

colnames(countdata)

countdata$genes <- rownames(countdata)

colnames(countdata)

sex.chr <- read.table(file="X_and_Y.Chr.csv", header=FALSE, sep=",")

head(sex.chr)

subset <- countdata %>%
  filter(!(genes %in% sex.chr$V1))

head(subset)

rownames(subset) <- subset$genes

subset$genes <- NULL


phenodata <- read.table("phenodata.cortex.csv", sep=",", header=TRUE)
head(phenodata)

ncol(subset)
nrow(phenodata)
nrow(subset)
# Convert to matrix
countdata <- as.matrix(subset)
head(subset)

# Assign condition according to colnames
condition <- factor(phenodata[,11])
condition
stage<- factor(phenodata[,7])
stage
batch <- factor(phenodata[,5])
batch

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
coldata <- data.frame(row.names=colnames(subset), stage, condition, batch)
coldata

dds <- DESeqDataSetFromMatrix(countData=subset, colData=coldata, design=~ batch + condition)

dds

dds <- dds[ rowSums(counts(dds)) > 1, ]

dds <- DESeq(dds)

resultsNames(dds)

res <- results(dds)


stage.all.wpc <- results(dds, name="condition_XY_vs_XX",independentFiltering = TRUE, pAdjustMethod = "BH", cooksCutoff = FALSE) # To remove extreme outliers and avoid padj to NA
stage.all.wpc
summary(stage.all.wpc)

write.csv(as.data.frame(stage.all.wpc), file="results.cortex.subset.ALL.no.sex.csv")


# Plot dispersions

pdf("qc-dispersions.pdf")
plotDispEsts(dds, main="Dispersion plot")
dev.off()

# Regularized log transformation for clustering/heatmaps, etc. rlog tends to work well on small datasets (n <30)
rld <- rlogTransformation(dds, blind = TRUE)
head(assay(rld),3)
hist(assay(rld))

# vst is faster to compute and is less sensitive to high count outliers. Is recommended for large datasets (hundreds).
vsd <- vst(dds, blind=TRUE)
head(assay(vsd),3)
hist(assay(vsd))

# Sample Distances
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$stage, rld$condition, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pdf("Sample-distances.no.sex.chr.pdf")
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors, fontsize = 8)
dev.off()

# Is better to perform both transformations and compare the meanSdPlot

# Colors for plots using RColorBrewer
library(RColorBrewer)
mycols <- c("#F781BF","#377EB8")
mycols



# Principal components analysis

pcaData <- plotPCA(rld, returnData = TRUE)
pcaData

plotPCA.2_3 <- function (object, intgroup = "condition", ntop = 500, 
                         returnData = FALSE) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                               drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC2 = pca$x[, 2], PC3 = pca$x[, 3], group = group, 
                  intgroup.df, name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[2:3]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC2", y = "PC3", color = "group")) + 
    geom_point(size = 3) + xlab(paste0("PC2: ", round(percentVar[2] * 
                                                        100), "% variance")) + ylab(paste0("PC3: ", round(percentVar[3] * 
                                                                                                            100), "% variance")) + coord_fixed()
}

rv <- rowVars(assay(rld))
select <- order(rv, decreasing = TRUE)[seq_len(min(500,length(rv)))]
pca <- prcomp(t(assay(rld)[select, ]))
aload <- abs(pca$rotation)
loadings <- sweep(aload, 2, colSums(aload), "/")
colSums(sweep(aload, 2, colSums(aload), "/"))
write.table(loadings, file="loadings.no.sex.chr.csv", sep=",")

pcaData.2_3 <- plotPCA.2_3(rld, returnData = TRUE)

percentVar <- round(100 * attr(pcaData.2_3, "percentVar"))

pca.2_3 <- ggplot(pcaData.2_3, aes(x = PC2, y = PC3, color = group, shape=stage,label=name)) +
  geom_point(size=3) +
  xlab(paste0("PC2: ", percentVar[1], "% variance")) +
  ylab(paste0("PC3: ", percentVar[2], "% variance")) +
  coord_fixed()



pdf("pca.rld.pdf")
plotPCA(rld,intgroup= "stage")
dev.off()

pdf("pca.vsd.pdf")
plotPCA(vsd, intgroup= "stage")
dev.off()



percentVar <- round(100 * attr(pcaData, "percentVar"))

pca<-ggplot(pcaData, aes(x = PC1, y = PC2, color = group, shape=stage,label=name)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()

pdf("pca.2_3.no.sex.chr.pdf")
pca.2_3 + scale_color_manual(values=c("#ff8080","#99bbff"))
dev.off()


pdf("pca.no.sex.chr.pdf")
pca + scale_color_manual(values=c("#ff8080","#99bbff"))
dev.off()

pdf("pca_labels.pdf")
pca + geom_text(aes(label=phenodata$Number),hjust=0, vjust=0)+ scale_color_manual(values=c("#ff8080","#99bbff"))
dev.off()

tissue <- phenodata$Origin
pca.2<-ggplot(pcaData, aes(x = PC1, y = PC2, color = group, shape=tissue,label=name)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()

pdf("pca.tissue.pdf")
pca.2  + scale_color_manual(values=c("#ff8080","#99bbff"))
dev.off()

pdf("pca.tissue_labels.pdf")
pca.2 + geom_text(aes(label=phenodata$Number),hjust=0, vjust=0) + scale_color_manual(values=c("#ff8080","#99bbff"))
dev.off()

# Counts Plot

pdf("SRY.pdf")
plotCounts(dds, gene ="SRY")
dev.off()

pdf("TTTY14.pdf")
plotCounts(dds, gene ="TTTY14")
dev.off()

pdf("KCNC1.pdf")
plotCounts(dds, gene ="KCNC1")
dev.off()

pdf("CD99.pdf")
plotCounts(dds, gene ="CD99")
dev.off()

pdf("STS.pdf")
plotCounts(dds, gene ="STS")
dev.off()

pdf("TDRD12.pdf")
plotCounts(dds, gene ="TDRD12")
dev.off()

library("ggbeeswarm")
xist.counts <- plotCounts(dds, gene = "XIST",returnData = TRUE)
xist.counts
xist <- ggplot(xist.counts, aes(x = condition, y = count, color=condition)) + scale_y_log10() +  geom_beeswarm(cex = 3)
pdf("XIST.counts.pdf")
xist + scale_color_manual(values=c("#F781BF","#377EB8"))
dev.off()

# MA-plot

ma <- lfcShrink(dds, contrast =c("condition","XX","XY"))
pdf("MAplot.pdf")
plotMA(ma)
topGene <- rownames(res.05)[which.min(res.05$padj)]
with(res.05[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})
dev.off()

# Histogram of p-values. Exlude genes with very small counts (mean normalized count larger than 1), which otherwise generate spikes in the histogram

res <- results(dds)

pdf("histogram.p-values.pdf")
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")
dev.off()

# Gene Clustering, selecting the 20 genes with the highest variance across samples

mat  <- assay(rld)[ c("USP9Y",
                      "PRKY",
                      "DDX3Y",
                      "KDM5D",
                      "RPS4Y1",
                      "EIF1AY",
                      "UTY",
                      "TTTY15",
                      "NLGN4Y",
                      "ZFY",
                      "TMSB4Y",
                      "KALP",
                      "TBL1Y",
                      "TTTY14",
                      "PCDH11Y",
                      "TTTY10",
                      "SRY",
                      "OTOR",
                      "OGN",
                      "ASPN",
                      "ENTPD2",
                      "NOS2P3",
                      "EBF2",
                      "FAM157C",
                      "PTPN20A",
                      "NPIPA5",
                      "NAIP",
                      "SMN2",
                      "HERC2P3",
                      "FAM21B",
                      "PLGLA",
                      "PWRN3",
                      "ZNF812",
                      "FAM226B",
                      "AMBN",
                      "H2BFM",
                      "GSTM1",
                      "TSIX",
                      "FGF19",
                      "XIST"
),]

mat  <- mat - rowMeans(mat)

anno <- as.data.frame(colData(rld)[, c("stage","condition")])
anno <- as.data.frame(phenodata$Sex)
row.names(anno) <- phenodata$Sample
colnames(anno) <- "condition"


condition <- c("#ff8080","#99bbff")
names(condition) <- c("XX", "XY")
anno_colors <- list(condition = condition)

pdf("heatmap.20-top.genes.pdf")
pheatmap(mat, annotation_col = anno, annotation_colors = anno_colors, fontsize_row=8, fontsize_col=4,cluster_cols = TRUE, cluster_rows = FALSE)
dev.off()

#clusterProfiler

mcols(stage.all.wpc)

resSig <- subset(stage.all.wpc, padj < 0.05)

sig <- as.data.frame((resSig))

universe <- as.data.frame(res)
nrow(universe)

keytypes(org.Hs.eg.db)

gene.df <- bitr(rownames(sig), fromType = "SYMBOL",toType = "ENTREZID", OrgDb = org.Hs.eg.db)

View(gene.df)

#ggo <- groupGO(gene = gene.df$ENTREZID , OrgDb    = org.Hs.eg.db, ont      = "CC", level    = 2, readable = TRUE)

#head(ggo)

ggo <- groupGO(gene = rownames(sig), OrgDb = org.Hs.eg.db, keytype = "SYMBOL", ont = "BP", level = 2)

head(ggo)

#GO overrepresentattion test

ego <- enrichGO(gene = rownames(sig), universe = rownames(universe), OrgDb = org.Hs.eg.db, keytype = "SYMBOL", ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)

head(ego)
ego.s <- simplify(ego)

pdf("Adrenal.XY.vs.XX.ego.pdf")
barplot(ego, drop=TRUE, showCategory=20)
dev.off()

pdf("Adrenal.XY.vs.XX.ggo.pdf")
barplot(ggo, drop=TRUE, showCategory=20)
dev.off()


pdf("Adrenal.XY.vs.XX.dotplot.pdf")
dotplot(ego.s)
dev.off()

enrichMap(simplify(ego))


save.image(file="Cortex.subset.XY.vs.XX.NO.sex.chr.RData")




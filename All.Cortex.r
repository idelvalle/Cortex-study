
rm(list=ls())
setwd("~/Desktop/RNA-Seq-Brain-June-2018/Cortex/All.Cortex/")

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
library("AnnotationDbi")
library("Homo.sapiens")
library("limma")
library("sva")

countdata <- read.table(file="counts.cortex.csv",header=TRUE, sep=",", row.names=1)

head(countdata)

colnames(countdata)

phenodata <- read.table("phenodata.cortex.csv", sep=",", header=TRUE, row.names = 1)
phenodata

ncol(countdata)
nrow(phenodata)

colnames(countdata) <- phenodata$Sequencing_ID

head(phenodata)

phenodata.ctx <- phenodata %>% filter(Tissue=="cortex")
phenodata.ctx
countdata.ctx <- countdata[,phenodata.ctx$Sequencing_ID]
  
# Convert to matrix
countdata.ctx <- as.matrix(countdata.ctx)
head(countdata.ctx)

# Assign condition according to colnames

sex <- factor(phenodata.ctx$Sex)
sex


# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
coldata <- data.frame(row.names=colnames(countdata.ctx), sex)
coldata


######### ALL CORTEX XY vs XX ########################

dds <- DESeqDataSetFromMatrix(countData=countdata.ctx, colData=coldata, design=~sex)

dds

dds <- dds[ rowSums(counts(dds)) > 10, ]

dds <- DESeq(dds)

#BUILDING RESULTS TABLE
resultsNames(dds)

res <- results(dds, name = "sex_XY_vs_XX", cooksCutoff = FALSE)


mcols(res, use.names=TRUE)
summary(res)

columns(Homo.sapiens)

res.symbol <- res

res.symbol$symbol <- mapIds(Homo.sapiens,keys=row.names(res),column="SYMBOL",keytype="ENSEMBL",multiVals="first")

res.symbol <- res.symbol[complete.cases(res.symbol[, 6:7]),] #Remove NA values from Gene Symbol & padj
nrow(res.symbol)

res.dup <- res.symbol[duplicated(res.symbol$symbol),] #Visualize duplicates

nrow(res.dup)

write.csv(res.symbol, file="results.Cortex.XY.vs.XX.csv")
write.csv(res.dup, file="results.Cortex.XY.vs.XX.dup.csv")


################### MA PLOT ##############################################################################################


DESeq2::plotMA(res, alpha = 0.05, ylim=c(-10,10)) #Points colored in red if padj is lower than 0.05.Points which fall out of the window are plotted as open triangles pointing either up or down.



########### RLD TRANSFORMATION ##################

# Regularized log transformation for clustering/heatmaps, etc. rlog tends to work well on small datasets (n <30)
rld <- rlogTransformation(dds, blind = TRUE)
head(assay(rld),3)
hist(assay(rld))
names(colData(dds))

plotPCA(rld, intgroup="sex")


pcaData.sex <- plotPCA(rld, intgroup="sex", returnData=TRUE)
percentVar <- round(100 * attr(pcaData.sex, "percentVar"))
ggplot(pcaData.sex, aes(PC1, PC2, color=sex, shape=sex)) +
  geom_point(size=3) + geom_text(aes(label=phenodata.ctx$Sequencing_ID),hjust=0.5, vjust=1.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()



#### SAMPLE DISTANCES #####

# Sample Distances
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- colnames(rld)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

# Plot dispersions

pdf("qc-dispersions.pdf")
plotDispEsts(dds, main="Dispersion plot")
dev.off()


##########################################################################################################################

########### SAMPLE CLUSTERING AND VISUALIZATION #################

data <- as.data.frame((assay(rld)))
nrow(data) #32276
data$symbol <- mapIds(Homo.sapiens,keys=row.names(data),column="SYMBOL",keytype="ENSEMBL",multiVals="first")
data <- data[complete.cases(data[, 68]),] #Remove NA values from Gene Symbol
nrow(data) #20596
head(data)
data <- data[!duplicated(data$symbol),]
nrow(data) #20553
rownames(data) <- data$symbol
data <- data[,-68]


data.adrenal <- data[,c(31:67)]

steroido <- c("LDLR","SCARB1","DHCR24","STAR","CYP11A1", "FDX1", "HSD3B2", "CYP17A1", "PAPSS2", "FOXO4", "MAP3K15", "INHA", "MGARP",
  "RMDN2", "GRAMD1B", "SLC16A9", "SLC8B1")
steroido <- c("CYP19A1","NR5A1", "HSD17B3", "CYP17A1", "CYP11A1", "HSD3B2", "STAR", "DHCR24", "MRAP", "POR", "CYB5A", "MC2R", "FDX1", "SULT2A1", "CYP11B1", "CYP21A2", "SOAT1", "LDLR", "SCARB1", "PAPSS2")
steroido.data <- subset(data.adrenal, row.names(data.adrenal) %in% steroido)

pheatmap(steroido.data, cluster_rows = TRUE, cluster_cols = FALSE, fontsize =8, width=7, height=6)





res.ordered <- res.symbol[order(-res.symbol$log2FoldChange),]
head(res.ordered)
selected <- res.ordered[1:20,]
rld.df <- as.data.frame(assay(rld))
                        
t <- subset(rld.df, row.names(rld.df) %in% row.names(selected))
t[1:5,]

t <- t[rownames(selected),] #### IMPORTANT TO MAINTAIN ORDER AFTER MATCHING



t <- t[, col_order]

t$symbol <- mapIds(Homo.sapiens,keys=row.names(t),column="SYMBOL",keytype="ENSEMBL",multiVals="asNA")
t
row.names(t) <- t$symbol
t <- t[,-68]


pheatmap(t, cluster_rows = FALSE, cluster_cols = FALSE, fontsize =8, width=7, height=6)


topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 20)
topVarGenes
mat  <- as.data.frame(assay(rld)[ topVarGenes, ])
mat$symbol <- mapIds(Homo.sapiens,keys=row.names(mat),column="SYMBOL",keytype="ENSEMBL",multiVals="first")
mat
row.names(mat) <- mat$symbol
mat <- mat[,-68]

mat  <- mat - rowMeans(mat)
mat
mat <- mat[, col_order]

anno <- as.data.frame(colData(vsd)[, c("batch","tissue")])

pdf("heatmap.20-top.genes.pdf")
pheatmap(mat, cluster_rows = TRUE, cluster_cols = TRUE, annotation_col = anno, fontsize =8, width=7, height=6)
dev.off()




















# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix

dds.time <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~batch + condition + stage + condition:stage)

dds.time

dds.time <- dds.time[ rowSums(counts(dds.time)) > 1, ]

dds.time <- DESeq(dds.time, test="LRT", reduced = ~batch + condition + stage)

resTC <- results(dds.time)
mcols(resTC, use.names=TRUE)
summary(resTC)

resTC$symbol <- mcols(dds.time)$symbol
head(resTC[order(resTC$padj),], 4)

fiss <- plotCounts(dds.time, which.min(resTC$log2FoldChange), 
                   intgroup = c("condition","stage"), returnData = TRUE)
ggplot(fiss,
       aes(x = stage, y = count, color = condition, group = condition)) + 
  geom_point() + geom_smooth(se = FALSE, method = "loess") + scale_y_log10()

#BUILDING RESULTS TABLE
resultsNames(dds.time)

betas <- coef(dds.time)
colnames(betas)

topGenes <- head(order(resTC$padj),20)
mat <- betas[topGenes, -c(1,2)]
thr <- 3 
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=FALSE)

res <- results(dds)
mcols(res, use.names=TRUE)
summary(res)

resLFC1 <- results(dds, lfcThreshold=2)
table(resLFC1$padj < 0.05) #1037 Genes pass the Filter

columns(Homo.sapiens)

res$symbol <- mapIds(Homo.sapiens,keys=row.names(res),column="SYMBOL",keytype="ENSEMBL",multiVals="first")

resOrdered <- res[order(res$padj),]
head(resOrdered)

resOrderedDF <- as.data.frame(resOrdered)
write.csv(resOrderedDF, file="results.Adrenal.vs.Control.csv")






















#clusterProfiler

mcols(res)

resSig <- subset(resOrdered, padj < 0.05 & log2FoldChange > 2 & baseMean>30)

sig <- as.data.frame((resSig))

universe <- as.data.frame(res)
nrow(universe)

keytypes(org.Hs.eg.db)

gene.df <- bitr(rownames(sig), fromType = "ENSEMBL",toType = "ENTREZID", OrgDb = org.Hs.eg.db)

View(gene.df)

ggo <- groupGO(gene = gene.df$ENTREZID , OrgDb    = org.Hs.eg.db, ont      = "CC", level    = 3, readable = TRUE)

head(ggo)

ggo <- groupGO(gene = rownames(sig), OrgDb = org.Hs.eg.db, keytype = "SYMBOL", ont = "BP", level = 2)

#GO overrepresentattion test

ego <- enrichGO(gene = rownames(sig), universe = rownames(universe), OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)

head(ego)

pdf("Adrenal.vs.Control.ego.pdf")
barplot(simplify(ego), drop=TRUE, showCategory=20)
dev.off()

pdf("Adrenal.vs.Control.ggo.pdf")
barplot(ggo, drop=TRUE, showCategory=20)
dev.off()

ego.s <- simplify(ego)

pdf("Adrenal.vs.Control.dotplot.pdf")
dotplot(ego.s)
dev.off()

pdf("Adrenal.vs.Control.enrichmap.pdf")
enrichMap(simplify(ego))
dev.off()







save.image(file="Adrenal.subset.RData")



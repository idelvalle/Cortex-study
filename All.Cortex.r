
rm(list=ls())
setwd("~/Desktop/Cortex-study/")

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

head(countdata)
head(phenodata)

# Assign condition according to colnames

sex <- factor(phenodata$Sex)
sex

stage <- factor(phenodata$Stage)
stage

tissue <- factor(phenodata$Tissue)
tissue

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix

coldata <- data.frame(row.names=colnames(countdata),sex,stage,tissue)
coldata


######### ALL SAMPLES XY vs XX ########################

dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design =~sex)

dds

dds <- dds[ rowSums(counts(dds)) > 10, ]

dds <- DESeq(dds)

colData(dds)

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

write.csv(res.symbol, file="results.ALL.XY.vs.XX.csv")
write.csv(res.dup, file="results.ALL.XY.vs.XX.dup.csv")

## MA PLOT ###


DESeq2::plotMA(res, alpha = 0.05, ylim=c(-10,10)) #Points colored in red if padj is lower than 0.05.Points which fall out of the window are plotted as open triangles pointing either up or down.



## RLD TRANSFORMATION ##

# Regularized log transformation for clustering/heatmaps, etc. rlog tends to work well on small datasets (n <30)
rld <- rlogTransformation(dds, blind = FALSE)
head(assay(rld),3)
hist(assay(rld))
names(colData(dds))

plotPCA(rld, intgroup="sex")


pcaData.sex <- plotPCA(rld, intgroup="sex", returnData=TRUE)
percentVar <- round(100 * attr(pcaData.sex, "percentVar"))
pca.Cortex <- ggplot(pcaData.sex, aes(PC1, PC2, color=sex, shape=tissue)) +
  geom_point(size=3) + geom_text(aes(label=phenodata$Stage),hjust=0.5, vjust=1.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

ggplot(pcaData.sex, aes(PC1, PC2, color=sex, shape=tissue)) +
  geom_point(size=3) +
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

all.sampleDistances <- pheatmap(sampleDistMatrix,
                                clustering_distance_rows = sampleDists,
                                clustering_distance_cols = sampleDists,
                                col = colors)

# Plot dispersions

pdf("qc-dispersions.pdf")
plotDispEsts(dds, main="Dispersion plot")
dev.off()


###############################################################################################################################




















######### ALL CORTEX XY vs XX ########################


phenodata.ctx <- phenodata %>% filter(Tissue=="cortex")
phenodata.ctx
countdata.ctx <- countdata[,phenodata.ctx$Sequencing_ID]
head(countdata.ctx)

# Convert to matrix
countdata.ctx <- as.matrix(countdata.ctx)
head(countdata.ctx)

# Assign condition according to colnames

sex <- factor(phenodata.ctx$Sex)
sex

stage <- factor(phenodata.ctx$Stage)
stage

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
coldata <- data.frame(row.names=colnames(countdata.ctx), sex, stage)
coldata



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

######### ALL FOREBRAIN XY vs XX ########################


phenodata.fb <- phenodata %>% filter(Tissue=="forebrain")
phenodata.fb
countdata.fb <- countdata[,phenodata.fb$Sequencing_ID]
head(countdata.fb)

# Convert to matrix
countdata.fb <- as.matrix(countdata.fb)
head(countdata.fb)

# Assign condition according to colnames

sex <- factor(phenodata.fb$Sex)
sex

stage <- factor(phenodata.fb$Stage)
stage

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
coldata <- data.frame(row.names=colnames(countdata.fb), sex, stage)
coldata



dds <- DESeqDataSetFromMatrix(countData=countdata.fb, colData=coldata, design=~sex)

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

write.csv(res.symbol, file="results.Forebrain.XY.vs.XX.csv")
write.csv(res.dup, file="results.Forebrain.2.XY.vs.XX.dup.csv")

######### ALL TELENCEPHALON XY vs XX ########################


phenodata.tl <- phenodata %>% filter(Tissue=="telencephalon")
phenodata.tl
countdata.tl <- countdata[,phenodata.tl$Sequencing_ID]

# Convert to matrix
countdata.tl <- as.matrix(countdata.tl)
head(countdata.tl)

# Assign condition according to colnames

sex <- factor(phenodata.tl$Sex)
sex

stage <- factor(phenodata.tl$Stage)
stage

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
coldata <- data.frame(row.names=colnames(countdata.tl), sex, stage)
coldata



dds <- DESeqDataSetFromMatrix(countData=countdata.tl, colData=coldata, design=~sex)

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

write.csv(res.symbol, file="results.Telencephalon.XY.vs.XX.csv")
write.csv(res.dup, file="results.Telencephalon.XY.vs.XX.dup.csv")

######### ALL BRAIN XY vs XX ########################


phenodata.br <- phenodata %>% filter(Tissue=="Brain")
phenodata.br
countdata.br <- countdata[,phenodata.br$Sequencing_ID]

# Convert to matrix
countdata.br <- as.matrix(countdata.br)
head(countdata.br)

# Assign condition according to colnames

sex <- factor(phenodata.br$Sex)
sex

stage <- factor(phenodata.br$Stage)
stage

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
coldata <- data.frame(row.names=colnames(countdata.br), sex, stage)
coldata



dds <- DESeqDataSetFromMatrix(countData=countdata.br, colData=coldata, design=~sex)

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

write.csv(res.symbol, file="results.Brain.XY.vs.XX.csv")
write.csv(res.dup, file="results.Brain.XY.vs.XX.dup.csv")




####### MF DESIGN ######

dds$group <- factor(paste0(dds$sex, dds$tissue))
design(dds) <- ~ group
dds <- DESeq(dds)
dds$group

#Forebrain XY vs XX

res.fb <- results(dds, contrast = c("group", "XYforebrain", "XXforebrain"), cooksCutoff = FALSE)

summary(res.fb)

res.fb.symbol <- res.fb
res.fb.symbol$symbol <- mapIds(Homo.sapiens,keys=row.names(res.fb),column="SYMBOL",keytype="ENSEMBL",multiVals="first")
res.fb.symbol <- res.fb.symbol[complete.cases(res.fb.symbol[, 6:7]),] #Remove NA values from Gene Symbol & padj
nrow(res.fb.symbol)
res.fb.dup <- res.fb.symbol[duplicated(res.fb.symbol$symbol),] #Visualize duplicates
nrow(res.fb.dup)

write.csv(res.fb.symbol, file="results.Forebrain.XY.vs.XX.csv")
write.csv(res.fb.dup, file="results.Forebrain.XY.vs.XX.dup.csv")

#Telencephalon XY vs XX

res.tl <- results(dds, contrast = c("group", "XYtelencephalon", "XXtelencephalon"), cooksCutoff = FALSE)

summary(res.tl)

res.tl.symbol <- res.tl
res.tl.symbol$symbol <- mapIds(Homo.sapiens,keys=row.names(res.tl),column="SYMBOL",keytype="ENSEMBL",multiVals="first")
res.tl.symbol <- res.tl.symbol[complete.cases(res.tl.symbol[, 6:7]),] #Remove NA values from Gene Symbol & padj
nrow(res.tl.symbol)
res.tl.dup <- res.tl.symbol[duplicated(res.tl.symbol$symbol),] #Visualize duplicates
nrow(res.tl.dup)

write.csv(res.tl.symbol, file="results.Telencephalon.XY.vs.XX.csv")
write.csv(res.tl.dup, file="results.Telencephalon.XY.vs.XX.dup.csv")


#Cortex XY vs XX

res.ctx <- results(dds, contrast = c("group", "XYcortex", "XXcortex"), cooksCutoff = FALSE)

summary(res.ctx)

res.ctx.symbol <- res.ctx
res.ctx.symbol$symbol <- mapIds(Homo.sapiens,keys=row.names(res.ctx),column="SYMBOL",keytype="ENSEMBL",multiVals="first")
res.ctx.symbol <- res.ctx.symbol[complete.cases(res.ctx.symbol[, 6:7]),] #Remove NA values from Gene Symbol & padj
nrow(res.ctx.symbol)
res.ctx.dup <- res.ctx.symbol[duplicated(res.ctx.symbol$symbol),] #Visualize duplicates
nrow(res.ctx.dup)

write.csv(res.ctx.symbol, file="results.Cortex.XY.vs.XX.csv")
write.csv(res.ctx.dup, file="results.Cortex.XY.vs.XX.dup.csv")


#Brain XY vs XX

res.br <- results(dds, contrast = c("group", "XYBrain", "XXBrain"), cooksCutoff = FALSE)

summary(res.br)

res.br.symbol <- res.br
res.br.symbol$symbol <- mapIds(Homo.sapiens,keys=row.names(res.br),column="SYMBOL",keytype="ENSEMBL",multiVals="first")
res.br.symbol <- res.br.symbol[complete.cases(res.br.symbol[, 6:7]),] #Remove NA values from Gene Symbol & padj
nrow(res.br.symbol)
res.br.dup <- res.br.symbol[duplicated(res.br.symbol$symbol),] #Visualize duplicates
nrow(res.br.dup)

write.csv(res.br.symbol, file="results.Brain.XY.vs.XX.csv")
write.csv(res.br.dup, file="results.Brain.XY.vs.XX.dup.csv")


save.image(file="Brain-study.RData")



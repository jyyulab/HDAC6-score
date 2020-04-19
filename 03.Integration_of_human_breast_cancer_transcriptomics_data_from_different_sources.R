## To integrate the 4 datasets of human breast cancer samples
## Coded by Qingfei Pan (Qingfei.Pan@stjude.org)
## R-3.6

## 0. Configuration
require(NetBID2)
setwd("/Users/qpan/Desktop/HDAC6Manuscript/CodesAndData/DATA")

## 1. read in the esets
# TCGA
load("./TCGA/TCGA_RNASeq.BRCA.log2CPM_geneLevel.eset")
tcga <- TCGA_RNASeq.BRCA.log2CPM_geneLevel.eset; rm(TCGA_RNASeq.BRCA.log2CPM_geneLevel.eset)

# METABRIC
load("./METABRIC/METABRIC_microArray.BRCA.log2Intensity_geneLevel.eset")
metabric <- METABRIC_microArray.BRCA.log2Intensity_geneLevel.eset; rm(METABRIC_microArray.BRCA.log2Intensity_geneLevel.eset)

# IBC
load("./IBC/GEP.197.quantileNormalized.log2.eset")
ibc <- GEP.197.quantileNormalized.log2.eset; rm(GEP.197.quantileNormalized.log2.eset)
fd <- fData(ibc)[,c("ID", "Gene Symbol")];
colnames(fd) <- c("ID", "geneSymbol"); head(fd)
d <- exprs(ibc)
exp <- merge(fd, d, by="row.names", all = T); exp <- exp[,-c(1:2)]
exp.gene <- aggregate(. ~ geneSymbol, exp, mean, na.rm = T)
exp.gene <- exp.gene[grepl("\\w", exp.gene$geneSymbol),]; dim(exp.gene)
row.names(exp.gene) <- exp.gene$geneSymbol
exp.ibc <- exp.gene[,-1]; dim(exp.ibc)

# Clinical
load("./ClinicalTrial/ClinicalTrial.log2CPM.genelevel.afterBatchEffectRemoval.eset")
clinical <- ClinicalTrial.log2CPM.genelevel.afterBatchEffectRemoval.eset; rm(ClinicalTrial.log2CPM.genelevel.afterBatchEffectRemoval.eset)

## 2. prepare the pd
pd.tcga <- pData(tcga); pd.tcga <- pd.tcga[,c(1,3,4)]; pd.tcga$Dataset <- "TCGA"
colnames(pd.tcga) <- c("SampleIndex", "Subtype1", "Subtype2", "Dataset")

pd.metabric <- pData(metabric); pd.metabric <- pd.metabric[,c(1,2,3)]; pd.metabric$Dataset <- "Metabric"
colnames(pd.metabric) <- c("SampleIndex", "Subtype1", "Subtype2", "Dataset")

pd.ibc <- pData(ibc); pd.ibc <- pd.ibc[,c(14,13)]; pd.ibc$Subtype2 <- paste0(pd.ibc$ER.IHC, ",", pd.ibc$PR.IHC); pd.ibc$Dataset <- "IBC"
colnames(pd.ibc) <- c("SampleIndex", "Subtype1", "Subtype2", "Dataset")

pd.trial <- pData(clinical); pd.trial <- pd.trial[,c(1,5,3)]; pd.trial$Dataset <- "Clinical"
colnames(pd.trial) <- c("SampleIndex", "Subtype1", "Subtype2", "Dataset")

pd.merge <- rbind(pd.tcga, pd.metabric, pd.ibc, pd.trial)

## 3. merge the exp
exp.tcga <- exprs(tcga); dim(exp.tcga)
exp.metabric <- exprs(metabric); dim(exp.metabric)
exp.clinical <- exprs(clinical); dim(exp.clinical)

exp1 <- merge(exp.tcga, exp.metabric, by = "row.names", all = F); dim(exp1)
row.names(exp1) <- exp1$Row.names; exp1 <- exp1[,-1]
exp2 <- merge(exp1, exp.ibc, by = "row.names", all = F); dim(exp2)
row.names(exp2) <- exp2$Row.names; exp2 <- exp2[,-1]
exp3 <- merge(exp2, exp.clinical, by = "row.names", all = F); dim(exp3)
row.names(exp3) <- exp3$Row.names; exp3 <- exp3[,-1]

# re-normalize
exp.merge <- RNASeqCount.normalize.scale(mat = exp3, total = 1000000, pseudoCount = 0); colSums(exp.merge)

## 4. prepare the fd
fd.merge <- data.frame(row.names = row.names(exp.merge), GeneSymbol = row.names(exp.merge))

## 5. generate the eset
eset <- new("ExpressionSet", phenoData = new("AnnotatedDataFrame", pd.merge), featureData = new("AnnotatedDataFrame", fd.merge),
            annotation = 'TCGA_METABRIC_IBC_CLINICAL', exprs = as.matrix(exp.merge))

FourSets.TCGA_METABRIC_IBC_CLINICAL.beforeBatchEffectRemoval.eset <- eset
save(FourSets.TCGA_METABRIC_IBC_CLINICAL.beforeBatchEffectRemoval.eset, file = "FourSets.TCGA_METABRIC_IBC_CLINICAL.beforeBatchEffectRemoval.eset")

## 6. remove batch effect
exp.final <- removeBatchEffect(exprs(eset), batch = pData(eset)$Dataset)
eset.clean <- new("ExpressionSet", phenoData = new("AnnotatedDataFrame", pData(eset)), featureData = new("AnnotatedDataFrame", fData(eset)), annotation = 'TCGA_METABRIC_IBC_CLINICAL', exprs = as.matrix(exp.final))

FourSets.TCGA_METABRIC_IBC_CLINICAL.afterBatchEffectRemoval.eset <- eset.clean
save(FourSets.TCGA_METABRIC_IBC_CLINICAL.afterBatchEffectRemoval.eset, file = "FourSets.TCGA_METABRIC_IBC_CLINICAL.afterBatchEffectRemoval.eset")

## 7. PCA
eset.sel <- FourSets.TCGA_METABRIC_IBC_CLINICAL.beforeBatchEffectRemoval.eset
eset.sel <- FourSets.TCGA_METABRIC_IBC_CLINICAL.afterBatchEffectRemoval.eset
pca <- prcomp(t(exprs(eset.sel)), scale = T) # no filtration on the eset
PC <- data.frame(row.names = row.names(pca$x), SampleIndex = row.names(pca$x), pca$x[, c("PC1", "PC2", "PC3")])
merged <- merge(pData(eset.sel), PC, by = "row.names")
xlab <- sprintf('PC1(%s%s variance)',format(summary(pca)$importance[2,'PC1']*100,digits=3),'%')
ylab <- sprintf('PC2(%s%s variance)',format(summary(pca)$importance[2,'PC2']*100,digits=3),'%')
p <- ggplot(merged, aes(x=PC1, y=PC2, color = Dataset))
p <- p + geom_point(shape=16, size = 1) + labs(x = xlab, y = ylab) +
  theme(legend.position = c(0.85,0.15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = "black", size = 1),
        axis.title = element_text(size = 15, face = "bold", colour = "black"),
        axis.text.x = element_text(size = 12, face = "bold", hjust = 0.5, color = "black"),
        axis.text.y = element_text(size = 12, face = "bold", hjust = 1, color = "black")) 
p



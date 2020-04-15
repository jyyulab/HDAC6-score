require(NetBID2)

## TCGA
load("/Users/qpan/Desktop/HDAC6Manuscript/clinicalTrial/Final/TCGA_BRCA.201503.26106_1221.log2counts.gene.eset")
tcga <- TCGA_BRCA.201503.26106_1221.log2counts.gene.eset; rm(TCGA_BRCA.201503.26106_1221.log2counts.gene.eset)

## METABRIC
load("/Users/qpan/Desktop/HDAC6Manuscript/clinicalTrial/Final/METABRIC.ALL_Tumor.24368_1904.Intensity.gene.eset")
metabric <- METABRIC.ALL_Tumor.24368_1904.Intensity.gene.eset; rm(METABRIC.ALL_Tumor.24368_1904.Intensity.gene.eset)

## IBC
load("/Users/qpan/Desktop/HDAC6Manuscript/clinicalTrial/Final/GEP.197.quantileNormalized.log2.eset")
ibc <- GEP.197.quantileNormalized.log2.eset; rm(GEP.197.quantileNormalized.log2.eset)
fd <- fData(ibc)[,c("ID", "Gene Symbol")];
colnames(fd) <- c("ID", "geneSymbol"); head(fd)
d <- exprs(ibc)
exp <- merge(fd, d, by="row.names", all = T); exp <- exp[,-c(1:2)]
exp.gene <- aggregate(. ~ geneSymbol, exp, mean, na.rm = T)
exp.gene <- exp.gene[grepl("\\w", exp.gene$geneSymbol),]; dim(exp.gene)
row.names(exp.gene) <- exp.gene$geneSymbol
exp.ibc <- exp.gene[,-1]; dim(exp.ibc)

## Clinical
load("/Users/qpan/Desktop/HDAC6Manuscript/clinicalTrial/Final/ClinicalTrial.log2TP50M.gene.50372_10.afterBatchEffectRemoval.eset")
clinical <- ClinicalTrial.log2TP50M.gene.50372_10.afterBatchEffectRemoval.eset; rm(ClinicalTrial.log2TP50M.gene.50372_10.afterBatchEffectRemoval.eset)

## pd
pd.tcga <- pData(tcga); pd.tcga <- pd.tcga[,c(1,3,4)]; pd.tcga$Dataset <- "TCGA"
colnames(pd.tcga) <- c("SampleIndex", "Subtype1", "Subtype2", "Dataset")

pd.metabric <- pData(metabric); pd.metabric <- pd.metabric[,c(3,1,2)]; pd.metabric$Dataset <- "Metabric"
colnames(pd.metabric) <- c("SampleIndex", "Subtype1", "Subtype2", "Dataset")

pd.ibc <- pData(ibc); pd.ibc <- pd.ibc[,c(14,13)]; pd.ibc$Subtype2 <- paste0(pd.ibc$ER.IHC, ",", pd.ibc$PR.IHC); pd.ibc$Dataset <- "IBC"
colnames(pd.ibc) <- c("SampleIndex", "Subtype1", "Subtype2", "Dataset")

pd.trial <- pData(clinical); pd.trial <- pd.trial[,c(1,5,6)]; pd.trial$Dataset <- "Clinical"
colnames(pd.trial) <- c("SampleIndex", "Subtype1", "Subtype2", "Dataset")

pd.merge <- rbind(pd.tcga, pd.metabric, pd.ibc, pd.trial)

## exp
exp.tcga <- exprs(tcga); dim(exp.tcga)
exp.metabric <- exprs(metabric); dim(exp.metabric)
#exp.ibc <- exprs(ibc); dim(exp.ibc)
exp.clinical <- exprs(clinical); dim(exp.clinical)

exp1 <- merge(exp.tcga, exp.metabric, by = "row.names", all = F); dim(exp1)
row.names(exp1) <- exp1$Row.names; exp1 <- exp1[,-1]
exp2 <- merge(exp1, exp.ibc, by = "row.names", all = F); dim(exp2)
row.names(exp2) <- exp2$Row.names; exp2 <- exp2[,-1]
exp3 <- merge(exp2, exp.clinical, by = "row.names", all = F); dim(exp3)
row.names(exp3) <- exp3$Row.names; exp3 <- exp3[,-1]

exp.merge <- RNASeqCount.normalize.scale(mat = exp3, total = 1000000, pseudoCount = 0); colSums(exp.merge)

## fd
fd.merge <- data.frame(row.names = row.names(exp.merge), GeneSymbol = row.names(exp.merge))

##
eset <- new("ExpressionSet", phenoData = new("AnnotatedDataFrame", pd.merge), featureData = new("AnnotatedDataFrame", fd.merge),
            annotation = 'TCGA_METABRIC_IBC_CLINICAL', exprs = as.matrix(exp.merge))

FourSets.TCGA_METABRIC_IBC_CLINICAL.beforeBatchEffectRemoval.14719_3332.eset <- eset
save(FourSets.TCGA_METABRIC_IBC_CLINICAL.beforeBatchEffectRemoval.14719_3332.eset, file = "/Users/qpan/Desktop/HDAC6Manuscript/clinicalTrial/Final/FourSets.TCGA_METABRIC_IBC_CLINICAL.beforeBatchEffectRemoval.14719_3332.eset")


exp.final <- removeBatchEffect(exprs(FourSets.TCGA_METABRIC_IBC_CLINICAL.beforeBatchEffectRemoval.14719_3332.eset), batch = pData(FourSets.TCGA_METABRIC_IBC_CLINICAL.beforeBatchEffectRemoval.14719_3332.eset)$Dataset)
eset <- new("ExpressionSet", phenoData = new("AnnotatedDataFrame", pData(FourSets.TCGA_METABRIC_IBC_CLINICAL.beforeBatchEffectRemoval.14719_3332.eset)), featureData = new("AnnotatedDataFrame", fData(FourSets.TCGA_METABRIC_IBC_CLINICAL.beforeBatchEffectRemoval.14719_3332.eset)),
            annotation = 'TCGA_METABRIC_IBC_CLINICAL', exprs = as.matrix(exp.final))

FourSets.TCGA_METABRIC_IBC_CLINICAL.afterBatchEffectRemoval.14719_3332.eset <- eset
save(FourSets.TCGA_METABRIC_IBC_CLINICAL.afterBatchEffectRemoval.14719_3332.eset, file = "/Users/qpan/Desktop/HDAC6Manuscript/clinicalTrial/Final/FourSets.TCGA_METABRIC_IBC_CLINICAL.afterBatchEffectRemoval.14719_3332.eset")


eset.sel <- FourSets.TCGA_METABRIC_IBC_CLINICAL.beforeBatchEffectRemoval.14719_3332.eset
eset.sel <- FourSets.TCGA_METABRIC_IBC_CLINICAL.afterBatchEffectRemoval.14719_3332.eset
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
ggsave("/Users/qpan/Desktop/HDAC6Manuscript/clinicalTrial/Final/03_PCA_4datasets.afterBatchEffectRemoval.pdf", p, width = 7, height = 7, units = "in", useDingbats = F)

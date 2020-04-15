require(NetBID2)

## TCGA
load("/Users/qpan/Desktop/HDAC6Manuscript/clinicalTrial/Final/TCGA_BRCA.201503.26106_1221.log2counts.gene.eset")
tcga <- TCGA_BRCA.201503.26106_1221.log2counts.gene.eset; rm(TCGA_BRCA.201503.26106_1221.log2counts.gene.eset)

## METABRIC
load("/Users/qpan/Desktop/HDAC6Manuscript/clinicalTrial/Final/METABRIC.ALL_Tumor.24368_1904.Intensity.gene.eset")
metabric <- METABRIC.ALL_Tumor.24368_1904.Intensity.gene.eset; rm(METABRIC.ALL_Tumor.24368_1904.Intensity.gene.eset)

## pd
pd.tcga <- pData(tcga); pd.tcga <- pd.tcga[,c(1,3,4)]; pd.tcga$Dataset <- "TCGA"
colnames(pd.tcga) <- c("SampleIndex", "Subtype1", "Subtype2", "Dataset")

pd.metabric <- pData(metabric); pd.metabric <- pd.metabric[,c(3,1,2)]; pd.metabric$Dataset <- "Metabric"
colnames(pd.metabric) <- c("SampleIndex", "Subtype1", "Subtype2", "Dataset")

pd.merge <- rbind(pd.tcga, pd.metabric)

## exp
exp.tcga <- exprs(tcga); dim(exp.tcga)
exp.metabric <- exprs(metabric); dim(exp.metabric)

exp1 <- merge(exp.tcga, exp.metabric, by = "row.names", all = F); dim(exp1)
row.names(exp1) <- exp1$Row.names; exp1 <- exp1[,-1]

exp.merge <- RNASeqCount.normalize.scale(mat = exp1, total = 1000000, pseudoCount = 0); colSums(exp.merge)

## fd
fd.merge <- data.frame(row.names = row.names(exp.merge), GeneSymbol = row.names(exp.merge))

##
eset <- new("ExpressionSet", phenoData = new("AnnotatedDataFrame", pd.merge), featureData = new("AnnotatedDataFrame", fd.merge),
            annotation = 'TCGA_METABRIC', exprs = as.matrix(exp.merge))

TwoSets.TCGA_METABRICL.beforeBatchEffectRemoval.17115_3125.eset <- eset
save(TwoSets.TCGA_METABRICL.beforeBatchEffectRemoval.17115_3125.eset, file = "/Users/qpan/Desktop/HDAC6Manuscript/clinicalTrial/Final/TwoSets.TCGA_METABRICL.beforeBatchEffectRemoval.17115_3125.eset")


exp.final <- removeBatchEffect(exprs(TwoSets.TCGA_METABRICL.beforeBatchEffectRemoval.17115_3125.eset), batch = pData(TwoSets.TCGA_METABRICL.beforeBatchEffectRemoval.17115_3125.eset)$Dataset)
eset <- new("ExpressionSet", phenoData = new("AnnotatedDataFrame", pData(TwoSets.TCGA_METABRICL.beforeBatchEffectRemoval.17115_3125.eset)), featureData = new("AnnotatedDataFrame", fData(TwoSets.TCGA_METABRICL.beforeBatchEffectRemoval.17115_3125.eset)),
            annotation = 'TCGA_METABRIC', exprs = as.matrix(exp.final))

TwoSets.TCGA_METABRIC.afterBatchEffectRemoval.17115_3125.eset <- eset
save(TwoSets.TCGA_METABRIC.afterBatchEffectRemoval.17115_3125.eset, file = "/Users/qpan/Desktop/HDAC6Manuscript/clinicalTrial/Final/TwoSets.TCGA_METABRIC.afterBatchEffectRemoval.17115_3125.eset")

## load the regulon
regulon.file <- read.table("/Volumes/yu3grp/solidTumor_JY/yu3grp/BRCA/HDAC6/Qingfei/04_Datasets/06_finalFigures/15_compare_regulon_new.txt", header = F, sep = "\t", stringsAsFactors = F)
regulon <- data.frame(row.names = regulon.file$V1, target = regulon.file$V1); head(regulon)
geneset <- list(); geneset[[1]] <- regulon; names(geneset) <- "HDAC6regulon"

## calculate the HDAC6 score
exp <- exprs(eset)
hdac6score <- cal.Activity(target_list = geneset, cal_mat = as.matrix(exp),
                           es.method = 'mean', # 'Weightedmean', 'mean', 'maxmean', 'absmean'
                           std = T, memory_constrain = F # if true, the calculation strategy will not use Matrix Cross Products, which is memory consuming.
)
hdac6score <- data.frame(t(hdac6score))

std <- function(x){
  tmp_mean <- mean(x, na.rm = T); tmp_sd <- sd(x, na.rm = T); (x - tmp_mean) / tmp_sd
}
hdac6score$HDAC6regulon <- std(hdac6score$HDAC6regulon)

## pd
pd <- pData(eset)
master <- merge(pd, hdac6score, by = "row.names", all = T); dim(master)
table(master$Subtype1); table(master$Subtype2);

#########################################################
for (i in 1:nrow(master)) {
  if (master$Subtype1[i] == "IBC") { master$Subtype2[i] <- "IBC" }
  else if (master$Subtype1[i] == "NIBC") { master$Subtype2[i] <- "NIBC" }
  else if (master$Subtype2[i] == "Non-responder") { master$Subtype1[i] <- "Non-responder" }
  else if (master$Subtype2[i] == "Responder") { master$Subtype1[i] <- "Responder" }
}
master$Subtype1[master$Subtype1 == "HRpositive_HP"] <- "HRpositive"
master$Subtype1[master$Subtype1 == "HRpositive_LP"] <- "HRpositive"

write.table(master, file = "/Users/qpan/Desktop/HDAC6Manuscript/clinicalTrial/Final/03_2datasets_tcgametabric_HDAC6score.masterTable.txt",
            col.names = T, row.names = F, sep = "\t", quote = F)

# HRpositive, HER2positive, TNBC, normalSamples
c1 <- master$HDAC6regulon[master$Subtype1 == "HRpositive"]
c2 <- master$HDAC6regulon[master$Subtype1 == "HER2positive"]
t.test(c1,c2)



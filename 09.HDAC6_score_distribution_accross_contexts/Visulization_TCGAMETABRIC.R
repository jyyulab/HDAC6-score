## To visulize the HDAC6 score distribution among datasets
## Coded by Qingfei Pan (Qingfei.Pan@stjude.org)
## R-3.6

## 0. configuration
require(NetBID2)
setwd("./DATA")

## 1. read in the esets
# TCGA
load("./TCGA/TCGA_RNASeq.BRCA.log2CPM_geneLevel.eset")
tcga <- TCGA_RNASeq.BRCA.log2CPM_geneLevel.eset; rm(TCGA_RNASeq.BRCA.log2CPM_geneLevel.eset)

# METABRIC
load("./METABRIC/METABRIC_microArray.BRCA.log2Intensity_geneLevel.eset")
metabric <- METABRIC_microArray.BRCA.log2Intensity_geneLevel.eset; rm(METABRIC_microArray.BRCA.log2Intensity_geneLevel.eset)

## 2. prepare the pd
pd.tcga <- pData(tcga); pd.tcga <- pd.tcga[,c(1,3,4)]; pd.tcga$Dataset <- "TCGA"
colnames(pd.tcga) <- c("SampleIndex", "Subtype1", "Subtype2", "Dataset")

pd.metabric <- pData(metabric); pd.metabric <- pd.metabric[,c(1,2,3)]; pd.metabric$Dataset <- "Metabric"
colnames(pd.metabric) <- c("SampleIndex", "Subtype1", "Subtype2", "Dataset")

pd.merge <- rbind(pd.tcga, pd.metabric)

## 3. merge the exp
exp.tcga <- exprs(tcga); dim(exp.tcga)
exp.metabric <- exprs(metabric); dim(exp.metabric)

exp1 <- merge(exp.tcga, exp.metabric, by = "row.names", all = F); dim(exp1)
row.names(exp1) <- exp1$Row.names; exp1 <- exp1[,-1]

# re-normalize
exp.merge <- RNASeqCount.normalize.scale(mat = exp1, total = 1000000, pseudoCount = 0); colSums(exp.merge)

## 4. prepare the fd
fd.merge <- data.frame(row.names = row.names(exp.merge), GeneSymbol = row.names(exp.merge))

## 5. generate the eset
eset <- new("ExpressionSet", phenoData = new("AnnotatedDataFrame", pd.merge), featureData = new("AnnotatedDataFrame", fd.merge),
            annotation = 'TCGA_METABRIC', exprs = as.matrix(exp.merge))

TwoSets.TCGA_METABRICL.beforeBatchEffectRemoval.eset <- eset
save(TwoSets.TCGA_METABRICL.beforeBatchEffectRemoval.eset, file = "TwoSets.TCGA_METABRICL.beforeBatchEffectRemoval.eset")


exp.final <- removeBatchEffect(exprs(eset), batch = pData(eset)$Dataset)
eset <- new("ExpressionSet", phenoData = new("AnnotatedDataFrame", pData(eset)), featureData = new("AnnotatedDataFrame", fData(eset)),
            annotation = 'TCGA_METABRIC', exprs = as.matrix(exp.final))

TwoSets.TCGA_METABRIC.afterBatchEffectRemoval.eset <- eset
save(TwoSets.TCGA_METABRIC.afterBatchEffectRemoval.eset, file = "TwoSets.TCGA_METABRIC.afterBatchEffectRemoval.eset")



## 6. load the regulon
regulon.file <- read.table("./HDAC6_Breast_Cancer_Regulon.txt", header = F, sep = "\t", stringsAsFactors = F)
regulon <- data.frame(row.names = regulon.file$V1, target = regulon.file$V1); head(regulon)
geneset <- list(); geneset[[1]] <- regulon; names(geneset) <- "HDAC6regulon"

## 7. calculate the HDAC6 score
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

## 8. generate the eset
d <- t(hdac6score)
pd <- pData(eset)[colnames(d),]
fd <- data.frame(row.names = row.names(d), Regulon = row.names(d))
eset.act <- generate.eset(exp_mat = d, phenotype_info = pd, feature_info = fd, annotation_info = 'TCGA_METABRIC, GeneLevel, Activity')
TCGA_METABRIC.activity.eset <- eset.act
save(TCGA_METABRIC.activity.eset, file = "TCGA_METABRIC.activity.eset")

## 9. visulization
eset <- TCGA_METABRIC.activity.eset; rm(TCGA_METABRIC.activity.eset)

d <- exprs(eset)
pd <- pData(eset)
master <- merge(pd, t(d), by = "row.names", all = T); dim(master)

# subtype1
table(master$Subtype1)
for (i in 1:nrow(master)) {
  if (master$Subtype1[i] == "HRpositive_HP") {
    master$Subtype1[i] <- "HRpositive"
  } else if (master$Subtype1[i] == "HRpositive_LP") {
    master$Subtype1[i] <- "HRpositive"
  } else if (master$Subtype1[i] == "tripleNegative") {
    master$Subtype1[i] <- "TNBC"
  } else if (master$Subtype1[i] == "normalSamples") {
    master$Subtype1[i] <- "Paired-Normal"
  }
}

master.ihc <- master[master$Subtype1 != "NA",]
models <- unique(master.ihc$Subtype1)
median_d <- data.frame()
for (i in 1:length(models)) {
  tag <- as.character(models[i])
  median_d[i,1] <- tag
  data.tmp <- master.ihc[master.ihc$Subtype1 == tag, ]
  median_d[i,2] <- median(data.tmp$HDAC6regulon, na.rm = T)
}
median_sorted <- median_d[order(-median_d$V2), ]
model_order1 <- median_sorted$V1
p <- ggplot(master.ihc, aes(x=Subtype1, y=HDAC6regulon))
p <- p + geom_boxplot(aes(fill=Subtype1), outlier.shape = NA) + labs(x = "", y = "HDAC6 Score") + ylim(-4,4) + 
  geom_jitter(shape=16, position=position_jitter(0.3), size = 0.6) +
  theme(legend.position = "",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = "black", size = 1),
        axis.title = element_text(size = 15, face = "bold", colour = "black"),
        axis.text.x = element_text(size = 15, face = "bold", hjust = 1, color = "black", angle = 45),
        axis.text.y = element_text(size = 15, face = "bold", hjust = 1, color = "black")) +
  scale_x_discrete(limits = model_order1, labels=c("Paired-Normal" = "Paired-Normal\n(n=111)", "HER2positive" = "HER2+/(HR+ & HR-)\n(n=296)",
                                                   "HRpositive" = "HR+/HER2-\n(n=1708)", "TNBC" = "TNBC\n(n=429)")) 
p

# subtype2
table(master$Subtype2)
master.mol <- master[master$Subtype2 != "NA",]
models <- unique(master.mol$Subtype2)
median_d <- data.frame()
for (i in 1:length(models)) {
  tag <- as.character(models[i])
  median_d[i,1] <- tag
  data.tmp <- master.mol[master.mol$Subtype2 == tag, ]
  median_d[i,2] <- median(data.tmp$HDAC6regulon, na.rm = T)
}
median_sorted <- median_d[order(-median_d$V2), ]
model_order1 <- median_sorted$V1
p <- ggplot(master.mol, aes(x=Subtype2, y=HDAC6regulon))
p <- p + geom_boxplot(aes(fill=Subtype2), outlier.shape = NA) + labs(x = "", y = "HDAC6 Score") + ylim(-4,4) + 
  geom_jitter(shape=16, position=position_jitter(0.3), size = 0.6) +
  theme(legend.position = "",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = "black", size = 1),
        axis.title = element_text(size = 15, face = "bold", colour = "black"),
        axis.text.x = element_text(size = 15, face = "bold", hjust = 1, color = "black", angle = 45),
        axis.text.y = element_text(size = 15, face = "bold", hjust = 1, color = "black")) +
  scale_x_discrete(limits = model_order1, labels=c("Paired-Normal" = "Paired-Normal\n(n=111)", "Claudin_Low" = "Claudin-low\n(n=199)",
                                                   "Basal" = "Basal-like\n(n=385)", "Normal" = "Normal-like\n(n=177)",
                                                   "Her2" = "Her2-enriched\n(n=298)", "LumA" = "Luminal-A\n(n=1217)", "LumB" = "Luminal-B\n(n=664)")) 
p

## To calculate the HDAC6 score of TCGA
## Coded by Qingfei Pan (Qingfei.Pan@stjude.org)
## R-3.6

## 0. configuration
library(NetBID2)
setwd("./DATA/TCGA")

## 1. load the eset
load("TCGA_RNASeq.BRCA.log2FPKM_geneLevel.eset")
eset <- TCGA_RNASeq.BRCA.log2FPKM_geneLevel.eset; rm(TCGA_RNASeq.BRCA.log2FPKM_geneLevel.eset)

## 2. read the regulon
regulon.file <- read.table("./HDAC6_Breast_Cancer_Regulon.txt", header = F, sep = "\t", stringsAsFactors = F)
regulon <- data.frame(row.names = regulon.file$V1, target = regulon.file$V1); head(regulon)
geneset <- list(); geneset[[1]] <- regulon; names(geneset) <- "HDAC6regulon"

## 3. calculate the HDAC6 score
exp <- exprs(eset)
hdac6score <- cal.Activity(target_list = geneset, cal_mat = as.matrix(exp),
                           es.method = 'mean', # 'Weightedmean', 'mean', 'maxmean', 'absmean'
                           std = T, memory_constrain = F # if true, the calculation strategy will not use Matrix Cross Products, which is memory consuming.
)
hdac6score <- data.frame(t(hdac6score))

# z-normalize the raw HDAC6 score
std <- function(x){
  tmp_mean <- mean(x, na.rm = T); tmp_sd <- sd(x, na.rm = T); (x - tmp_mean) / tmp_sd
}
hdac6score$HDAC6regulon <- std(hdac6score$HDAC6regulon)

## 4. generate the eset
d <- t(hdac6score)
pd <- pData(eset)[colnames(d),]
fd <- data.frame(row.names = row.names(d), Regulon = row.names(d))
eset.act <- generate.eset(exp_mat = d, phenotype_info = pd, feature_info = fd, annotation_info = 'TCGA, GeneLevel, Activity')
TCGA.activity.eset <- eset.act
save(TCGA.activity.eset, file = "TCGA.activity.eset")

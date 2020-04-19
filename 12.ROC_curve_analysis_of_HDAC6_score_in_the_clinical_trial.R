## ROC curve fitting
## Coded by Qingfei Pan (Qingfei.Pan@stjude.org)
## R-3.6

## 0. configuration
require(NetBID2)
library(pROC)
require(ggplot2)
setwd("./DATA/ClinicalTrial")

## 1. load the eset
load("ClinicalTrial.log2CP50M.genelevel.afterBatchEffectRemoval.eset")
eset <- ClinicalTrial.log2CP50M.genelevel.afterBatchEffectRemoval.eset; rm(ClinicalTrial.log2CP50M.genelevel.afterBatchEffectRemoval.eset)

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

## 4. prepare the master table
pd <- pData(eset)
master <- merge(pd, hdac6score, by = "row.names", all = T); dim(master)
head(master)

## 5. ROC curve fitting
roc1 <- roc(master$Response2, master$HDAC6regulon, plot=TRUE, legacy.axes = T, #percent = T,
            xlab = "1 - Specificity", ylab = "Sensitivity", col="red", lwd=4, print.auc=T, ci = T)
coords(roc1, "all", ret=c("threshold", "sens", "spec", "accuracy")) # accuracy, precision, recall


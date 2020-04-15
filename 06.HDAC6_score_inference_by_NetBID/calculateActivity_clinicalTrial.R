library(NetBID2)
load("/Users/qpan/Desktop/HDAC6Manuscript/clinicalTrial/Final/ClinicalTrial.log2TP50M.gene.50372_10.afterBatchEffectRemoval.eset")
eset <- ClinicalTrial.log2TP50M.gene.50372_10.afterBatchEffectRemoval.eset; rm(ClinicalTrial.log2TP50M.gene.50372_10.afterBatchEffectRemoval.eset)

######### 1 calculate HDAC6 score
## prepare regulon
regulon.file <- read.table("/Users/qpan/Desktop/HDAC6Manuscript/clinicalTrial/15_compare_regulon_new.txt", header = F, sep = "\t", stringsAsFactors = F)
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
head(master)

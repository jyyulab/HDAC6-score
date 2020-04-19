## To calculate the HDAC6 score of IBC
## Coded by Qingfei Pan (Qingfei.Pan@stjude.org)
## R-3.6

## 0. configuration
library(NetBID2)
setwd("./DATA/IBC")

## 1. load the eset
load("GEP.195.quantileNormalized.log2.eset")
eset <- GEP.195.quantileNormalized.log2.eset; rm(GEP.195.quantileNormalized.log2.eset)

## 2. read the regulon
regulon.file <- read.table("./HDAC6_Breast_Cancer_Regulon.txt", header = F, sep = "\t", stringsAsFactors = F)
regulon <- data.frame(row.names = regulon.file$V1, target = regulon.file$V1); head(regulon)
geneset <- list(); geneset[[1]] <- regulon; names(geneset) <- "HDAC6regulon"

## 3. calculate the HDAC6 score
fd <- fData(eset)[,c("ID", "Gene Symbol")];
colnames(fd) <- c("ID", "geneSymbol"); head(fd)
d <- exprs(eset)
exp <- merge(fd, d, by="row.names", all = T); exp <- exp[,-c(1:2)]
exp.gene <- aggregate(. ~ geneSymbol, exp, mean, na.rm = T)
exp.gene <- exp.gene[grepl("\\w", exp.gene$geneSymbol),]; dim(exp.gene)
row.names(exp.gene) <- exp.gene$geneSymbol
exp <- exp.gene[,-1]


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
eset.act <- generate.eset(exp_mat = d, phenotype_info = pd, feature_info = fd, annotation_info = 'IBC, GeneLevel, Activity')
IBC.activity.eset <- eset.act
save(IBC.activity.eset, file = "IBC.activity.eset")

require(dplyr)
require(openxlsx)

dir <- "/Volumes/yu3grp/solidTumor_JY/yu3grp/BRCA/HDAC6/Qingfei/02_subNetworks/01_TCGABRCA"
hubGenes_list <- "/Volumes/yu3grp/Network_JY/yu3grp/GeneAnnotation/TF_SIG/human/tf_sigs.updated.201806.xlsx"

## load the eset object
load(paste0(dir, "/00_inputEsets/TCGA.BRCA.RNASeq.isoform.950samples.20130326.RSEM.TP50M_normalized.eset"))
raw <- TCGA.BRCA.RNASeq.isoform.950samples.20130326.RSEM.TP50M_normalized.eset; rm(TCGA.BRCA.RNASeq.isoform.950samples.20130326.RSEM.TP50M_normalized.eset)

hubGenes <- read.xlsx(hubGenes_list, cols = c(1:2))
hubGenes_sig <- hubGenes[grepl("SIG", hubGenes$funcType),]
hubGenes_tf <- hubGenes[grepl("TF", hubGenes$funcType),]

## mark the status of ER, PR and HER2

d <- exprs(raw)
fd <- data.frame(isoformId = row.names(fData(raw)), geneSymbol = fData(raw)$geneSymbol)
fd <- fd[!grepl("\\?", fd$geneSymbol),] ## remove the rows with '?' as the gene symbol
fd <- fd[grepl("\\w", fd$geneSymbol),] ## remove the rows with 'BLANK' as the gene symbol
pd <- data.frame(sampleName = row.names(pData(raw)),
                 sampleType = pData(raw)$sampleType,
                 ER.IHC = pData(raw)$breastcarcinomaestrogenreceptorstatus,
                 PR.IHC = pData(raw)$breastcarcinomaprogesteronereceptorstatus,
                 HER2.IHC = pData(raw)$labprocher2neuimmunohistochemistryreceptorstatus)

## 1 clinical subgroup: Normal
pd_normal <- filter(pd, sampleType == '11') ## '11' represents NORMAL samples
dim(pd_normal)
d_normal <- select(data.frame(d), pd_normal$sampleName)
d_normal <- d_normal[rowSums(d_normal) > 0,] ## remove the rows with all '0'
dim(d_normal)
df_normal <- data.frame(row.names = row.names(d_normal), isoformId = row.names(d_normal))
df_normal <- merge(df_normal, fd, by = "isoformId", all = F); row.names(df_normal) <- df_normal$isoformId
df_normal <- merge(df_normal, d_normal, by = "row.names", all = F);
df_normal <- df_normal[, -1]

outdir <- paste0(dir, "/TCGABRCA_Normal"); dir.create(outdir, recursive = T)
outfile <- paste0(outdir, "/TCGABRCA.Normal_", length(unique(df_normal$isoformId)), "_", length(unique(df_normal$geneSymbol)), "_", dim(pd_normal)[1], ".exp")
write.table(df_normal, file = outfile, quote = F, row.names = F, col.names = T, sep = "\t")

sig <- df_normal[df_normal$geneSymbol%in%hubGenes_sig$geneSymbol,]
sig <- as.matrix(sig[,1])
outdir <- paste0(dir, "/TCGABRCA_Normal/sig"); dir.create(outdir, recursive = T)
outfile <- paste0(outdir, "/TCGABRCA.Normal_", length(unique(df_normal$isoformId)), "_", length(unique(df_normal$geneSymbol)), "_", dim(pd_normal)[1], "_sig.txt")
write.table(sig, file = outfile, quote = F, row.names = F, col.names = F, sep = "\t")

tf <- df_normal[df_normal$geneSymbol%in%hubGenes_tf$geneSymbol,]
tf <- as.matrix(tf[,1])
outdir <- paste0(dir, "/TCGABRCA_Normal/tf"); dir.create(outdir, recursive = T)
outfile <- paste0(outdir, "/TCGABRCA.Normal_", length(unique(df_normal$isoformId)), "_", length(unique(df_normal$geneSymbol)), "_", dim(pd_normal)[1], "_tf.txt")
write.table(tf, file = outfile, quote = F, row.names = F, col.names = F, sep = "\t")

## 2 clinical subgroup: HR
pd_HR <- filter(pd, sampleType == '01' & ER.IHC == "positive" & HER2.IHC == "negative") ## '01' represents Primary Solid Tumor
dim(pd_HR)
d_HR <- select(data.frame(d), pd_HR$sampleName)
d_HR <- d_HR[rowSums(d_HR) > 0,] ## remove the rows with all '0'
dim(d_HR)
df_HR <- data.frame(row.names = row.names(d_HR), isoformId = row.names(d_HR))
df_HR <- merge(df_HR, fd, by = "isoformId", all = F); row.names(df_HR) <- df_HR$isoformId
df_HR <- merge(df_HR, d_HR, by = "row.names", all = F);
df_HR <- df_HR[, -1]

outdir <- paste0(dir, "/TCGABRCA_HR"); dir.create(outdir, recursive = T)
outfile <- paste0(outdir, "/TCGABRCA.HR_", length(unique(df_HR$isoformId)), "_", length(unique(df_HR$geneSymbol)), "_", dim(pd_HR)[1], ".exp")
write.table(df_HR, file = outfile, quote = F, row.names = F, col.names = T, sep = "\t")

sig <- df_normal[df_normal$geneSymbol%in%hubGenes_sig$geneSymbol,]
sig <- as.matrix(sig[,1])
outdir <- paste0(dir, "/TCGABRCA_HR/sig"); dir.create(outdir, recursive = T)
outfile <- paste0(outdir, "/TCGABRCA.HR_", length(unique(df_normal$isoformId)), "_", length(unique(df_normal$geneSymbol)), "_", dim(pd_normal)[1], "_sig.txt")
write.table(sig, file = outfile, quote = F, row.names = F, col.names = F, sep = "\t")

tf <- df_normal[df_normal$geneSymbol%in%hubGenes_tf$geneSymbol,]
tf <- as.matrix(tf[,1])
outdir <- paste0(dir, "/TCGABRCA_HR/tf"); dir.create(outdir, recursive = T)
outfile <- paste0(outdir, "/TCGABRCA.HR_", length(unique(df_normal$isoformId)), "_", length(unique(df_normal$geneSymbol)), "_", dim(pd_normal)[1], "_tf.txt")
write.table(tf, file = outfile, quote = F, row.names = F, col.names = F, sep = "\t")

## 3 clinical subgroup: TN
pd_TN <- filter(pd, sampleType == '01' & ER.IHC == "negative" & HER2.IHC == "negative")
dim(pd_TN)
d_TN <- select(data.frame(d), pd_TN$sampleName)
d_TN <- d_TN[rowSums(d_TN) > 0,] ## remove the rows with all '0'
dim(d_TN)
df_TN <- data.frame(row.names = row.names(d_TN), isoformId = row.names(d_TN))
df_TN <- merge(df_TN, fd, by = "isoformId", all = F); row.names(df_TN) <- df_TN$isoformId
df_TN <- merge(df_TN, d_TN, by = "row.names", all = F);
df_TN <- df_TN[, -1]

outdir <- paste0(dir, "/TCGABRCA_TN"); dir.create(outdir, recursive = T)
outfile <- paste0(outdir, "/TCGABRCA.TN_", length(unique(df_TN$isoformId)), "_", length(unique(df_TN$geneSymbol)), "_", dim(pd_TN)[1], ".exp")
write.table(df_TN, file = outfile, quote = F, row.names = F, col.names = T, sep = "\t")

sig <- df_normal[df_normal$geneSymbol%in%hubGenes_sig$geneSymbol,]
sig <- as.matrix(sig[,1])
outdir <- paste0(dir, "/TCGABRCA_TN/sig"); dir.create(outdir, recursive = T)
outfile <- paste0(outdir, "/TCGABRCA.TN_", length(unique(df_normal$isoformId)), "_", length(unique(df_normal$geneSymbol)), "_", dim(pd_normal)[1], "_sig.txt")
write.table(sig, file = outfile, quote = F, row.names = F, col.names = F, sep = "\t")

tf <- df_normal[df_normal$geneSymbol%in%hubGenes_tf$geneSymbol,]
tf <- as.matrix(tf[,1])
outdir <- paste0(dir, "/TCGABRCA_TN/tf"); dir.create(outdir, recursive = T)
outfile <- paste0(outdir, "/TCGABRCA.TN_", length(unique(df_normal$isoformId)), "_", length(unique(df_normal$geneSymbol)), "_", dim(pd_normal)[1], "_tf.txt")
write.table(tf, file = outfile, quote = F, row.names = F, col.names = F, sep = "\t")

## 4 clinical subgroup: HER2
pd_HER2 <- filter(pd, sampleType == '01' & HER2.IHC == "positive")
dim(pd_HER2)
d_HER2 <- select(data.frame(d), pd_HER2$sampleName)
d_HER2 <- d_HER2[rowSums(d_HER2) > 0,] ## remove the rows with all '0'
dim(d_HER2)
df_HER2 <- data.frame(row.names = row.names(d_HER2), isoformId = row.names(d_HER2))
df_HER2 <- merge(df_HER2, fd, by = "isoformId", all = F); row.names(df_HER2) <- df_HER2$isoformId
df_HER2 <- merge(df_HER2, d_HER2, by = "row.names", all = F);
df_HER2 <- df_HER2[, -1]

outdir <- paste0(dir, "/TCGABRCA_HER2"); dir.create(outdir, recursive = T)
outfile <- paste0(outdir, "/TCGABRCA.HER2_", length(unique(df_HER2$isoformId)), "_", length(unique(df_HER2$geneSymbol)), "_", dim(pd_HER2)[1], ".exp")
write.table(df_HER2, file = outfile, quote = F, row.names = F, col.names = T, sep = "\t")

sig <- df_normal[df_normal$geneSymbol%in%hubGenes_sig$geneSymbol,]
sig <- as.matrix(sig[,1])
outdir <- paste0(dir, "/TCGABRCA_HER2/sig"); dir.create(outdir, recursive = T)
outfile <- paste0(outdir, "/TCGABRCA.HER2_", length(unique(df_normal$isoformId)), "_", length(unique(df_normal$geneSymbol)), "_", dim(pd_normal)[1], "_sig.txt")
write.table(sig, file = outfile, quote = F, row.names = F, col.names = F, sep = "\t")

tf <- df_normal[df_normal$geneSymbol%in%hubGenes_tf$geneSymbol,]
tf <- as.matrix(tf[,1])
outdir <- paste0(dir, "/TCGABRCA_HER2/tf"); dir.create(outdir, recursive = T)
outfile <- paste0(outdir, "/TCGABRCA.HER2_", length(unique(df_normal$isoformId)), "_", length(unique(df_normal$geneSymbol)), "_", dim(pd_normal)[1], "_tf.txt")
write.table(tf, file = outfile, quote = F, row.names = F, col.names = F, sep = "\t")

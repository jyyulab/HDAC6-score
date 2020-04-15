## To prepared the gene expression profiles of clinical trial samples
## Coded by Qingfei Pan (Qingfei.Pan@stjude.org)
## R-3.6


## 0. Configuration
library(NetBID2)
library(limma)
setwd("/Users/qpan/Desktop/HDAC6Manuscript/CodesAndData/DATA/ClinicalTrial")

## 2. Read the salmon outputs
samples <- c("pt03", "pt04", "pt05", "pt06", "pt07", "pt08", "pt09", "pt13", "pt14", "pt15rep1", "pt15rep2")
for (i in 1:length(samples)) {
  input_file <- paste0(samples[i], ".quantbySalmon.txt")
  input <- read.table(input_file, header = T)
  tmp <- input[, c(1,5)]; colnames(tmp) <- c("GeneID", samples[i])
  if (i == 1) { d <- tmp }
  else { d <- dplyr::full_join(d, tmp) }
}
row.names(d) <- d$GeneID; d <- d[,-1]
d <- d[rowSums(d) > 0,]

## normalization and transformation
d <- RNASeqCount.normalize.scale(mat = as.matrix(d), total = 1000000, pseudoCount = 0)
d <- log2(d + 1)

## merge the replicates of pt15
for (i in 1:nrow(d)) {
  d$pt15[i] <- mean(d$pt15rep1[i], d$pt15rep2[i])
}
d <- d[,c("pt03", "pt04", "pt05", "pt06", "pt07", "pt08", "pt09", "pt13", "pt14", "pt15")]

## prepare fd
fd <- data.frame(row.names = row.names(d), GeneSymbol = row.names(d), stringsAsFactors = F)

## prepare pd
pd <- data.frame(row.names = colnames(d), SampleID = colnames(d))
patients <- c("pt03", "pt04", "pt05", "pt06", "pt07", "pt08", "pt09", "pt13", "pt14", "pt15")
pd$Batch <- c("batch 1", "batch 1", "batch 1", "batch 1", "batch 1", "batch 1", "batch 3", "batch 2", "batch 2", "batch 3")
pd$Subtype <- c("HR+/HER2-", "HR+/HER2-", "TNBC", "TNBC", "TNBC", "HR+/HER2-", "HR+/HER2-", "HR+/HER2-", "HR+/HER2-", "HR+/HER2-")
pd$Response1 <- c("PD", "SD", "PD", "PD", "SD", "SD", "SD", "PR", "SD", "SD")
pd$Response2 <- c("Non-responder", "Responder", "Non-responder", "Non-responder", "Responder", "Responder", "Responder", "Responder", "Responder", "Responder")
pd$PFStime <- c(2.53, 5.03, 1.84, 1.08, 20.28, 11.41, 6.51, 4.47, 2.76, 5.19)
pd$PFSindicator <- c(1,0,1,1,1,1,1,1,0,1)

## create eset and qc
eset <- generate.eset(exp_mat = d, phenotype_info = pd[colnames(d), ], feature_info = fd[row.names(d), ], annotation_info = 'RNASeqRibo0, Salmon, GeneLevel, log2(CPM+1)')
dim(eset)

## remove non-informative genes
fData(eset)$IQR <- apply(exprs(eset), 1, IQR, na.rm = TRUE)  ## Calculate the IQR of each genes and save into fData
eset <- eset[order(-fData(eset)$IQR),] ## Sort the list by IQR in Decreasing order
eset <- eset[fData(eset)$IQR >= quantile(fData(eset)$IQR, 0.05)] ## Remove the genes with IQRs of bottom 5%


draw.eset.QC(eset, outdir = "./QC", do.logtransform = FALSE, prefix = 'Salmon.geneLevel_log2CPM.',
             intgroup = NULL, choose_plot = c("heatmap", "pca", "density", "correlation", "meansd"), generate_html = TRUE, correlation_strategy = "pearson", plot_all_point = FALSE,
             emb_plot_type='2D.text'
)

## batch effect correlation
exp.limma <- removeBatchEffect(exprs(eset), batch = pData(eset)$Batch)
eset.limma <- generate.eset(exp_mat = as.matrix(exp.limma), phenotype_info = pData(eset), feature_info = fData(eset), annotation_info = 'RNASeqRibo0, Salmon, GeneLevel, log2(CPM+1), BatchEffectCorrelatedBylimma')
draw.eset.QC(eset.limma, outdir = "./QC", do.logtransform = FALSE, prefix = 'Salmon.geneLevel_log2CPM.batchEffectCorrelatedBylimma',
             intgroup = NULL, choose_plot = c("heatmap", "pca", "density", "correlation", "meansd"), generate_html = TRUE, correlation_strategy = "pearson", plot_all_point = FALSE,
             emb_plot_type='2D.text' # "2D", "2D.interactive", "2D.ellipse", "2D.text" or "3D" 
)

## save eset
ClinicalTrial.log2CPM.genelevel.afterBatchEffectRemoval.eset <- eset.limma;
save(ClinicalTrial.log2CPM.genelevel.afterBatchEffectRemoval.eset, file = "./ClinicalTrial.log2CPM.genelevel.afterBatchEffectRemoval.eset")









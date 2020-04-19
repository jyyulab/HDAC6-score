## To prepare the gene expression profiles of cells from CCLE
## Coded by Qingfei Pan (Qingfei.Pan@stjude.org)
## R-3.6

## 0. Configuration
require(openxlsx)
require(Biobase)
require(stringr)
require(NetBID2)
setwd("./DATA/CCLE")

## 1. Prepare pd
pd <- read.xlsx("DepMap-2019q1-celllines.xlsx"); dim(pd) # Downloaded from: https://depmap.org/portal/download/api/download/external?file_name=processed_portal_downloads%2Fdepmap-public-cell-line-metadata-183e.4%2FDepMap-2019q1-celllines_v2.csv
pd <- pd[!grepl("^\\[MERGED", pd$CCLE_Name),]; dim(pd)
pd$CCLE_Name <- gsub("^(\\d.+)", "X\\1", pd$CCLE_Name)
row.names(pd) <- pd$DepMap_ID; dim(pd)

## 2. Prepare exp
exp <- read.csv("CCLE_depMap_19Q1_TPM.csv", header = T) # Downloaded from: https://depmap.org/portal/download/api/download/external?file_name=ccle%2Fdepmap-rnaseq-expression-data-ccd0.12%2FCCLE_depMap_19Q1_TPM.csv
row.names(exp) <- exp$X; exp <- exp[,-1]

table(row.names(pd) %in% row.names(exp)); table(row.names(exp) %in% row.names(pd))
exp.tmp <- merge(pd[,c(1,2)], exp, by = "row.names", all.y = T)
row.names(exp.tmp) <- exp.tmp$CCLE_Name; exp.tmp <- exp.tmp[,-c(1,2,3)]
exp <- t(exp.tmp)
row.names(exp) <- gsub("\\.\\..+$", "", row.names(exp))
exp <- exp[rowSums(exp)>0,]
exp <- exp[!duplicated(row.names(exp)),]

## 3. Prepare fd
fd <- data.frame(row.names = row.names(exp), GeneSymbol = row.names(exp))

## 4. Prepare eset
pd <- pd[pd$CCLE_Name %in% colnames(exp),];row.names(pd) <- pd$CCLE_Name
exp <- exp[row.names(fd),]; exp <- exp[, row.names(pd)] ## reorder the rows and columes of the expression matrix
CCLE_RNASeq.log2TPM_geneLevel.all.eset <- new("ExpressionSet", phenoData = new("AnnotatedDataFrame", pd), featureData=new("AnnotatedDataFrame", fd), annotation='Cell Line RNASeq', exprs=as.matrix(exp))
save(CCLE_RNASeq.log2TPM_geneLevel.all.eset, file = "CCLE_RNASeq.log2TPM_geneLevel.all.eset")


#draw.eset.QC(CCLE_RNASeq.log2TPM_geneLevel.all.eset, outdir = "./", do.logtransform = FALSE, prefix = 'CCLE_RNASeq.log2TPM_geneLevel.all.',
#             intgroup = NULL, choose_plot = c("heatmap", "pca", "density", "correlation", "meansd"), generate_html = TRUE, correlation_strategy = "pearson", plot_all_point = FALSE,
#             emb_plot_type='2D.ellipse' # "2D", "2D.interactive", "2D.ellipse", "2D.text" or "3D" 
#)

## 5. Extract the breast cancer cell lines
pd.breast <- pd[pd$Primary.Disease == "Breast Cancer", ]
exp.breast <- exp[, colnames(exp) %in% row.names(pd.breast)]
fd.breast <- data.frame(row.names = row.names(exp.breast), GeneSymbol = row.names(exp.breast))
exp.breast <- exp.breast[row.names(fd.breast),]; exp.breast <- exp.breast[, row.names(pd.breast)] ## reorder the rows and columes of the expression matrix
CCLE_RNASeq.log2TPM_geneLevel.breast.eset <- new("ExpressionSet", phenoData = new("AnnotatedDataFrame", pd.breast), featureData=new("AnnotatedDataFrame", fd.breast), annotation='Cell Line RNASeq', exprs=as.matrix(exp.breast))
save(CCLE_RNASeq.log2TPM_geneLevel.breast.eset, file = "CCLE_RNASeq.log2TPM_geneLevel.breast.eset")

draw.eset.QC(CCLE_RNASeq.log2TPM_geneLevel.breast.eset, outdir = "./", do.logtransform = FALSE, prefix = 'CCLE_RNASeq.log2TPM_geneLevel.breast.',
             intgroup = NULL, choose_plot = c("heatmap", "pca", "density", "correlation", "meansd"), generate_html = TRUE, correlation_strategy = "pearson", plot_all_point = FALSE,
             emb_plot_type='2D.ellipse' # "2D", "2D.interactive", "2D.ellipse", "2D.text" or "3D" 
)

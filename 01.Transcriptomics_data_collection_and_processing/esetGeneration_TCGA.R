## To prepared the gene expression profiles of TCGA samples
## Coded by Qingfei Pan (Qingfei.Pan@stjude.org)
## R-3.6

## 0. Configureation
require(openxlsx)
require(NetBID2)
setwd("/Users/qpan/Desktop/HDAC6Manuscript/CodesAndData/DATA/TCGA")


## 1. Load the TCGA dataset
load("TCGA.Oncoland.201503.FPKM.eset") # curated from Qiagen Oncoland TCGA dataset
eset <- TCGA.Oncoland.201503.FPKM.eset; rm(TCGA.Oncoland.201503.FPKM.eset)

## 2. Extract the BRCA eset
brca_rmpd <- eset[, pData(eset)$TumorType=="BRCA"]; dim(brca_rmpd) ## remove the non-BRCA samples
fd_new <- exprs(brca_rmpd)[rowSums(exprs(brca_rmpd)) > 0,]; dim(fd_new) ## remove the non-expressed
brca <- brca_rmpd[fData(brca_rmpd)$TranscriptIndex%in%row.names(fd_new),]

## 3. Read the subtype information of TCGA BRCA samples
subtype <- read.xlsx("13058_2016_724_MOESM2_ESM.xlsx", sheet = "RNA-Seq 1148") # Downloaded from: Additional file 2, PMID: 27386846 
colnames(subtype)[1] <- "newSampleID"

# Edit the phenotype data
pd_old <- pData(brca)
pd_old$newSampleID <- gsub("(.*)\\w", "\\1", pd_old$SampleID) ## Add the column of 'newSampleID' for merging
pd_new <- merge(subtype, pd_old, by = "newSampleID", all.x = F, all.y = T)

# Define the columns of 'SampleSubtype.Clinical' and 'SampleSubtype.Molecular'
for (i in 1:nrow(pd_new)) {
    if (pd_new$SampleType[i] == "Solid Tissue Normal") {
        pd_new$SampleSubtype.Clinical[i] <- "normalSamples"; pd_new$SampleSubtype.Molecular[i] <- "normalSamples"
    } else {
        if (is.na(pd_new$PAM50[i])) {
            pd_new$SampleSubtype.Molecular[i] <- "NA"
        } else {
            pd_new$SampleSubtype.Molecular[i] <- pd_new$PAM50[i]
        }
        
        if (is.na(pd_new$Her2.Status[i])) {
            pd_new$SampleSubtype.Clinical[i] <- "NA"
        } else if (pd_new$Her2.Status[i] == "Positive") {
            pd_new$SampleSubtype.Clinical[i] <- "HER2positive"
        } else if (pd_new$Her2.Status[i] == "Negative") {
            if (is.na(pd_new$ER.Status[i])) {
                pd_new$SampleSubtype.Clinical[i] <- "NA"
            } else if (pd_new$ER.Status[i] == "Negative") {
                pd_new$SampleSubtype.Clinical[i] <- "tripleNegative"
            } else if (pd_new$ER.Status[i] == "Positive") {
                pd_new$SampleSubtype.Clinical[i] <- "HRpositive"
            }
        }
    }
}

# Reorder the columns
cols.top <- c("SampleIndex", "SampleID", "newSampleID", "SampleType", "SampleSubtype.Clinical", "SampleSubtype.Molecular")
pd_new <- pd_new[, c(cols.top, colnames(pd_new)[!colnames(pd_new)%in%cols.top])]
row.names(pd_new) <- pd_new$SampleIndex
pData(brca) <- pd_new
exprs(brca) <- log2(exprs(brca) + 0.1)

## 4. Save the eset of isoform level
dim(brca)
TCGA_RNASeq.BRCA.log2FPKM_isoformLevel.73707_1221.eset <- brca
save(TCGA_RNASeq.BRCA.log2FPKM_isoformLevel.73707_1221.eset, file = "TCGA_RNASeq.BRCA.log2FPKM_isoformLevel.73707_1221.eset")

## 5. Create and Save the eset of gene level
pd <- pData(brca)[,c(1,3,5,6)]
pd$SampleID <- paste0("X", pd$SampleIndex)
row.names(pd) <- pd$SampleID
fd <- fData(brca)[,c(1,5)]

exp <- exprs(brca); colnames(exp) <- paste0("X", colnames(exp))
exp.tmp <- merge(fd, exp, by = "row.names", all.y = T)
exp.tmp <- exp.tmp[,-c(1,2)]
exp <- aggregate(. ~ GeneID, exp.tmp, max, na.rm = T) ## use the max isoform expression value to represent the gene's expression
row.names(exp) <- exp$GeneID; exp <- exp[, -1]

fd <- data.frame(row.names = row.names(exp), GeneSymbol = row.names(exp))

require(Biobase)
exp <- exp[row.names(fd),]; exp <- exp[, row.names(pd)] ## reorder the rows and columes of the expression matrix
TCGA_RNASeq.BRCA.log2FPKM_geneLevel.eset <- new("ExpressionSet", phenoData = new("AnnotatedDataFrame", pd), featureData=new("AnnotatedDataFrame", fd), annotation='TCGA BRCA', exprs=as.matrix(exp))
save(TCGA_RNASeq.BRCA.log2FPKM_geneLevel.eset, file = "TCGA_RNASeq.BRCA.log2FPKM_geneLevel.eset")

draw.eset.QC(TCGA_RNASeq.BRCA.log2FPKM_geneLevel.eset, outdir = dir.qc, do.logtransform = FALSE, prefix = 'TCGA_RNASeq.BRCA.log2FPKM_geneLevel.',
             intgroup = NULL, choose_plot = c("heatmap", "pca", "density", "correlation", "meansd"), generate_html = TRUE, correlation_strategy = "pearson", plot_all_point = FALSE,
             emb_plot_type='2D.ellipse' # "2D", "2D.interactive", "2D.ellipse", "2D.text" or "3D" 
)

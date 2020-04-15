## To prepared the gene expression profiles of IBC/non-IBC samples (PMID: 21339811)
## Coded by Qingfei Pan (Qingfei.Pan@stjude.org)
## R-3.6

## 0: Configuration
require(NetBID2)
require(limma)
setwd("/Users/qpan/Desktop/HDAC6Manuscript/CodesAndData/DATA/IBC")


## 1. read the expression data from GEO
eset.raw <- load.exp.GEO(out.dir = ".", GSE = 'GSE23720', GPL = 'GPL570', getGPL=TRUE)
exp.raw <- exprs(eset.raw)
fd.raw <- fData(eset.raw)

## 2. update the pd
pd.raw <- read.csv("pdata.197.csv", header = T) # pdata was curated from the paper: PMID: 21339811 
row.names(pd.raw) <- pd.raw$geo_accession

## 3. create the eset
eset <- new("ExpressionSet", phenoData = new("AnnotatedDataFrame", pd.raw), featureData=new("AnnotatedDataFrame", fd.raw), annotation='MicroArray', exprs=as.matrix(exp.raw))
dim(eset)
GEP.197.quantileNormalized.log2.eset <- eset
save(GEP.197.quantileNormalized.log2.eset, file = "GEP.197.quantileNormalized.log2.eset")

draw.eset.QC(GEP.197.quantileNormalized.log2.eset, outdir = dir.qc, do.logtransform = FALSE, prefix = 'GEP.197.quantileNormalized.log2.',
             intgroup = NULL, choose_plot = c("heatmap", "pca", "density", "correlation", "meansd"), generate_html = TRUE, correlation_strategy = "pearson", plot_all_point = FALSE,
             emb_plot_type='2D.ellipse' # "2D", "2D.interactive", "2D.ellipse", "2D.text" or "3D" 
)

## 4. remove outliers
eset.sel <- eset[, !(pData(eset)$sampleName %in% c("T60", "T61"))]; dim(eset.sel)
GEP.195.quantileNormalized.log2.eset <- eset.sel
save(GEP.195.quantileNormalized.log2.eset, file = "GEP.195.quantileNormalized.log2.eset")

draw.eset.QC(GEP.195.quantileNormalized.log2.eset, outdir = dir.qc, do.logtransform = FALSE, prefix = 'GEP.195.quantileNormalized.log2.',
             intgroup = NULL, choose_plot = c("heatmap", "pca", "density", "correlation", "meansd"), generate_html = TRUE, correlation_strategy = "pearson", plot_all_point = FALSE,
             emb_plot_type='2D.ellipse' # "2D", "2D.interactive", "2D.ellipse", "2D.text" or "3D" 
)






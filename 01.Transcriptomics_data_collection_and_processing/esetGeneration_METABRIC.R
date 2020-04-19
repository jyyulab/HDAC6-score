## To prepare the gene expression profiles of METABRIC samples
## Coded by Qingfei Pan (Qingfei.Pan@stjude.org)
## R-3.6

## 0. Configureation
require(openxlsx)
require(NetBID2)
require(Biobase)
setwd("./DATA/METABRIC")

## 1. Input the Rdata
load("brca_metabric_source.Rdata") # Downloaded from cBioPortal

## 2. Prepare the exp
exp <- all_samples_expr;
exp[is.na(exp)] <- 0
fd <- data.frame(row.names = row.names(exp), GeneID = row.names(exp))

## 3. Prepare the pd
pd <- patient_data
pd <- pd[colnames(exp),]

# Define the columns of 'SampleSubtype.Clinical' and 'SampleSubtype.Molecular'
for (i in 1:nrow(pd)) {
  if (pd$THREEGENE[i] == "HER2+") {
    pd$SampleSubtype.Clinical[i] <- "HER2positive"
  } else if (pd$THREEGENE[i] == "ER+/HER2- Low Prolif") {
    pd$SampleSubtype.Clinical[i] <- "HRpositive_LP"
  } else if (pd$THREEGENE[i] == "ER+/HER2- High Prolif") {
    pd$SampleSubtype.Clinical[i] <- "HRpositive_HP"
  } else if (pd$THREEGENE[i] == "ER-/HER2-") {
    pd$SampleSubtype.Clinical[i] <- "tripleNegative"
  } else {
    pd$SampleSubtype.Clinical[i] <- "NA"
  }
  
  if (pd$CLAUDIN_SUBTYPE[i] == "claudin-low") {
    pd$SampleSubtype.Molecular[i] <- "Claudin_Low"
  } else if (pd$CLAUDIN_SUBTYPE[i] == "Basal") {
    pd$SampleSubtype.Molecular[i] <- "Basal"
  } else if (pd$CLAUDIN_SUBTYPE[i] == "Normal") {
    pd$SampleSubtype.Molecular[i] <- "Normal"
  } else if (pd$CLAUDIN_SUBTYPE[i] == "Her2") {
    pd$SampleSubtype.Molecular[i] <- "Her2"
  } else if (pd$CLAUDIN_SUBTYPE[i] == "LumA") {
    pd$SampleSubtype.Molecular[i] <- "LumA"
  } else if (pd$CLAUDIN_SUBTYPE[i] == "LumB") {
    pd$SampleSubtype.Molecular[i] <- "LumB"
  } else {
    pd$SampleSubtype.Molecular[i] <- "NA"
  }
}

pd$SampleID <- row.names(pd)

# Reorder the columns
cols.top <- c("SampleID", "SampleSubtype.Clinical", "SampleSubtype.Molecular")
pd <- pd[, c(cols.top, colnames(pd)[!colnames(pd)%in%cols.top])]

## 4. Prepare eset
eset <- new("ExpressionSet", phenoData = new("AnnotatedDataFrame", pd), featureData=new("AnnotatedDataFrame", fd), annotation='MicroArray', exprs=as.matrix(exp))
METABRIC_microArray.BRCA.log2Intensity_geneLevel.eset <- eset
save(METABRIC_microArray.BRCA.log2Intensity_geneLevel.eset, file = "METABRIC_microArray.BRCA.log2Intensity_geneLevel.eset")

draw.eset.QC(METABRIC_microArray.BRCA.log2Intensity_geneLevel.eset, outdir = "./", do.logtransform = FALSE, prefix = 'METABRIC_microArray.BRCA.log2Intensity_geneLevel.',
             intgroup = NULL, choose_plot = c("heatmap", "pca", "density", "correlation", "meansd"), generate_html = TRUE, correlation_strategy = "pearson", plot_all_point = FALSE,
             emb_plot_type='2D.ellipse' # "2D", "2D.interactive", "2D.ellipse", "2D.text" or "3D" 
)

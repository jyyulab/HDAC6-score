## To prepare the gene expression profiles of GEMM from PMID: 24220145
## Coded by Qingfei Pan (Qingfei.Pan@stjude.org)
## R-3.6

## 0. Configureation
require(openxlsx)
require(GEOquery)
require(Biobase)
require(NetBID2)
setwd("./DATA/GEMM")

## 1. Prepare the raw pd
# get the sample information
sampleInfo <- read.xlsx("13059_2013_3180_MOESM1_ESM.xlsx", sheet = "Supplemental Table S1") # Downloaded from PMID: 24220145
table(sampleInfo$Model); table(sampleInfo$Expression.Class.Name)
pd <- data.frame(row.names = sampleInfo$GEO.Number, GEO.Number = sampleInfo$GEO.Number, Model.Name = sampleInfo$Model, Class.Name = sampleInfo$Expression.Class.Name,
                 Array.Type = sampleInfo$Expression.Array.Type, Sample.Note = sampleInfo$Expression.Sample, Strain = sampleInfo$Strain, Transgene1 = sampleInfo$Transgene.1,
                 Transgene2 = sampleInfo$Transgene.2, Comment = sampleInfo$Comment)

## 2. Prepare the raw exp
# Download data from GEO
dir.create("GSE42640"); gse42640 <- getGEO("GSE42640", GSEMatrix =TRUE, getGPL=TRUE, destdir = "GSE42640")
dir.create("GSE3165"); gse3165 <- getGEO("GSE3165", GSEMatrix =TRUE, getGPL=TRUE, destdir = "GSE3165")
dir.create("GSE8516"); gse8516 <- getGEO("GSE8516", GSEMatrix =TRUE, getGPL=TRUE, destdir = "GSE8516")
dir.create("GSE9343"); gse9343 <- getGEO("GSE9343", GSEMatrix =TRUE, getGPL=TRUE, destdir = "GSE9343")
dir.create("GSE14457"); gse14457 <- getGEO("GSE14457", GSEMatrix =TRUE, getGPL=TRUE, destdir = "GSE14457")
dir.create("GSE15263"); gse15263 <- getGEO("GSE15263", GSEMatrix =TRUE, getGPL=TRUE, destdir = "GSE15263")
dir.create("GSE17916"); gse17916 <- getGEO("GSE17916", GSEMatrix =TRUE, getGPL=TRUE, destdir = "GSE17916")
dir.create("GSE27101"); gse27101 <- getGEO("GSE27101", GSEMatrix =TRUE, getGPL=TRUE, destdir = "GSE27101")
dir.create("GSE35722"); gse35722 <- getGEO("GSE35722", GSEMatrix =TRUE, getGPL=TRUE, destdir = "GSE35722")

# Extract the overlapped probes
probe <- fData(gse42640[[1]])$SEQUENCE[grepl("\\w",fData(gse42640[[1]])$SEQUENCE)]; length(unique(probe))
probe <- intersect(probe, fData(gse42640[[2]])$SEQUENCE[grepl("\\w",fData(gse42640[[2]])$SEQUENCE)]); length(unique(probe))
probe <- intersect(probe, fData(gse42640[[3]])$SEQUENCE[grepl("\\w",fData(gse42640[[3]])$SEQUENCE)]); length(unique(probe))
probe <- intersect(probe, fData(gse42640[[4]])$SEQUENCE[grepl("\\w",fData(gse42640[[4]])$SEQUENCE)]); length(unique(probe))
probe <- intersect(probe, fData(gse3165[[4]])$SEQUENCE[grepl("\\w",fData(gse3165[[4]])$SEQUENCE)]); length(unique(probe))
probe <- intersect(probe, fData(gse8516[[1]])$SEQUENCE[grepl("\\w",fData(gse8516[[1]])$SEQUENCE)]); length(unique(probe))
probe <- intersect(probe, fData(gse9343[[1]])$SEQUENCE[grepl("\\w",fData(gse9343[[1]])$SEQUENCE)]); length(unique(probe))
probe <- intersect(probe, fData(gse14457[[1]])$SEQUENCE[grepl("\\w",fData(gse14457[[1]])$SEQUENCE)]); length(unique(probe))
probe <- intersect(probe, fData(gse14457[[3]])$SEQUENCE[grepl("\\w",fData(gse14457[[3]])$SEQUENCE)]); length(unique(probe))
probe <- intersect(probe, fData(gse15263[[1]])$SEQUENCE[grepl("\\w",fData(gse15263[[1]])$SEQUENCE)]); length(unique(probe))
probe <- intersect(probe, fData(gse17916[[1]])$SEQUENCE[grepl("\\w",fData(gse17916[[1]])$SEQUENCE)]); length(unique(probe))
probe <- intersect(probe, fData(gse27101[[1]])$SEQUENCE[grepl("\\w",fData(gse27101[[1]])$SEQUENCE)]); length(unique(probe))
probe <- intersect(probe, fData(gse27101[[3]])$SEQUENCE[grepl("\\w",fData(gse27101[[3]])$SEQUENCE)]); length(unique(probe))
probe <- intersect(probe, fData(gse35722[[1]])$SEQUENCE[grepl("\\w",fData(gse35722[[1]])$SEQUENCE)]); length(unique(probe))
probe <- intersect(probe, fData(gse35722[[2]])$SEQUENCE[grepl("\\w",fData(gse35722[[2]])$SEQUENCE)]); length(unique(probe))
probe <- intersect(probe, fData(gse35722[[3]])$SEQUENCE[grepl("\\w",fData(gse35722[[3]])$SEQUENCE)]); length(unique(probe))

## 3. Prepare the raw fd
fd.tmp <- fData(gse42640[[1]])[!duplicated(fData(gse42640[[1]])$SEQUENCE),]
fd <- fd.tmp[fd.tmp$SEQUENCE %in% probe,]
row.names(fd) <- fd$Agilent_probe_name; fd <- fd[,-1]

## 4. Refine the exp
arrays <- c(gse42640[[1]], gse42640[[2]], gse42640[[3]], gse42640[[4]], gse3165[[4]], gse8516[[1]], gse9343[[1]], gse14457[[1]], gse14457[[3]], gse15263[[1]], gse17916[[1]],
            gse27101[[1]], gse27101[[3]], gse35722[[1]], gse35722[[2]], gse35722[[3]])

for (i in 1:length(arrays)) {
  exp.tmp <- exprs(arrays[i][[1]]); exp.tmp <- exp.tmp[, colnames(exp.tmp) %in% pd$GEO.Number]
  fd.tmp <- data.frame(ID = fData(arrays[i][[1]])$ID, Sequence = fData(arrays[i][[1]])$SEQUENCE)
  exp.tmp <- merge(fd.tmp, exp.tmp, by.x = "ID", by.y = "row.names", all.x = T)
  
  exp.tmp <- exp.tmp[exp.tmp$Sequence %in% fd$SEQUENCE, ]
  exp.tmp <- aggregate(exp.tmp, list(SEQUENCE = exp.tmp$Sequence), mean)
  row.names(exp.tmp) <- exp.tmp$SEQUENCE; exp.tmp <- exp.tmp[, -c(1,2,3)]
  if (i == 1) { exp <- exp.tmp }
  else {
    exp <- merge(exp, exp.tmp, by = "row.names");
    row.names(exp) <- exp$Row.names;
    exp <- exp[, -1]
  }
}

# exp here includes duplicated columes which were set for normalization
# Check the cross-platform normalization
pd.neu <- pd[pd$Model.Name == "MMTV Neu",]; head(pd.neu)
exp.neu <- exp[, colnames(exp) %in% c(pd.neu$GEO.Number, paste0(pd.neu$GEO.Number, ".x"), paste0(pd.neu$GEO.Number, ".y"))]
boxplot(exp.neu)
exp.neu[,c(1,11)]

pd.c3tag <- pd[pd$Model.Name == "C3 Tag",]
exp.c3tag <- exp[, colnames(exp) %in% c(pd.c3tag$GEO.Number, paste0(pd.c3tag$GEO.Number, ".x"), paste0(pd.c3tag$GEO.Number, ".y"))]; dim(exp.c3tag)
boxplot(exp.c3tag)
exp.c3tag[,c(1,9)]
# Note: The same sample from different platforms share the same expression ratio values, which means the data has been normalized.

colnames(exp) <- sub("(.+).x", "\\1", colnames(exp)) ## remove the '.x' from the column names
exp <- exp[, !grepl("\\.y$", colnames(exp))] ## remove the duplicated columns marked with '.y'

fd.tmp <- data.frame(ProbeID = fd$Agilent_probe_name, Sequence = fd$SEQUENCE)
exp <- merge(fd.tmp, exp, by.x = "Sequence", by.y = "row.names")
table(grepl("^A_", exp$ProbeID)) # check if the spike probe removed or not
row.names(exp) <- exp$ProbeID; exp <- exp[,-c(1,2)]

table(colnames(as.matrix(exp)) %in% row.names(pd)); table(row.names(pd) %in% colnames(as.matrix(exp)));
table(row.names(as.matrix(exp)) %in% row.names(fd)); table(row.names(fd) %in% row.names(as.matrix(exp)));

## 5. Refine the pd based on Table 2 of PMID: 24220145
pd$PredictedHumanCounterpart <- "NA"
for (i in 1:nrow(pd)) {
  if (is.na(pd$Class.Name[i])) {
    next;
  }
  else if (pd$Class.Name[i] == "Erbb2-likeEx") {
    pd$PredictedHumanCounterpart[i] <- "HER2-enriched"
  }
  else if (pd$Class.Name[i] == "MycEx") {
    pd$PredictedHumanCounterpart[i] <- "Basal-like and Luminal B"
  }
  else if (pd$Class.Name[i] == "MycEx") {
    pd$PredictedHumanCounterpart[i] <- "Basal-like and Luminal B"
  }
  else if (pd$Class.Name[i] == "NeuEx") {
    pd$PredictedHumanCounterpart[i] <- "Luminal A"
  }
  else if (pd$Class.Name[i] == "Normal-likeEx") {
    pd$PredictedHumanCounterpart[i] <- "Normal-like"
  }
  else if (pd$Class.Name[i] == "p53null-BasalEx") {
    pd$PredictedHumanCounterpart[i] <- "Basal-like"
  }
  else if (pd$Class.Name[i] == "Class14Ex") {
    pd$PredictedHumanCounterpart[i] <- "Normal-like"
  }
  else if (pd$Class.Name[i] == "C3TagEx") {
    pd$PredictedHumanCounterpart[i] <- "Basal-like"
  }
  else if (pd$Class.Name[i] == "Claudin-lowEx") {
    pd$PredictedHumanCounterpart[i] <- "Claudin-low"
  }
}

## 6. eset of probes
exp <- exp[row.names(fd),]; exp <- exp[, row.names(pd)] ## reorder the rows and columes of the expression matrix
table(is.na(exp))
GEMM_microArray.probeLevel.eset <- new("ExpressionSet", phenoData = new("AnnotatedDataFrame", pd), featureData = new("AnnotatedDataFrame", fd), annotation='Mouse Model MicroArray', exprs = as.matrix(exp))
save(GEMM_microArray.probeLevel.eset, file = "GEMM_microArray.probeLevel.eset")

draw.eset.QC(GEMM_microArray.probeLevel.eset, outdir = "./", do.logtransform = FALSE, prefix = 'GEMM_microArray.probeLevel.',
             intgroup = NULL, choose_plot = c("heatmap", "pca", "density", "correlation", "meansd"), generate_html = TRUE, correlation_strategy = "pearson", plot_all_point = FALSE,
             emb_plot_type='2D.ellipse' # "2D", "2D.interactive", "2D.ellipse", "2D.text" or "3D" 
)


## 7. eset of genes
eset <- GEMM_microArray.probeLevel.eset; rm(GEMM_microArray.probeLevel.eset)
fd.tmp <- data.frame(row.names = fd$Agilent_probe_name, geneSymbol = fd$`Gene Symbol`)
exp.gene <- merge(fd.tmp, exp, by = "row.names")
for (i in 1:nrow(exp.gene)) { exp.gene$Row.names[i] <- unlist(strsplit(as.character(exp.gene$geneSymbol[i]), "\\|"))[1] }
exp.gene.b <- aggregate(.~Row.names, exp.gene, mean, na.rm = T)

exp.gene <- exp.gene.b[, -2] ## remove the unnecessary columns
exp.gene <- exp.gene[grepl("\\w", exp.gene$Row.names), ] ## remove the rows with no gene symbol
row.names(exp.gene) <- exp.gene$Row.names; exp.gene <- exp.gene[, -1]
dim(exp.gene)
fd.gene <- data.frame(row.names = row.names(exp.gene), geneSymbol = row.names(exp.gene))
pd.gene <- pData(eset)[colnames(exp.gene),]

GEMM_microArray.geneLevel.eset <- new("ExpressionSet", phenoData = new("AnnotatedDataFrame", pd.gene), featureData=new("AnnotatedDataFrame", fd.gene), annotation='Mouse Model MicroArray', exprs=as.matrix(exp.gene))
save(GEMM_microArray.geneLevel.eset, file = "GEMM_microArray.geneLevel.eset")

draw.eset.QC(GEMM_microArray.geneLevel.eset, outdir = "./", do.logtransform = FALSE, prefix = 'GEMM_microArray.geneLevel.',
             intgroup = NULL, choose_plot = c("heatmap", "pca", "density", "correlation", "meansd"), generate_html = TRUE, correlation_strategy = "pearson", plot_all_point = FALSE,
             emb_plot_type='2D.ellipse' # "2D", "2D.interactive", "2D.ellipse", "2D.text" or "3D" 
)



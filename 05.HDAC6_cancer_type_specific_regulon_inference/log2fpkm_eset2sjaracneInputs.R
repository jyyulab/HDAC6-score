############################################################
# This script is used to prepare SJARACNe inputs from TCGA eset.
# Qingfei Pan, 20190801
#
# NOTE:
#   1) cancer types, of which the sample size is less than 20, were discarded;
#   2) samples, of which the PC1 is greater than 3x standard deviation, were defined as outliers and removed;
#   3) isoforms, of which the expression is in bottom 5% among over 90% samples, were defined as low expression isoforms and removed;
#   4) isoforms, of which the IQR is lower than than 50% percentile of total IQR, were defined as low expression variation isoforms and removed;
#   5) isoforms of hub/driver genes, of which the IQR is between 10% and 50% percentile of total IQR, were rescued to improve the hug/driver coverage of networks.
############################################################

outdir <- "/Volumes/yu3grp/Network_JY/yu3grp/generatedNetworks/TCGA/sjaracne.OL201905/log2fpkm"

##### 1. read the new TF_SIG list #####
hub <- read.table("/Volumes/yu3grp/Network_JY/yu3grp/GeneAnnotation/TF_SIG/v201903/human/HubGenes_2008TFs_9659SIGs.human_v201903.txt", header = T, sep = "\t", stringsAsFactors = F, quote = "", fill = T)
tf <- unique(hub[hub$funcType != "SIG",]$HGNCsymbol); length(tf)
sig <- unique(hub[hub$funcType != "TF",]$HGNCsymbol); length(sig)
tf_sig <- unique(c(tf, sig)); length(tf_sig)

##### 2. load the eset object #####
load("/Volumes/yu3grp/Network_JY/yu3grp/DataSets/TCGA/OL201905/TCGA.Oncoland_201905.isoformLevel_log2FPKM.lowMappingSampleExcluded_97621_11125.eset")
eset <- TCGA.Oncoland_201905.isoformLevel_log2FPKM.lowMappingSampleExcluded_97621_11125.eset; rm(TCGA.Oncoland_201905.isoformLevel_log2FPKM.lowMappingSampleExcluded_97621_11125.eset)
dim(eset)

##### 3. define the sample groups #####
d <- exprs(eset); pd <- pData(eset)[, 1:7]; fd <- fData(eset)[, 1:8]
## 3.1 add group info to the pd
table(pd$TumorType); table(pd$SampleType)
for (i in 1:nrow(pd)) {
    if (grepl("Normal", pd$SampleType[i])) {
        pd$Group[i] <- paste0(pd$TumorType[i], "_", "N") # 'N' for NORMAL samples
    } else if (grepl("Primary", pd$SampleType[i])) {
        pd$Group[i] <- paste0(pd$TumorType[i], "_", "T") # 'T' for PRIMARY TUMOR samples
    } else {
        pd$Group[i] <- paste0(pd$TumorType[i], "_", "M") # 'M' for METASTASIS or RECURRENT TUMOR samples
    }
}
table(pd$Group)
grps <- names(table(pd$Group)[table(pd$Group)>=20]) # groups with <20 samples were discarded
pd <- pd[pd$Group %in% grps,]

## 3.2 remove the version number of Ensembl transcript ids

##### 4. prepare the SJARACNe inputs #####
for (i in 1:length(grps)) {
    cat(i,":", grps[i], "is in process...\n")
    
    ## 4.1 remove the outliers of samples
    pd.sel <- pd[pd$Group == grps[i],];
    cat("\t", dim(pd.sel)[1], "raw samples found...\n")
    d.sel <- d[,colnames(d) %in% pd.sel$SampleIndex];
    
    # Outlier is defined as the samples of which the PC1 is greater than 3x standard deviations. From: https://www.biostars.org/p/281767/
    pca <- prcomp(t(d.sel))
    PC <- data.frame(row.names = row.names(pca$x), SampleIndex = row.names(pca$x), pca$x[, c("PC1", "PC2")])
    pc1_sd <- sd(PC$PC1)*sqrt((length(PC$PC1) - 1) / (length(PC$PC1)))
    pc1_mean <- mean(PC$PC1)
    PC$zScore <- sapply(PC$PC1, function(x){(x - pc1_mean) / pc1_sd})
    pd.sel <- merge(pd.sel, PC, by = "SampleIndex", all = T)
    pd.sel.pca <- pd.sel[abs(pd.sel$zScore) <= 3,]
    cat("\t", dim(pd.sel.pca)[1], "qualified samples left after outlier removal...\n")
    
    ## 4.2 remove the low informative genes
    d.sel.pca <- d.sel[, colnames(d.sel) %in% pd.sel.pca$SampleIndex];
    
    # remove the genes of low expression, of which the expression is in the bottom 5% among > 90% samples
    d.sel.pca <- d.sel.pca[apply(d.sel.pca <= quantile(d.sel.pca, probs = 0.05), 1, sum) <= ncol(d.sel.pca) * 0.90, ]
    cat("\t", dim(d.sel.pca)[1], "isoforms left after low-expression gene removal...\n")
    
    # remove the genes of low expression variation, of which the IQR is of bottom 40%
    fd.sel <- fd[fd$TranscriptIndex %in% row.names(d.sel.pca),]; dim(fd.sel)
    fd.sel$IQR <- apply(d.sel.pca, 1, IQR)
    th1 <- quantile(fd.sel$IQR, 0.5); th2 <- quantile(fd.sel$IQR, 0.1); # change the cutoffs accordingly.
    fd.sel.iqr1 <- fd.sel[(fd.sel$IQR >= th1), ]; dim(fd.sel.iqr1)
    
    # rescue the hub genes of which the IQR is greater that bottom 10%
    fd.sel.iqr2 <- fd.sel[(fd.sel$GeneName %in% tf_sig) & (fd.sel$IQR > th2) & (fd.sel$IQR < th1), ]; dim(fd.sel.iqr2)
    
    fd.sel.iqr <- rbind(fd.sel.iqr1, fd.sel.iqr2)
    cat("\t", dim(fd.sel.iqr)[1], "isoforms left after low expression variation gene removal...\n")
    
    ## 4.3 write the sjaracne input files
    d.sel.pca.iqr <- d.sel.pca[row.names(d.sel.pca) %in% fd.sel.iqr$TranscriptIndex,]
    cat("\t", dim(d.sel.pca.iqr)[1], "isoforms cross", dim(d.sel.pca.iqr)[2], "samples will be used to build the network...\n")
    output_exp <- data.frame(row.names = fd.sel.iqr$TranscriptIndex, isoformId = fd.sel.iqr$TranscriptID, geneSymbol = fd.sel.iqr$GeneName)
    output_exp <- merge(output_exp, d.sel.pca.iqr, by = "row.names", all = T); output_exp <- output_exp[, -1]
    
    # expression matrix
    n_isoform <- length(unique(output_exp$isoformId)); n_gene <- length(unique(output_exp$geneSymbol)); n_sample <- dim(output_exp)[2] - 2;
    tag <- paste0(grps[i], ".", n_isoform, "_", n_gene, "_", n_sample)
    dir.create(paste0(outdir, "/", tag, "/sig"), recursive = T); dir.create(paste0(outdir, "/", tag, "/tf"), recursive = T);
    write.table(output_exp, file = paste0(outdir, "/", tag, "/", tag, ".exp"), row.names = F, col.names = T, quote = F, sep = "\t")
    
    # tf list
    tf.sel <- output_exp[output_exp$geneSymbol %in% tf,]
    tf.sel.nIsoform <- length(unique(tf.sel$isoformId)); tf.sel.nGene <- length(unique(tf.sel$geneSymbol))
    tf.sel.tag <- paste0(grps[i], ".", tf.sel.nIsoform, "_", tf.sel.nGene, "_", n_sample, ".tf.txt")
    write.table(tf.sel[,1], file = paste0(outdir, "/", tag, "/tf/", tf.sel.tag), row.names = F, col.names = F, quote = F, sep = "\t")
    
    # sig list
    sig.sel <- output_exp[output_exp$geneSymbol %in% sig,]
    sig.sel.nIsoform <- length(unique(sig.sel$isoformId)); sig.sel.nGene <- length(unique(sig.sel$geneSymbol))
    sig.sel.tag <- paste0(grps[i], ".", sig.sel.nIsoform, "_", sig.sel.nGene, "_", n_sample, ".sig.txt")
    write.table(sig.sel[,1], file = paste0(outdir, "/", tag, "/sig/", sig.sel.tag), row.names = F, col.names = F, quote = F, sep = "\t")
}

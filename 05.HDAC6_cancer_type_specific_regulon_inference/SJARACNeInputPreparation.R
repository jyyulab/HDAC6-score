## To prepare the SJARACNe inputs for subtype-specific regulon from TCGA
## Coded by Qingfei Pan (Qingfei.Pan@stjude.org)
## R-3.6

## 0. Configuration
require(NetBID2)
require(dplyr)
require(openxlsx)
setwd("./DATA")

## HubGene List
hub_sig <- read.table("./TCGA/NetBID_isoformID_geneID.SIG.txt", header = T, sep = "\t")
hub_sig <- hub_sig$TranscriptID
hub_tf <- read.table("./TCGA/NetBID_isoformID_geneID.TF.txt", header = T, sep = "\t")
hub_tf <- hub_tf$TranscriptID

## log2FPKM
load("./TCGA/TCGA.Oncoland.201503.FPKM.eset")
eset <- TCGA.Oncoland.201503.FPKM.eset

# modify the pd
table(pData(eset)$TumorType); table(pData(eset)$SampleType)
for (i in 1:nrow(pData(eset))) {
    if (grepl("Normal", pData(eset)$SampleType[i])) {
        pData(eset)$Group[i] <- paste0(pData(eset)$TumorType[i], "_", "N") # 'N' for NORMAL samples
    } else if (grepl("Primary", pData(eset)$SampleType[i])) {
        pData(eset)$Group[i] <- paste0(pData(eset)$TumorType[i], "_", "T") # 'T' for PRIMARY TUMOR samples
    } else {
        pData(eset)$Group[i] <- paste0(pData(eset)$TumorType[i], "_", "M") # 'M' for METASTASIS or RECURRENT TUMOR samples
    }
}
table(pData(eset)$Group)

# modify the exp
exprs(eset) <- log2(exprs(eset) + 0.1)

# refine the eset
grps <- names(table(pData(eset)$Group)[table(pData(eset)$Group)>=20]) # groups with <20 samples were discarded
eset <- eset[, pData(eset)$Group %in% grps]

## clinical subtypes
setwd("./DATA/subtypeSpecific")
type <- unique(pData(eset)$Group)
for (i in 1:length(type)) {
    if (type[i] == "NA") {next}
    eset.sel <- eset[, pData(eset)$Group == type[i]]
    
    expMatrix <- exprs(eset.sel)
    fd <- fData(eset.sel); fd$IQR <- apply(expMatrix, 1, IQR) # Calculate IQR for the Isoforms/Genes
    
    IQR.cutoff <- 0.6; IQR.rescue <- 0.3
    fd.cutoff <- fd[fd$IQR > quantile(fd$IQR, IQR.cutoff),]
    fd.cutoff.sig <- intersect(fd.cutoff$TranscriptID, hub_sig)
    fd.cutoff.tf <- intersect(fd.cutoff$TranscriptID, hub_tf)
    
    fd.rescue <- fd[fd$IQR >= quantile(fd$IQR, IQR.rescue) & fd$IQR <= quantile(fd$IQR, IQR.cutoff),]
    fd.rescue.sig <- intersect(fd.rescue$TranscriptID, hub_sig)
    fd.rescue.tf <- intersect(fd.rescue$TranscriptID, hub_tf)
    
    fd.sel.sig <- union(fd.cutoff.sig, fd.rescue.sig)
    fd.sel.tf <- union(fd.cutoff.tf, fd.rescue.tf)
    n.sig.total <- length(hub_sig); n.sig.sel <- length(fd.sel.sig); r.sig = n.sig.sel/n.sig.total; n.sig.cutoff <- length(fd.cutoff.sig); n.sig.rescue <- n.sig.sel - n.sig.cutoff
    n.tf.total <- length(hub_tf); n.tf.sel <- length(fd.sel.tf); r.tf = n.tf.sel/n.tf.total; n.tf.cutoff <- length(fd.cutoff.tf); n.tf.rescue <- n.tf.sel - n.tf.cutoff
    cat("Hub\t#Total\t#Filtered\t%Filtered\t#Cutoff\t#Rescue\nSIG\t", n.sig.total, "\t", n.sig.sel, "\t", r.sig, "\t", n.sig.cutoff, "\t", n.sig.rescue, "\nTF\t", n.tf.total, "\t", n.tf.sel, "\t", r.tf, "\t", n.tf.cutoff, "\t", n.tf.rescue, "\n")
    
    # Print outputs
    fd.filtered.index <- union(fd.cutoff$TranscriptIndex, fd.rescue$TranscriptIndex);
    fd.filtered <- fd[fd$TranscriptIndex%in%fd.filtered.index,];
    expMatrix.filtered <- expMatrix[row.names(expMatrix)%in%fd.filtered.index,]
    expMatrix.output <- data.frame(isoformId = fd.filtered$TranscriptID, geneSymbol = fd.filtered$GeneID, expMatrix.filtered)
    
    dir1 <- paste0(type[i], "/log2FPKM_0.1"); dir.create(dir1, recursive = T) #####
    write.table(expMatrix.output, file = paste0(dir1, "/GeneExpression.txt"), quote = F, row.names = F, col.names = T, sep = "\t")
    
    dir2 <- paste0(dir1, "/sig"); dir.create(dir2, recursive = T)
    write.table(as.matrix(fd.sel.sig), file = paste0(dir2, "/HubGene_sig.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
    
    dir3 <- paste0(dir1, "/tf"); dir.create(dir3, recursive = T)
    write.table(as.matrix(fd.sel.tf), file = paste0(dir3, "/HubGene_tf.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
    
}

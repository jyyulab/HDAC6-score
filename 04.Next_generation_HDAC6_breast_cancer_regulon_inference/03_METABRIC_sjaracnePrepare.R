require(dplyr)
require(openxlsx)
require(Biobase)
source("/Volumes/yu3grp/scRNASeq/yu3grp/qpan/Scripts/include_R/Functions_Included_QP.R")

## HubGene List
hub_sig <- read.table("/Volumes/yu3grp/scRNASeq/yu3grp/qpan/Database/GeneSymbol_List/NetBID_isoformID_geneID.SIG.txt", header = T, sep = "\t")
hub_sig <- hub_sig$GeneID
hub_tf <- read.table("/Volumes/yu3grp/scRNASeq/yu3grp/qpan/Database/GeneSymbol_List/NetBID_isoformID_geneID.TF.txt", header = T, sep = "\t")
hub_tf <- hub_tf$GeneID

indir <- "/Volumes/yu3grp/solidTumor_JY/yu3grp/BRCA/HDAC6/Qingfei/02_subNetworks/02_METABRIC/00_inputEsets"
outdir <- "/Volumes/yu3grp/solidTumor_JY/yu3grp/BRCA/HDAC6/Qingfei/02_subNetworks/02_METABRIC/01_sjaracne"

## logIntensity
load(paste0(indir, "/METABRIC.BRCA.PM50_Normal.24368_140.logIntensity.eset"))
expMatrix <- exprs(rmeset) ##### No more scale and log-transformation
fd <- fData(rmeset); fd$IQR <- apply(expMatrix, 1, IQR) # Calculate IQR for the Isoforms/Genes

# Check the IQR cutoffs
IQR.cutoff <- 0.6; IQR.rescue <- 0.3
fd.cutoff <- fd[fd$IQR > quantile(fd$IQR, IQR.cutoff),]
fd.cutoff.sig <- intersect(fd.cutoff$GeneID, hub_sig)
fd.cutoff.tf <- intersect(fd.cutoff$GeneID, hub_tf)

fd.rescue <- fd[fd$IQR >= quantile(fd$IQR, IQR.rescue) & fd$IQR <= quantile(fd$IQR, IQR.cutoff),]
fd.rescue.sig <- intersect(fd.rescue$GeneID, hub_sig)
fd.rescue.tf <- intersect(fd.rescue$GeneID, hub_tf)

fd.sel.sig <- union(fd.cutoff.sig, fd.rescue.sig)
fd.sel.tf <- union(fd.cutoff.tf, fd.rescue.tf)
n.sig.total <- length(unique(hub_sig)); n.sig.sel <- length(fd.sel.sig); r.sig = n.sig.sel/n.sig.total; n.sig.cutoff <- length(fd.cutoff.sig); n.sig.rescue <- n.sig.sel - n.sig.cutoff
n.tf.total <- length(unique(hub_tf)); n.tf.sel <- length(fd.sel.tf); r.tf = n.tf.sel/n.tf.total; n.tf.cutoff <- length(fd.cutoff.tf); n.tf.rescue <- n.tf.sel - n.tf.cutoff
cat("Hub\t#Total\t#Filtered\t%Filtered\t#Cutoff\t#Rescue\nSIG\t", n.sig.total, "\t", n.sig.sel, "\t", r.sig, "\t", n.sig.cutoff, "\t", n.sig.rescue, "\nTF\t", n.tf.total, "\t", n.tf.sel, "\t", r.tf, "\t", n.tf.cutoff, "\t", n.tf.rescue)

# Print outputs
fd.filtered.GeneID <- union(fd.cutoff$GeneID, fd.rescue$GeneID);
fd.filtered <- fd[fd$GeneID%in%fd.filtered.GeneID,];
expMatrix.filtered <- expMatrix[row.names(expMatrix)%in%fd.filtered.GeneID,]
expMatrix.output <- data.frame(GeneId = fd.filtered$GeneID, geneSymbol = fd.filtered$GeneID, expMatrix.filtered)

dir1 <- paste0(outdir, "/PM50_Normal/logIntensity"); dir.create(dir1, recursive = T) #####
write.table(expMatrix.output, file = paste0(dir1, "/GeneExpression.txt"), quote = F, row.names = F, col.names = T, sep = "\t")

dir2 <- paste0(dir1, "/sig"); dir.create(dir2, recursive = T)
write.table(as.matrix(fd.sel.sig), file = paste0(dir2, "/HubGene_sig.txt"), quote = F, row.names = F, col.names = F, sep = "\t")

dir3 <- paste0(dir1, "/tf"); dir.create(dir3, recursive = T)
write.table(as.matrix(fd.sel.tf), file = paste0(dir3, "/HubGene_tf.txt"), quote = F, row.names = F, col.names = F, sep = "\t")

require(NetBID2)
require(limma)

############### 1. HDAC6 regulons ###############
net.dir <- "/Volumes/yu3grp/Network_JY/yu3grp/generatedNetworks/TCGA/sjaracne/log2FPKM_0.1"
samples <- list.dirs(net.dir, full.names = F, recursive = F)
samples <- samples[grepl(".T_", samples)] # Remove Normal samples

HDAC6regulon <- list()
for (i in 1:length(samples)) {
  # SIG
  sig <- read.table(paste0(net.dir, "/", samples[i], "/sig/sjaracne2_out_.final/consensus_network_ncol_sjaracne2_out_.txt"), header = T, sep = "\t", stringsAsFactors = F)
  sig.sel <- sig[sig$source.symbol == "HDAC6",]; c.sig <- unique(sig.sel$target.symbol)
  regulon.sig <- data.frame(row.names = c.sig, target = c.sig, stringsAsFactors = F)
  name.sig <- gsub(".T_\\d+_\\d+_\\d+", "", samples[i])
  HDAC6regulon[[name.sig]] <- regulon.sig
  }

hdac6 <- read.table("/Volumes/project_space/yu3grp/solidTumor_JY/yu3grp/BRCA/HDAC6/Qingfei/04_Datasets/06_finalFigures/15_compare_regulon_new.txt", header = F, stringsAsFactors = F)
hdac6.new <- data.frame(row.names = hdac6$V1, target = hdac6$V1)
HDAC6regulon$BRCA <- hdac6.new

##### regulon lists

for (i in 1:length(HDAC6regulon)) {
  name <- names(HDAC6regulon)[i]
  vec <- unique(HDAC6regulon[[i]]$target)
  length <- length(vec)
  mat <- matrix(c(name, length, vec), nrow = 1, byrow = T)
  write.table(mat, file = paste0("/Volumes/project_space/yu3grp/solidTumor_JY/yu3grp/BRCA/HDAC6/Qingfei/04_Datasets/06_finalFigures/20191116/02_HDAC6Score_IC50_DenpendencyScore/RegulonOverlap/", i, ".txt"),
              col.names = F, row.names = F, sep = "\t")
}

##### fisher exact test
output <- c()
for (i in 1:(length(HDAC6regulon)-1)) {
  for (k in (i + 1):length(HDAC6regulon)) {
    name1 <- names(HDAC6regulon)[i]; name2 <- names(HDAC6regulon)[k]
    regulon1 <- unique(HDAC6regulon[[i]]$target); size1 <- length(regulon1)
    regulon2 <- unique(HDAC6regulon[[k]]$target); size2 <- length(regulon2)
    overlap <- intersect(regulon1, regulon2); size3 <- length(overlap)
    size4 <- 19799
    fisher <- fisher.test(matrix(c((size4 - (size1 - size3)), (size1 - size3), (size2 - size3), size3), nrow = 2, byrow = T))
    pvalue <- fisher$p.value
    output <- c(output, name1, name2, size1, size2, size3, size4, pvalue)
  }
}

output.matrix <- matrix(output, ncol = 7, byrow = T)
colnames(output.matrix) <- c("Regulon1", "Regulon2", "SizeOfRegulon1", "SizeOfRegulon2", "SizeOfOverlap", "SizeOfBackground", "Pvalue")
write.table(output.matrix, file = "/Volumes/project_space/yu3grp/solidTumor_JY/yu3grp/BRCA/HDAC6/Qingfei/04_Datasets/06_finalFigures/20191116/02_HDAC6Score_IC50_DenpendencyScore/RegulonOverlap/01_regulonOverlap.txt",
            quote = F, row.names = F, col.names = T, sep = "\t")

##### visualization

library("ggcorrplot")
library("reshape2")

## full all
output <- c()
for (i in 1:length(HDAC6regulon)) {
  for (k in 1:length(HDAC6regulon)) {
    name1 <- names(HDAC6regulon)[i]; name2 <- names(HDAC6regulon)[k]
    regulon1 <- unique(HDAC6regulon[[i]]$target); size1 <- length(regulon1)
    regulon2 <- unique(HDAC6regulon[[k]]$target); size2 <- length(regulon2)
    overlap <- intersect(regulon1, regulon2); size3 <- length(overlap)
    size4 <- 19799
    fisher <- fisher.test(matrix(c((size4 - (size1 - size3)), (size1 - size3), (size2 - size3), size3), nrow = 2, byrow = T))
    pvalue <- fisher$p.value
    output <- c(output, name1, name2, size1, size2, size3, size4, pvalue)
  }
}

output.matrix <- matrix(output, ncol = 7, byrow = T)

## full triangel
output <- c()
for (i in 1:(length(HDAC6regulon)-1)) {
  for (k in (i+1):length(HDAC6regulon)) {
    name1 <- names(HDAC6regulon)[i]; name2 <- names(HDAC6regulon)[k]
    regulon1 <- unique(HDAC6regulon[[i]]$target); size1 <- length(regulon1)
    regulon2 <- unique(HDAC6regulon[[k]]$target); size2 <- length(regulon2)
    overlap <- intersect(regulon1, regulon2); size3 <- length(overlap)
    size4 <- 19799
    fisher <- fisher.test(matrix(c((size4 - (size1 - size3)), (size1 - size3), (size2 - size3), size3), nrow = 2, byrow = T))
    pvalue <- fisher$p.value
    output <- c(output, name1, name2, size1, size2, size3, size4, pvalue)
  }
}

output.matrix <- matrix(output, ncol = 7, byrow = T)

df <- data.frame(Regulon1 = output.matrix[,1], Regulon2 = output.matrix[,2], Overlap = output.matrix[,5], Pvalue = output.matrix[,7], stringsAsFactors = F)
df$Pvalue <- as.numeric(df$Pvalue);
df$Overlap <- as.numeric(df$Overlap)

for (i in 1:nrow(df)) {
  if (df$Pvalue[i] > 0.01) {
    df$Pvalue.new[i] <- (df$Pvalue[i] - 0.01)/0.99*0.01 + 0.01
  } else {
    df$Pvalue.new[i] <- df$Pvalue[i]
  }
}
for (m in 1:nrow(df)) {
  if (df$Pvalue[m] <= 0.01) {
    df$Significant[m] <- "<=0.01"
    } else {
      df$Significant[m] <- ">0.01"
    }
}

df$Regulon1[df$Regulon1 == "COADREAD"] <- "CRC"; df$Regulon2[df$Regulon2 == "COADREAD"] <- "CRC"

sc <- scale_colour_gradient2(low = "red", mid = "grey", high = "blue", midpoint = 0.01, space = "Lab", limits = c(0,0.02))
p <- ggplot(df, aes(x=Regulon1, y = Regulon2, size = Overlap, color = Pvalue)) + geom_point(alpha=0.7) + sc +
  theme(legend.position = c(0.9,0.9), legend.box = "horizontal",
        legend.title = element_text(face = "bold"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_line(colour = "grey", size = 1),
        panel.background = element_blank(), #axis.line = element_line(colour = "black", size = 1),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 12, face = "bold", hjust = 1, color = "black", angle = 45),
        axis.text.y = element_text(size = 12, face = "bold", hjust = 1, color = "black")) +
  scale_y_discrete(limits = rev(unique(df$Regulon2)))
p
ggsave(filename = "/Volumes/project_space/yu3grp/solidTumor_JY/yu3grp/BRCA/HDAC6/Qingfei/04_Datasets/06_finalFigures/20191116/02_HDAC6Score_IC50_DenpendencyScore/RegulonOverlap/01_regulonOverlap_bubble.pdf", p, width = 10, height = 10, units = "in", dpi = 320)

sc <- scale_fill_gradient2(low = "red", mid = "grey95", high = "blue", midpoint = 0.01, space = "Lab", limits = c(0,0.02))
p <- ggplot(df, aes(x=Regulon1, y = Regulon2, fill = Pvalue)) + geom_tile(alpha=0.7, color = "white", size = 1) + sc +
  theme(legend.position = c(0.9,0.9), legend.box = "horizontal",
        legend.title = element_text(face = "bold"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_line(colour = "grey", size = 1),
        panel.background = element_blank(), #axis.line = element_line(colour = "black", size = 1),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 12, face = "bold", hjust = 1, color = "black", angle = 45),
        axis.text.y = element_text(size = 12, face = "bold", hjust = 1, color = "black")) +
  scale_y_discrete(limits = rev(unique(df$Regulon2)))
p
ggsave(filename = "/Volumes/project_space/yu3grp/solidTumor_JY/yu3grp/BRCA/HDAC6/Qingfei/04_Datasets/06_finalFigures/20191116/02_HDAC6Score_IC50_DenpendencyScore/RegulonOverlap/01_regulonOverlap_heatmap.pdf", p, width = 10, height = 10, units = "in", dpi = 320)


### clustering
regulon.all <- c()
for (i in 1:length(HDAC6regulon)) {
  if (i == 1) {
    regulon.all <- HDAC6regulon[[i]]$target
    } else {
      regulon.all <- c(regulon.all, HDAC6regulon[[i]]$target)
      }
  }
regulon.final <- unique(regulon.all); length(regulon.final)

df.final <- data.frame(row.names = regulon.final, stringsAsFactors = F)
for (i in 1:length(HDAC6regulon)) {
  df.final$name <- regulon.final %in% HDAC6regulon[[i]]$target
  colnames(df.final)[i] <- names(HDAC6regulon)[i]
}

d <- dist(as.matrix(t(df.final)))
hc <- hclust(d)
plot(hc)

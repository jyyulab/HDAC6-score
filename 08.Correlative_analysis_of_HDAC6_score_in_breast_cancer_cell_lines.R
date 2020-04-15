library(tidyr)
library(ggplot2)

outdir = "/Volumes/yu3grp/solidTumor_JY/yu3grp/BRCA/HDAC6/Qingfei/04_Datasets/06_finalFigures"

##### breastCellline_scaled #####
load("/Volumes/yu3grp/solidTumor_JY/yu3grp/BRCA/HDAC6/Qingfei/04_Datasets/02_cellLine/02_Dataset2_byDepMap/02_RNASeq_log2TPM.cellLine.breast.49417_56.activity_scaled.gene.eset")
eset <- RNASeq_log2TPM.cellLine.breast.49417_56.activity_scaled.gene.eset; rm(RNASeq_log2TPM.cellLine.breast.49417_56.activity_scaled.gene.eset)

d <- t(exprs(eset))
pd <- pData(eset)[,c(4,7)]
pd$Primary.Disease <- gsub("fibroblast", "Fibroblast", pd$Primary.Disease)
d <- merge(pd, d, by = "row.names")
d <- d[!is.na(d$IC50_of_ACY1215),]
data.scaled <- gather(d, Network, activityScore, HDAC6.ALLNormal_count_sig_425_373:HDAC6.MERGED.METABRICpan_IHCHRpositiveHP.TCGA_IHCHRpositive_fpkm.sig_294_269, factor_key = TRUE)

data1.scaled <- data.scaled[data.scaled$Network == "HDAC6.MERGED.METABRICpan_IHCHRpositiveHP.TCGA_IHCHRpositive_fpkm.sig_294_269",]
data2.scaled <- data1.scaled[data1.scaled$Primary.Disease == "Breast Cancer",]
data2.scaled <- data2.scaled[, c(1,2,5)]
##data2.scaled$IC50_of_ACY1215 <- log2(data2.scaled$IC50_of_ACY1215)
data2.scaled$activityScore <- 2^(data2.scaled$activityScore)
data2.scaled$Row.names <- gsub("_BREAST", "", data2.scaled$Row.names)

## 14 samples
p1 <- ggplot(data2.scaled, aes(x = IC50_of_ACY1215, y = activityScore)) + ylim(0,5)
p1 <- p1 + geom_point(size = 3) + labs(x = "IC50 of ACY-1215", y = "HDAC6 Score") +
  theme(legend.position = "",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 12, face = "bold", colour = "black"),
        axis.text.x = element_text(size = 10, face = "bold", hjust = 1, color = "black"),
        axis.text.y = element_text(size = 10, face = "bold", hjust = 1, color = "black"))
p1 <- p1 + stat_smooth(method = 'lm', formula = y ~ log2(x), se = F, fullrange = T)
p1
summary(lm(formula = activityScore ~ log2(IC50_of_ACY1215), data = data2.scaled))
cor.test(data2.scaled$IC50_of_ACY1215, data2.scaled$activityScore, method = "spearman", alternative = "less", exact = T)

ggsave(filename = paste0(outdir, "/13_CELLLINEDataset2_HDAC6regulon_BRCA_14samples/04_CCLE1165_HDAC6regulon_Breast_bycellline_smooth_14samples_Rsquared0.20_spearmance0.51pvalue0.0320.pdf"), p1, width = 7, height = 7, units = "in", dpi = 320)
ggsave(filename = paste0(outdir, "/13_CELLLINEDataset2_HDAC6regulon_BRCA_14samples/04_CCLE1165_HDAC6regulon_Breast_bycellline_smooth_14samples_Rsquared0.20_spearmance0.51pvalue0.0320.png"), p1, width = 7, height = 7, units = "in", dpi = 320)

p1 <- ggplot(data2.scaled, aes(x = IC50_of_ACY1215, y = activityScore, label = data2.scaled$Row.names)) + ylim(0,5) + xlim(0.05,35)
p1 <- p1 + geom_point(size = 3) + labs(x = "IC50 of ACY-1215", y = "HDAC6 Score") +
  theme(legend.position = "",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 12, face = "bold", colour = "black"),
        axis.text.x = element_text(size = 10, face = "bold", hjust = 1, color = "black"),
        axis.text.y = element_text(size = 10, face = "bold", hjust = 1, color = "black"))
p1 <- p1 + stat_smooth(method = 'lm', formula = y ~ log2(x), se = F, fullrange = T)
p1 <- p1 + geom_text(hjust = 0.5, nudge_x = 0, nudge_y = 0.1, size = 3, fontface = "bold")
p1
summary(lm(formula = activityScore ~ IC50_of_ACY1215, data = data2.scaled))
cor.test(data2.scaled$IC50_of_ACY1215, data2.scaled$activityScore, method = "spearman", alternative = "less", exact = T)

ggsave(filename = paste0(outdir, "/13_CELLLINEDataset2_HDAC6regulon_BRCA_14samples/04_CCLE1165_HDAC6regulon_Breast_bycellline_smooth_14samples_Rsquared0.20_spearmance0.51pvalue0.0320_lm2.pdf"), p1, width = 7, height = 7, units = "in", dpi = 320, useDingbats = F)
ggsave(filename = paste0(outdir, "/13_CELLLINEDataset2_HDAC6regulon_BRCA_14samples/04_CCLE1165_HDAC6regulon_Breast_bycellline_smooth_14samples_Rsquared0.20_spearmance0.51pvalue0.0320_lm2.png"), p1, width = 7, height = 7, units = "in", dpi = 320)


## 10 samples
data3.scaled <- data2.scaled[-c(2,6,8,13),]
p2 <- ggplot(data3.scaled, aes(x = IC50_of_ACY1215, y = activityScore)) + ylim(0,3)
p2 <- p2 + geom_point(size = 3) + labs(x = "IC50 of ACY-1215", y = "HDAC6 Score") +
  theme(legend.position = "",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 12, face = "bold", colour = "black"),
        axis.text.x = element_text(size = 10, face = "bold", hjust = 1, color = "black"),
        axis.text.y = element_text(size = 10, face = "bold", hjust = 1, color = "black"))
p2 <- p2 + stat_smooth(method = 'lm', formula = y ~ log2(x), se = F, fullrange = T)
p2
summary(lm(formula = activityScore ~ log2(IC50_of_ACY1215), data = data3.scaled))_
cor.test(data3.scaled$IC50_of_ACY1215, data3.scaled$activityScore, method = "spearman", alternative = "less", exact = T)

ggsave(filename = paste0(outdir, "/13_CELLLINEDataset2_HDAC6regulon_BRCA_10samples/04_CCLE1165_HDAC6regulon_Breast_bycellline_smooth_10samples_Rsquared0.70_spearmance0.62pvalue0.0301.pdf"), p2, width = 10, height = 7, units = "in", dpi = 320)
ggsave(filename = paste0(outdir, "/13_CELLLINEDataset2_HDAC6regulon_BRCA_10samples/04_CCLE1165_HDAC6regulon_Breast_bycellline_smooth_10samples_Rsquared0.70_spearmance0.62pvalue0.0301.png"), p2, width = 10, height = 7, units = "in", dpi = 320)
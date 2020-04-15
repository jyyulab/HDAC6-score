require(NetBID2)
require(ggplot2)

## load eset
load("/Users/qpan/Desktop/HDAC6Manuscript/clinicalTrial/Final/FourSets.TCGA_METABRIC_IBC_CLINICAL.afterBatchEffectRemoval.14719_3332.eset")
eset <- FourSets.TCGA_METABRIC_IBC_CLINICAL.afterBatchEffectRemoval.14719_3332.eset; rm(FourSets.TCGA_METABRIC_IBC_CLINICAL.afterBatchEffectRemoval.14719_3332.eset)
exp <- exprs(eset)

## load the regulon
regulon.file <- read.table("/Volumes/yu3grp/solidTumor_JY/yu3grp/BRCA/HDAC6/Qingfei/04_Datasets/06_finalFigures/15_compare_regulon_new.txt", header = F, sep = "\t", stringsAsFactors = F)
regulon <- data.frame(row.names = regulon.file$V1, target = regulon.file$V1); head(regulon)
geneset <- list(); geneset[[1]] <- regulon; names(geneset) <- "HDAC6regulon"

## calculate the HDAC6 score
hdac6score <- cal.Activity(target_list = geneset, cal_mat = as.matrix(exp),
                           es.method = 'mean', # 'Weightedmean', 'mean', 'maxmean', 'absmean'
                           std = T, memory_constrain = F # if true, the calculation strategy will not use Matrix Cross Products, which is memory consuming.
)
hdac6score <- data.frame(t(hdac6score))

std <- function(x){
  tmp_mean <- mean(x, na.rm = T); tmp_sd <- sd(x, na.rm = T); (x - tmp_mean) / tmp_sd
}
hdac6score$HDAC6regulon <- std(hdac6score$HDAC6regulon)

## pd
pd <- pData(eset)
master <- merge(pd, hdac6score, by = "row.names", all = T); dim(master)
table(master$Subtype1); table(master$Subtype2);

#########################################################
## TCGA & METABRIC
tcgametabric <- master[master$Dataset == "TCGA" | master$Dataset == "Metabric", ]; dim(tcgametabric)
tcgametabric$Subtype1[tcgametabric$Subtype1 == "HRpositive_HP"] <- "HRpositive"
tcgametabric$Subtype1[tcgametabric$Subtype1 == "HRpositive_LP"] <- "HRpositive"
table(tcgametabric$Subtype1)

## 1) IHC
tcgametabric.ihc <- tcgametabric[tcgametabric$Subtype1 != "NA",]
models <- unique(tcgametabric.ihc$Subtype1)
median_d <- data.frame()
for (i in 1:length(models)) {
  tag <- as.character(models[i])
  median_d[i,1] <- tag
  data.tmp <- tcgametabric.ihc[tcgametabric.ihc$Subtype1 == tag, ]
  median_d[i,2] <- median(data.tmp$HDAC6regulon, na.rm = T)
}
median_sorted <- median_d[order(-median_d$V2), ]
model_order1 <- median_sorted$V1
p1 <- ggplot(tcgametabric.ihc, aes(x=Subtype1, y=HDAC6regulon))
p1 <- p1 + geom_boxplot(aes(fill=Subtype1), outlier.shape = NA) + labs(x = "", y = "HDAC6 Score") + ylim(-4,4) + 
  geom_jitter(shape=16, position=position_jitter(0.3), size = 0.6) +
  theme(legend.position = "",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = "black", size = 1),
        axis.title = element_text(size = 15, face = "bold", colour = "black"),
        axis.text.x = element_text(size = 15, face = "bold", hjust = 1, color = "black", angle = 45),
        axis.text.y = element_text(size = 15, face = "bold", hjust = 1, color = "black")) +
  scale_x_discrete(limits = model_order1, labels=c("normalSamples" = "Paired-Normal\n(n=111)", "HRpositive" = "HRpositive\n(n=1708)",
                                                   "HER2positive" = "HER2positive\n(n=296)", "tripleNegative" = "TNBC\n(n=429)")) 
p1

dir <- "/Users/qpan/Desktop/HDAC6Manuscript/clinicalTrial/Final"
ggsave(filename = paste0(dir, "/03_TCGAMETABRIC_Clinical.pdf"), p1, width = 7, height = 7, units = "in", useDingbats = F)
ggsave(filename = paste0(dir, "/03_TCGAMETABRIC_Clinical.png"), p1, width = 7, height = 7, units = "in", dpi = 320)

## 2) Molecular
tcgametabric.mol <- tcgametabric[tcgametabric$Subtype2 != "NA",]
models <- unique(tcgametabric.mol$Subtype2)
median_d <- data.frame()
for (i in 1:length(models)) {
  tag <- as.character(models[i])
  median_d[i,1] <- tag
  data.tmp <- tcgametabric.mol[tcgametabric.mol$Subtype2 == tag, ]
  median_d[i,2] <- median(data.tmp$HDAC6regulon, na.rm = T)
}
median_sorted <- median_d[order(-median_d$V2), ]
model_order1 <- median_sorted$V1
p2 <- ggplot(tcgametabric.mol, aes(x=Subtype2, y=HDAC6regulon))
p2 <- p2 + geom_boxplot(aes(fill=Subtype2), outlier.shape = NA) + labs(x = "", y = "HDAC6 Score") + ylim(-4,4) + 
  geom_jitter(shape=16, position=position_jitter(0.3), size = 0.6) +
  theme(legend.position = "",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = "black", size = 1),
        axis.title = element_text(size = 15, face = "bold", colour = "black"),
        axis.text.x = element_text(size = 15, face = "bold", hjust = 1, color = "black", angle = 45),
        axis.text.y = element_text(size = 15, face = "bold", hjust = 1, color = "black")) +
  scale_x_discrete(limits = model_order1, labels=c("normalSamples" = "Paired-Normal\n(n=111)", "Normal" = "Normal-like\n(n=177)",
                                                   "Claudin_Low" = "Claudin-low\n(n=199)", "Basal" = "Basal-like\n(n=385)",
                                                   "LumA" = "Luminal-A\n(n=1217)", "LumB" = "Luminal-B\n(n=664)", "Her2" = "Her2-enriched\n(n=298)")) 
p2

dir <- "/Users/qpan/Desktop/HDAC6Manuscript/clinicalTrial/Final/"
ggsave(filename = paste0(dir, "/03_TCGAMETABRIC_Molecular.pdf"), p2, width = 7, height = 7, units = "in", useDingbats = F)
ggsave(filename = paste0(dir, "/03_TCGAMETABRIC_Molecular.png"), p2, width = 7, height = 7, units = "in", dpi = 320)

#########################################################
## all 4 groups
for (i in 1:nrow(master)) {
  if (master$Subtype1[i] == "IBC") { master$Subtype2[i] <- "IBC" }
  else if (master$Subtype1[i] == "NIBC") { master$Subtype2[i] <- "NIBC" }
  else if (master$Subtype2[i] == "Non-responder") { master$Subtype1[i] <- "Non-responder" }
  else if (master$Subtype2[i] == "Responder") { master$Subtype1[i] <- "Responder" }
}
master$Subtype1[master$Subtype1 == "HRpositive_HP"] <- "HRpositive"
master$Subtype1[master$Subtype1 == "HRpositive_LP"] <- "HRpositive"

write.table(master, file = "/Users/qpan/Desktop/HDAC6Manuscript/clinicalTrial/Final/03_4datasets_HDAC6score.masterTable.txt",
            col.names = T, row.names = F, sep = "\t", quote = F)

# HRpositive, HER2positive, tripleNegative, normalSamples
# LumB, Her2, LumA, Basal, Normal, Claudin_Low, normalSamples
c1 <- master$HDAC6regulon[master$Subtype2 == "Basal"]
c2 <- master$HDAC6regulon[master$Subtype2 == "LumA"]
x <- t.test(c1,c2); x$p.value



###### all4datasets: clinical
master.ihc <- master[master$Subtype1 != "NA",]
models <- unique(master.ihc$Subtype1)
median_d <- data.frame()
for (i in 1:length(models)) {
  tag <- as.character(models[i])
  median_d[i,1] <- tag
  data.tmp <- master.ihc[master.ihc$Subtype1 == tag, ]
  median_d[i,2] <- median(data.tmp$HDAC6regulon, na.rm = T)
}
median_sorted <- median_d[order(-median_d$V2), ]
model_order1 <- median_sorted$V1
p3 <- ggplot(master.ihc, aes(x=Subtype1, y=HDAC6regulon))
p3 <- p3 + geom_boxplot(aes(fill=Subtype1), outlier.shape = NA) + labs(x = "", y = "HDAC6 Score") + ylim(-4,4) + 
  geom_jitter(shape=16, position=position_jitter(0.3), size = 0.6) +
  theme(legend.position = "",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = "black", size = 1),
        axis.title = element_text(size = 15, face = "bold", colour = "black"),
        axis.text.x = element_text(size = 15, face = "bold", hjust = 1, color = "black", angle = 45),
        axis.text.y = element_text(size = 15, face = "bold", hjust = 1, color = "black")) +
  scale_x_discrete(limits = model_order1, labels=c("normalSamples" = "Paired-Normal\n(n=111)", "HER2positive" = "HER2positive\n(n=296)",
                                                   "HRpositive" = "HRpositive\n(n=1708)", "IBC" = "IBC\n(n=63)", "NIBC" = "NIBC\n(n=134)",
                                                   "Non-responder" = "Non-responder\n(n=3)", "Responder" = "Responder\n(n=7)", "tripleNegative" = "TNBC\n(n=429)")) 
p3

dir <- "/Users/qpan/Desktop/HDAC6Manuscript/clinicalTrial/Final/"
ggsave(filename = paste0(dir, "/03_all4groups_IHC.pdf"), p3, width = 10, height = 7, units = "in", useDingbats = F)
ggsave(filename = paste0(dir, "/03_all4groups_IHC.png"), p3, width = 10, height = 7, units = "in", dpi = 320)


c1 <- master.ihc$HDAC6regulon[master.ihc$Subtype2 == "Non-responder"]
c2 <- master.ihc$HDAC6regulon[master.ihc$Subtype2 == "Responder"]
t.test(c1,c2)
wilcox.test(c1, c2)


###### all4datasets: molecular
master.mol <- master[master$Subtype2 != "NA",]
models <- unique(master.mol$Subtype2)
median_d <- data.frame()
for (i in 1:length(models)) {
  tag <- as.character(models[i])
  median_d[i,1] <- tag
  data.tmp <- master.mol[master.mol$Subtype2 == tag, ]
  median_d[i,2] <- median(data.tmp$HDAC6regulon, na.rm = T)
}
median_sorted <- median_d[order(-median_d$V2), ]
model_order1 <- median_sorted$V1
p4 <- ggplot(master.mol, aes(x=Subtype2, y=HDAC6regulon))
p4 <- p4 + geom_boxplot(aes(fill=Subtype2), outlier.shape = NA) + labs(x = "", y = "HDAC6 Score") + ylim(-4,4) + 
  geom_jitter(shape=16, position=position_jitter(0.3), size = 0.6) +
  theme(legend.position = "",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = "black", size = 1),
        axis.title = element_text(size = 15, face = "bold", colour = "black"),
        axis.text.x = element_text(size = 15, face = "bold", hjust = 1, color = "black", angle = 45),
        axis.text.y = element_text(size = 15, face = "bold", hjust = 1, color = "black")) +
  scale_x_discrete(limits = model_order1, labels=c("normalSamples" = "Paired-Normal\n(n=111)", "Normal" = "Normal-like\n(n=177)",
                                                   "Claudin_Low" = "Claudin-low\n(n=199)", "Basal" = "Basal-like\n(n=385)",
                                                   "LumA" = "Luminal-A\n(n=1217)", "LumB" = "Luminal-B\n(n=664)", "Her2" = "Her2-enriched\n(n=298)",
                                                   "IBC" = "IBC\n(n=63)", "NIBC" = "NIBC\n(n=134)", "Non-responder" = "Non-responder\n(n=3)", "Responder" = "Responder\n(n=7)")) 
p4

dir <- "/Users/qpan/Desktop/HDAC6Manuscript/clinicalTrial/Final/"
ggsave(filename = paste0(dir, "/03_all4groups_PAM50.pdf"), p4, width = 10, height = 7, units = "in", useDingbats = F)
ggsave(filename = paste0(dir, "/03_all4groups_PAM50.png"), p4, width = 10, height = 7, units = "in", dpi = 320)


c1 <- master.ihc$HDAC6regulon[master.ihc$Subtype2 == "Non-responder"]
c2 <- master.ihc$HDAC6regulon[master.ihc$Subtype2 == "Responder"]
t.test(c1,c2)
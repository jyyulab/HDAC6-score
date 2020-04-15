library(NetBID2)
load("/Users/qpan/Desktop/HDAC6Manuscript/clinicalTrial/Final/ClinicalTrial.log2TP50M.gene.50372_10.afterBatchEffectRemoval.eset")
eset <- ClinicalTrial.log2TP50M.gene.50372_10.afterBatchEffectRemoval.eset; rm(ClinicalTrial.log2TP50M.gene.50372_10.afterBatchEffectRemoval.eset)

load("ClinicalTrial.log2CPM.genelevel.afterBatchEffectRemoval.eset")
eset <- ClinicalTrial.log2CPM.genelevel.afterBatchEffectRemoval.eset

######### 1 calculate HDAC6 score
## prepare regulon
regulon.file <- read.table("/Users/qpan/Desktop/HDAC6Manuscript/clinicalTrial/Final/15_compare_regulon_new.txt", header = F, sep = "\t", stringsAsFactors = F)
regulon <- data.frame(row.names = regulon.file$V1, target = regulon.file$V1); head(regulon)
geneset <- list(); geneset[[1]] <- regulon; names(geneset) <- "HDAC6regulon"

## calculate the HDAC6 score
exp <- exprs(eset)
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
head(master)

## boxplot
p <- ggplot(master, aes(x=Response2, y=HDAC6regulon, label = Row.names))
p <- p + geom_boxplot() + geom_point(aes(x=Response2, y=HDAC6regulon, color = Subtype), shape=16, size = 5) +
  labs(x = "", y = "HDAC6 Score") +
  geom_text(nudge_x = 0.08, nudge_y = 0.05, fontface = "bold")+
  theme(legend.position = c(0.15,0.9),
        legend.title = element_text(colour = "black", size = 12, face = "bold"),
        legend.text = element_text(colour="black", size=10, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = "black", size = 1),
        axis.title = element_text(size = 15, face = "bold", colour = "black"),
        axis.text.x = element_text(size = 15, face = "bold", hjust = 0.5, color = "black"),
        axis.text.y = element_text(size = 12, face = "bold", hjust = 1, color = "black")) 
p

c1 = master$HDAC6regulon[master$Response.new == "Non-responder"]
c2 = master$HDAC6regulon[master$Response.new == "Responder"]
t.test(c1, c2)
wilcox.test(c1, c2)
ggsave("/Users/qpan/Desktop/HDAC6Manuscript/clinicalTrial/Final/02_boxplot.NonrespVSresp_ttest0.03764_ranktest0.01667.pdf", p, width = 7, height = 7, units = "in", useDingbats = F)


## ROC curve
library(pROC)
roc1 <- roc(master$Response2, master$HDAC6regulon, plot=TRUE, legacy.axes = T, #percent = T,
            xlab = "1 - Specificity", ylab = "Sensitivity", col="red", lwd=4, print.auc=T, ci = T)
coords(roc1, "all", ret=c("threshold", "sens", "spec", "accuracy")) # accuracy, precision, recall


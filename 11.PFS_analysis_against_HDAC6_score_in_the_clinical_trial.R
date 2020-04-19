## PFS analysis
## Coded by Qingfei Pan (Qingfei.Pan@stjude.org)
## R-3.6

## 0. configuration
require(NetBID2)
require(survminer)
require(survival)
require(ggplot2)
setwd("./DATA/ClinicalTrial")

## 1. load the eset
load("ClinicalTrial.log2CP50M.genelevel.afterBatchEffectRemoval.eset")
eset <- ClinicalTrial.log2CP50M.genelevel.afterBatchEffectRemoval.eset; rm(ClinicalTrial.log2CP50M.genelevel.afterBatchEffectRemoval.eset)

## 2. read the regulon
regulon.file <- read.table("./HDAC6_Breast_Cancer_Regulon.txt", header = F, sep = "\t", stringsAsFactors = F)
regulon <- data.frame(row.names = regulon.file$V1, target = regulon.file$V1); head(regulon)
geneset <- list(); geneset[[1]] <- regulon; names(geneset) <- "HDAC6regulon"

## 3. calculate the HDAC6 score
exp <- exprs(eset)
hdac6score <- cal.Activity(target_list = geneset, cal_mat = as.matrix(exp),
                           es.method = 'mean', # 'Weightedmean', 'mean', 'maxmean', 'absmean'
                           std = T, memory_constrain = F # if true, the calculation strategy will not use Matrix Cross Products, which is memory consuming.
)
hdac6score <- data.frame(t(hdac6score))

# z-normalize the raw HDAC6 score
std <- function(x){
  tmp_mean <- mean(x, na.rm = T); tmp_sd <- sd(x, na.rm = T); (x - tmp_mean) / tmp_sd
}
hdac6score$HDAC6regulon <- std(hdac6score$HDAC6regulon)

## 4. prepare the master table
pd <- pData(eset)
master <- merge(pd, hdac6score, by = "row.names", all = T); dim(master)
head(master)

## 5. boxplot
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

## 6. Survival curve
master$Status <- c("Low", "High", "Low", "Low", "High", "High", "High", "High","High", "High") # based on the ROC curve fitting

fit <- survfit(Surv(PFStime) ~ Status, data = master)
ggsurvplot(fit, data = master)

p <- ggsurvplot(
  fit, 
  data = master, 
  size = 1,                 # change line size
  palette = c("red", "blue"),# custom color palettes
  conf.int = T,          # Add confidence interval
  pval = TRUE,              # Add p-value
  pval.size = 5, pval.coord = c(0.1,0.1),
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  xlab = "Time (Days)", ylab = "Survival Probability",
  legend.title = "HDAC6 Score",
  legend.labs = c("High", "Low"),    # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_rect(fill = "white", color = "black", size = 1),
                  legend.title = element_text(colour="black", size=15, face="bold"),
                  legend.text = element_text(colour="black", size=15, face="bold"),
                  axis.title = element_text(size = 15, face = "bold", colour = "black"),
                  axis.text.x = element_text(size = 12, face = "bold", hjust = 0.5, color = "black", angle = 0),
                  axis.text.y = element_text(size = 12, face = "bold", hjust = 1, color = "black"))      # Change ggplot2 theme
  
)
p

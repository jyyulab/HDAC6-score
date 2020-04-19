## To visulize the HDAC6 score distribution among datasets
## Coded by Qingfei Pan (Qingfei.Pan@stjude.org)
## R-3.6

## 0. configuration
require(NetBID2)
require(ggplot2)
setwd("./DATA")

## 1. load activity eset
load("MergedHumanBRCA.activity.eset")
eset <- MergedHumanBRCA.activity.eset; rm(MergedHumanBRCA.activity.eset)

d <- exprs(eset)
pd <- pData(eset)
master <- merge(pd, t(d), by = "row.names", all = T); dim(master)

#########################################################
table(master$Subtype1)
for (i in 1:nrow(master)) {
  if (master$Subtype1[i] == "HRpositive_HP") {
    master$Subtype1[i] <- "HRpositive"
  } else if (master$Subtype1[i] == "HRpositive_LP") {
    master$Subtype1[i] <- "HRpositive"
  } else if (master$Subtype1[i] == "tripleNegative") {
    master$Subtype1[i] <- "TNBC"
  } else if (master$Subtype1[i] == "normalSamples") {
    master$Subtype1[i] <- "Paired-Normal"
  }
}


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
p <- ggplot(master.ihc, aes(x=Subtype1, y=HDAC6regulon))
p <- p + geom_boxplot(aes(fill=Subtype1), outlier.shape = NA) + labs(x = "", y = "HDAC6 Score") + ylim(-4,4) + 
  geom_jitter(shape=16, position=position_jitter(0.3), size = 0.6) +
  theme(legend.position = "",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = "black", size = 1),
        axis.title = element_text(size = 15, face = "bold", colour = "black"),
        axis.text.x = element_text(size = 15, face = "bold", hjust = 1, color = "black", angle = 45),
        axis.text.y = element_text(size = 15, face = "bold", hjust = 1, color = "black")) +
  scale_x_discrete(limits = model_order1, labels=c("Paired-Normal" = "Paired-Normal\n(n=111)", "HER2positive" = "HER2+/(HR+ & HR-)\n(n=296)",
                                                   "HRpositive" = "HR+/HER2-\n(n=1708)", "IBC" = "IBC\n(n=63)", "NIBC" = "NIBC\n(n=134)",
                                                   "Non-responder" = "Non-responder\n(n=3)", "Responder" = "Responder\n(n=7)", "tripleNegative" = "TNBC\n(n=429)")) 
p
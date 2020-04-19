## For the correlative analysis between HDAC6 score and IC50
## Coded by Qingfei Pan (Qingfei.Pan@stjude.org)
## R-3.6

## 0. configuration
library(NetBID2)
library(tidyr)
library(ggplot2)
setwd("./DATA/CCLE")

## 1. load the eset
load("CCLE_BRCA.activity.eset")
eset <- CCLE_BRCA.activity.eset; rm(CCLE_BRCA.activity.eset)

## 2. prepare the master table
d <- exprs(eset)
pd <- pData(eset)
master <- merge(pd, t(d), by = "row.names", all = T)

## 14 samples
p1 <- ggplot(master, aes(x = log2(IC50_of_ACY1215), y = HDAC6regulon)) + xlim(0.1,5) + ylim(-2,3)
p1 <- p1 + geom_point(size = 3) + labs(x = "IC50 of ACY-1215", y = "HDAC6 Score") +
  theme(legend.position = "",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 12, face = "bold", colour = "black"),
        axis.text.x = element_text(size = 10, face = "bold", hjust = 1, color = "black"),
        axis.text.y = element_text(size = 10, face = "bold", hjust = 1, color = "black"))
p1 <- p1 + stat_smooth(method = 'lm', formula = y ~ log2(x), se = F, fullrange = T)
p1
summary(lm(formula = HDAC6regulon ~ log2(IC50_of_ACY1215), data = master))
cor.test(master$IC50_of_ACY1215, master$HDAC6regulon, method = "spearman", exact = T)

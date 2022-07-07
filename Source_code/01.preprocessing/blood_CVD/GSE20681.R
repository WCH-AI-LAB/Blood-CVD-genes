rm(list=ls())
options(stringsAsFactors = F)
setwd("c:/data/CVD_gene_DC/")

library(data.table)
library(limma)
library(plyr)
library(WGCNA)

#GSE20681
dat1 <- fread("C:/download/back_up_220414/GSE20681_series_matrix.txt", fill = T)
dat1 <- data.frame(dat1, stringsAsFactors = F)
dat1 <- dat1[-nrow(dat1), ]
sample = unlist(dat1[42, ])
table(sample)
#
dat1 <- dat1[-(1:70), ]
colnames(dat1)[1] <- "ID"
#45015

#99#99
idx.c = which(sample == "disease state: Control (0)")
idx.d = which(sample == "disease state: Case (1)")

dat1 <- dat1[, c(1, idx.c, idx.d)]
write.table(dat1, "temp.txt", col.names = F, row.names = F, sep = "\t")

#45015#199
dat1 = fread("temp.txt", header = T, stringsAsFactors = F)
dat1 = data.frame(dat1)
#
dat1[, -1] <- normalizeQuantiles(as.matrix(dat1[, -1]))
boxplot(dat1[sample(1:nrow(dat1), 2000), -1])

#
source("c:/src/BI/GPL4133.r")

#high_variance
library(matrixStats)
sd = rowSds(as.matrix(dat1[, -1]))
cut_off = quantile(sd, 0.4)
dat1 = dat1[sd > cut_off, ]

#selecting a probe with maxmean value
rownames(dat1) = paste0("s", (1:nrow(dat1)))
dat2 = collapseRows(dat1[, -1], rowID = rownames(dat1), rowGroup = dat1[, 1])
datExpr = dat2$datETcollapsed
cn = paste0("CN_", c(1:length(idx.c)))
dx = paste0("CAD_", c(1:length(idx.d)))
colnames(datExpr) = c(cn, dx)
#10748
write.table(datExpr, "GSE20681_maxmean.txt", sep = "\t")

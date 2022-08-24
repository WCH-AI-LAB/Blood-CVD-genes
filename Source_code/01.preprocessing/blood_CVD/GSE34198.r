rm(list=ls()); options(stringsAsFactors = F)
setwd("c:/data/CVD_gene_DC/")
library(data.table); library(limma)
library(plyr); library(WGCNA)

#GSE34198
dat1 <- fread("c:/download/GSE34198_series_matrix.txt", fill = T)
dat1 <- dat1[-nrow(dat1), ]
dat1 <- data.frame(dat1, stringsAsFactors = F)
#
sample = t(dat1[c(44:50), -1])
table(sample[, 1])

#48#45
idx.c = which(sample[, 1] == "group: control")
idx.d = which(sample[, 1] == "group: AIM")
dat1 = dat1[, c(1, idx.c, idx.d)]

#
dat1 <- dat1[-(1:83), ]
write.table(dat1, "temp.txt", col.names = F, row.names = F, sep = "\t")

#41093#89
dat1 = fread("temp.txt")
dat1 = data.frame(dat1)
#boxplot(dat1[sample(1:nrow(dat1), 1000), -1])
min = min(dat1[, -1])
min = as.numeric(min)
dat1[, -1] = dat1[, -1]-min+1
dat1[, -1] = log2(dat1[, -1])
dat1[, -1] = normalizeQuantiles(dat1[, -1])
#boxplot(dat1[sample(1:nrow(dat1), 1000), -1])
#high_variance
library(matrixStats)
sd = rowSds(as.matrix(dat1[, -1]))
cut_off = quantile(sd, 0.4)
dat1 = dat1[sd > cut_off, ]
source("c:/src/BI/GPL6102.R")
#30936
dat1 = dat1[, -1]

#selecting a probe with maxmean value
rownames(dat1) = paste0("s", (1:nrow(dat1)))
dat2 = collapseRows(dat1[, -1], rowID = rownames(dat1), rowGroup = dat1[, 1])
datExpr = dat2$datETcollapsed
#17823
#write.table(datExpr, "GSE34198_maxmean.txt", sep = "\t")
#write.table(sample, "GSE34198_sample.txt", sep = "\t", row.names = F, col.names = F)

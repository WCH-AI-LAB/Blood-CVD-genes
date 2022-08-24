rm(list=ls()); options(stringsAsFactors = F)
setwd("c:/data/CVD_gene_DC/")
library(data.table); library(limma)
library(plyr); library(WGCNA)

#GSE90074
dat1 <- fread("c:/download/GSE90074_series_matrix.txt", fill = T)
dat1 <- dat1[-nrow(dat1), ]
dat1 <- data.frame(dat1, stringsAsFactors = F)

sample = t(dat1[c(55:59), ])
nrow(sample)
#
table(sample[, 1])
idx.c = which(sample[, 1] == "obstructive_cad: N")
idx.d = which(sample[, 1] == "obstructive_cad: Y")

#50#93
dat1 <- dat1[-(1:87), c(1, idx.c, idx.d)]
write.table(dat1, "temp.txt", col.names = F, row.names = F, sep = "\t")

#41093#89
dat1 <- fread("temp.txt")
dat1 = data.frame(dat1)
dat1[, -1] = 2^(dat1[, -1])
dat1[, -1] = log2(dat1[, -1] + 1)
dat1[, -1] = normalizeQuantiles(dat1[, -1])

#high_variance
library(matrixStats)
sd = rowSds(as.matrix(dat1[, -1]))
cut_off = quantile(sd, 0.4)
dat1 = dat1[sd > cut_off, ]

source("c:/src/BI/GPL6480.r")
#30936


#selecting a probe with maxmean value
rownames(dat1) = paste0("s", (1:nrow(dat1)))
dat2 = collapseRows(dat1[, -1], rowID = rownames(dat1), rowGroup = dat1[, 1])
datExpr = dat2$datETcollapsed
nrow(datExpr)
#12598
write.table(datExpr, "GSE90074_maxmean.txt", sep = "\t")
write.table(sample, "GSE90074_sample.txt", sep = "\t", row.names = F, col.names = F)

rm(list = ls())

library(data.table)
library(plyr);library(limma);library(WGCNA)
setwd("c:/data/CVD_gene_DC/")

#GSE59867
dat1 <- fread("c:/download/back_up_220414/GSE59867_series_matrix.txt", fill = T)
dat1 <- dat1[-nrow(dat1), ]
dat1 <- data.frame(dat1, stringsAsFactors = F)

sample  = t(dat1[c(36, 46, 47), ])
dat1 <- dat1[-(1:67), ]
colnames(dat1)[1] <- "ID"
#33297
#46#111
idx1 = grep("Control", sample[, 1])
idx2 = grep("Patient", sample[, 1])
idx2 = setdiff(idx2, idx1)
temp_idx = grep("sampling 1", sample[, 1])
idx2 = intersect(idx2, temp_idx)
idx2 = setdiff(idx2, idx1)
idx.c = idx1
idx.d = idx2
dat1 <- dat1[, c(1, idx1, idx2)]

#46#111
write.table(dat1, "temp.txt", col.names = F, row.names = F, sep = "\t")

#[HuGene-1_0-st] Affymetrix Human Gene 1.0 ST Array
dat1 = fread("temp.txt", header = T, stringsAsFactors = F)
dat1 = data.frame(dat1)
boxplot(dat1[sample(1:nrow(dat1), 2000), -1])

#
source("c:/src/BI/GPL6244.r")

#high_variance
library(matrixStats)
sd = rowSds(as.matrix(dat1[, -1]))
cut_off = quantile(sd, 0.4)
dat1 = dat1[sd > cut_off, ]
#15069

#selecting a probe with maxmean value
rownames(dat1) = paste0("s", (1:nrow(dat1)))
dat2 = collapseRows(dat1[, -1], rowID = rownames(dat1), rowGroup = dat1[, 1])
datExpr = dat2$datETcollapsed
nrow(datExpr)
cn = paste0("CAD_", c(1:length(idx.c)))
dx = paste0("ACS_", c(1:length(idx.d)))
colnames(datExpr) = c(cn, dx)
#9700
write.table(datExpr, "GSE59867_maxmean.txt", sep = "\t")

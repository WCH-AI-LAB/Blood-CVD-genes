rm(list = ls())
setwd("c:/data/CVD_gene_DC/")
library(data.table);library(limma)

#
data = "GSE34198"
storage = paste0(data, "_maxmean.txt")
dat1 = fread(storage)
dat1 = data.frame(dat1, stringsAsFactors = F)
rownames(dat1) = dat1[, 1]
dat1 = dat1[, -1]

#48#45
Disease = c(rep(0, 48), rep(1, 45))
design = model.matrix(~ Disease)
fit = lmFit(dat1, design)
fit = eBayes(fit, trend = TRUE)
result = topTable(fit, coef=2, n = nrow(dat1), adjust="fdr", p = 1)
result = result[rownames(dat1), ]

#
storage = paste0(data, "_deg.txt")
write.table(result, storage, sep = "\t")

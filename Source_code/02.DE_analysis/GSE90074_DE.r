rm(list = ls())
setwd("c:/data/CVD_gene_DC/")
library(data.table);library(limma)

#
data = "GSE90074"
storage = paste0(data, "_maxmean.txt")
dat1 = fread(storage)
dat1 = data.frame(dat1, stringsAsFactors = F)
rownames(dat1) = dat1[, 1]
dat1 = dat1[, -1]

#50#93
Disease = c(rep(0, 50), rep(1, 93))
design = model.matrix(~ Disease)
fit = lmFit(dat1, design)
fit = eBayes(fit, trend = TRUE)
result = topTable(fit, coef=2, n = nrow(dat1), adjust="fdr", p = 1)
result = result[rownames(dat1), ]

#
storage = paste0(data, "_deg.txt")
write.table(result, storage, sep = "\t")

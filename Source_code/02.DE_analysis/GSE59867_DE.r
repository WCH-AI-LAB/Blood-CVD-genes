rm(list = ls())
setwd("c:/data/CVD_gene_DC/")
library(data.table);library(limma)

#
data = "GSE59867"
storage = paste0(data, "_maxmean.txt")
dat1 = fread(storage)
dat1 = data.frame(dat1, stringsAsFactors = F)
rownames(dat1) = dat1[, 1]
dat1 = dat1[, -1]

#62#59
status1 = names(dat1)
locate1 = str_locate(status1, "_")[, 1]
status1 = substr(status1, start = 1, stop = locate1-1)
#
table(status1)
#
Disease = c(rep(0, 46), rep(1, 111))
design = model.matrix(~ Disease)
fit = lmFit(dat1, design)
fit = eBayes(fit, trend = TRUE)
result = topTable(fit, coef=2, n = nrow(dat1), adjust="fdr", p = 1)
result = result[rownames(dat1), ]

#
storage = paste0(data, "_deg.txt")
write.table(result, storage, sep = "\t")

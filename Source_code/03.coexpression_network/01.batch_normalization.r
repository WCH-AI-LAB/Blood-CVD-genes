rm(list=ls())
setwd("c:/data/CVD_gene_DC/")

library(data.table)
library(sva)
library(stringr)

#
dat1 = fread("GSE20681_maxmean.txt")
dat2 = fread("GSE59867_maxmean.txt")
dat3 = fread("GSE60993_maxmean.txt")

#
dat1 = data.frame(dat1, stringsAsFactors = F)
dat2 = data.frame(dat2, stringsAsFactors = F)
dat3 = data.frame(dat3, stringsAsFactors = F)

#
gene_common = Reduce(intersect, x = list(dat1[, 1], dat2[, 1], dat3[, 1]))
idx1 = match(gene_common, dat1[, 1])
idx2 = match(gene_common, dat2[, 1])
idx3 = match(gene_common, dat3[, 1])
dat1 = dat1[idx1, -1]
dat2 = dat2[idx2, -1]
dat3 = dat3[idx3, -1]

#
dat_all = cbind(dat1, dat2, dat3)
pdf("batch.before.pdf", width = 9, height = 5)
boxplot(dat_all)
dev.off()

#
status1 = names(dat1)
status2 = names(dat2)
status3 = names(dat3)
locate1 = str_locate(status1, "_")[, 1]
locate2 = str_locate(status2, "_")[, 1]
locate3 = str_locate(status3, "_")[, 1]
status1 = substr(status1, start = 1, stop = locate1-1)
status2 = substr(status2, start = 1, stop = locate2-1)
status3 = substr(status3, start = 1, stop = locate3-1)

#
table(disease)
disease = c(status1, status2, status3)
mod = model.matrix(~disease)
batch = c(rep(1, ncol(dat1)), rep(2, ncol(dat2)), rep(3, ncol(dat3)))
batch = as.factor(batch)

#
datExpr = ComBat(dat_all, batch=batch, mod=mod, par.prior=TRUE)
pdf("batch.after.pdf", width = 9, height = 5)
boxplot(datExpr)
dev.off()

rownames(datExpr) = gene_common
write.table(datExpr, "data_integ_batch.txt", sep = "\t")

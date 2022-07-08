rm(list=ls())
setwd("c:/data/CVD_gene_DC/")

#
library(WGCNA);library(data.table)
dat = fread("data_integ_batch.txt")
dat = data.frame(dat, stringsAsFactors = F)
rownames(dat) = dat[, 1]
dat = dat[, -1]

#
powers = c(seq(1,14,by=1),seq(15,30,by=2))
sft = pickSoftThreshold(data= t(dat), 
                        networkType = "signed", 
                        corFnc="bicor",
                        verbose=5,
                        powerVector=powers)

#
pdf("FS.SoftThreshold.pdf",width=12,height=6)
par(mfrow=c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     xlab="Soft Thresh Power", ylab="Scale free R^2",type="n")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     labels = powers, cex = 0.7, col="red",  xlab="Soft Thresh Power", ylab="Scale free R^2")
abline(h=0.8, col="black")
plot(sft$fitIndices[,1], sft$fitIndices[,5], 
     xlab = "Soft threshold power", 
     ylab = "Mean connectivity", type = "n")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = 0.7, col="black")
dev.off()


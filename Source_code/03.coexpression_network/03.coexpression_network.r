rm(list=ls())
setwd("c:/data/CVD_gene_DC/")

#
library(WGCNA);library(data.table)
library(flashClust)
dat = fread("data_integ_batch.txt")
dat = data.frame(dat, stringsAsFactors = F)
rownames(dat) = dat[, 1]
dat = dat[, -1]

#SoftThreshold:15
adjacencyA1 = adjacency(t(dat),
                        power=15,
                        type="signed", 
                        corFnc = "bicor")
diag(adjacencyA1) = 0
dissTOMA1 = 1-TOMsimilarity(adjacencyA1, TOMType="signed")
geneTreeA1 = flashClust(as.dist(dissTOMA1), method="average")

save(file = "CoExpression_network_220203.rdata", dat, dissTOMA1, geneTreeA1)


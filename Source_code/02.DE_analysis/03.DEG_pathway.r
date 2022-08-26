rm(list = ls())
setwd("c:/data/CVD_gene_DC/")
library(data.table);library(limma)

#
GSE90074 = read.table("GSE90074_deg.txt")
GSE20681 = read.table("GSE20681_deg.txt")
GSE34198 = read.table("GSE34198_deg.txt")
GSE60993 = read.table("GSE60993_deg.txt")
gene_common = Reduce(intersect, x = list(rownames(GSE90074),
                                         rownames(GSE20681),
                                         rownames(GSE34198),
                                         rownames(GSE60993)))


#
GSE90074 = GSE90074[gene_common, ]
GSE20681 = GSE20681[gene_common, ]
GSE34198 = GSE34198[gene_common, ]
GSE60993 = GSE60993[gene_common, ]

#
GSE90074 = GSE90074[GSE90074$P.Value < 0.05, ]
GSE20681 = GSE20681[GSE20681$P.Value < 0.05, ]
GSE34198 = GSE34198[GSE34198$P.Value < 0.05, ]
GSE60993 = GSE60993[GSE60993$P.Value < 0.05, ]
n_data = c("GSE90074", "GSE20681", "GSE34198", "GSE60993")

#
candidate.gene = rownames(GSE90074)
module = "GSE90074_DEG"
background = 20000
p_cut_off = 0.05
source("c:/src/BI/pathway.r")

#
candidate.gene = rownames(GSE20681)
module = "GSE20681_DEG"
background = 18000
p_cut_off = 0.05
source("c:/src/BI/pathway.r")

#
candidate.gene = rownames(GSE34198)
module = "GSE34198_DEG"
background = 20000
p_cut_off = 0.05
source("c:/src/BI/pathway.r")

#
candidate.gene = rownames(GSE60993)
module = "GSE60993_DEG"
background = 15000
p_cut_off = 0.05
source("c:/src/BI/pathway.r")

#
path_gse90074 = read.table("GSE90074_DEG.pathway.csv", header = T, sep = ",")
path_gse20681 = read.table("GSE20681_DEG.pathway.csv", header = T, sep = ",")
path_gse34198 = read.table("GSE34198_DEG.pathway.csv", header = T, sep = ",")
path_gse60993 = read.table("GSE60993_DEG.pathway.csv", header = T, sep = ",")

#
path_common = Reduce(union, x = list(path_gse90074$ID,
                                     path_gse20681$ID,
                                     path_gse34198$ID,
                                     path_gse60993$ID))
path_common = unique(path_common)
res = data.frame(ID = path_common)
#
idx1 = match(res$ID, path_gse90074$ID)
path_gse90074 = path_gse90074[idx1, c("Num_pathway", "Num_common_genes")]
colnames(path_gse90074)[2] = "GSE90074"
res = cbind(res, path_gse90074)

#
idx1 = match(path_gse20681$ID, res$ID)
res$Num_pathway[idx1] = path_gse20681$Num_pathway
res[idx1, "GSE20681"] = path_gse20681$Num_common_genes

#
idx1 = match(path_gse34198$ID, res$ID)
res$Num_pathway[idx1] = path_gse34198$Num_pathway
res[idx1, "GSE34198"] = path_gse34198$Num_common_genes

#
idx1 = match(path_gse60993$ID, res$ID)
res$Num_pathway[idx1] = path_gse60993$Num_pathway
res[idx1, "GSE60993"] = path_gse60993$Num_common_genes

#
res$group = 0
idx1 = which(!is.na(res$GSE90074))
res$group[idx1] = res$group[idx1]+1
idx1 = which(!is.na(res$GSE20681))
res$group[idx1] = res$group[idx1]+1
idx1 = which(!is.na(res$GSE34198))
res$group[idx1] = res$group[idx1]+1
idx1 = which(!is.na(res$GSE60993))
res$group[idx1] = res$group[idx1]+1
write.table(res, "pathway_res.csv", sep = ",", row.names = F)

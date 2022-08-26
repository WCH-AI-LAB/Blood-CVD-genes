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
deg = list()
deg[[1]] = GSE90074[GSE90074$P.Value < 0.05, ]
deg[[2]] = GSE20681[GSE20681$P.Value < 0.05, ]
deg[[3]] = GSE34198[GSE34198$P.Value < 0.05, ]
deg[[4]] = GSE60993[GSE60993$P.Value < 0.05, ]
n_data = c("GSE90074", "GSE20681", "GSE34198", "GSE60993")

#4579#1179
gene_union = Reduce(union, x = list(rownames(GSE90074),
                                    rownames(GSE20681),
                                    rownames(GSE34198),
                                    rownames(GSE60993)))

#
mat = matrix("", ncol = 4, nrow = 4)
colnames(mat) = rownames(mat) = c("GSE90074", "GSE20681", "GSE34198", "GSE60993")

#
res = NULL
for(i in 1:nrow(mat)){
  for(j in 1:nrow(mat)){
    if(i <= j){
      gene1 = rownames(deg[[i]])
      gene2 = rownames(deg[[j]])
      n_common.genes = length(intersect(gene1, gene2))
      mat[n_data[i], n_data[j]] = n_common.genes  
    } else {
      gene1 = rownames(deg[[i]])
      gene2 = rownames(deg[[j]])
      n_common.genes = length(intersect(gene1, gene2))
      n_path.genes = length(gene1)
      n_candidate.genes = length(gene2)
      p = phyper(q=n_common.genes - 1, 
                 m=n_path.genes, 
                 n=6364 - n_path.genes, 
                 k=n_candidate.genes, lower.tail=F)
      p = format(p, digits = 3)
      mat[n_data[i], n_data[j]] = p  
      
    }
  }  
}

#
write.table(mat, "DEG_number.csv", sep = ",")


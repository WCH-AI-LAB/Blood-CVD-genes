rm(list = ls())
setwd("c:/data/CVD_gene_DC/")
library(data.table);library(limma)
library(ggplot2)

#
GSE90074 = read.table("GSE90074_deg.txt")
GSE20681 = read.table("GSE20681_deg.txt")
GSE34198 = read.table("GSE34198_deg.txt")
GSE60993 = read.table("GSE60993_deg.txt")
deg_list = list()
deg_list[[1]] = GSE90074
deg_list[[2]] = GSE20681
deg_list[[3]] = GSE34198
deg_list[[4]] = GSE60993
n_data = c("GSE90074", "GSE20681", "GSE34198", "GSE60993")

#
res = list()
for(i in 1:length(deg_list)){
  for(j in 1:length(deg_list)){
    k = 4*(i-1) + j
    print(k)
    if(i > j){
      fc.i = deg_list[[i]]
      fc.j = deg_list[[j]]
      #
      genes_common = intersect(rownames(fc.i), rownames(fc.j))
      fc.i = fc.i[genes_common, "logFC"]
      fc.j = fc.j[genes_common, "logFC"]
      dat.t = data.frame(x = fc.i, y = fc.j)
      res[[k]] = ggplot(dat.t, aes(x = x, y = y))+
        geom_point(color = "grey70")+
        geom_smooth(method = lm, se = T, linetype = "dashed", color = "red")+
        theme_bw()+
        xlab(n_data[i])+
        ylab(n_data[j])
    } else if(i == j){
      fc.i = deg_list[[i]]
      #
      res[[k]]= ggplot(fc.i, aes(x = logFC))+ 
        geom_histogram(color="black", fill="grey80", bins = 30)+
        theme_bw()+
        xlab(paste0(n_data[i], " (logFC)"))+
        ggtitle(paste0("Number of genes: ", nrow(fc.i)))
    } else {
      fc.i = deg_list[[i]]
      fc.j = deg_list[[j]]
      #
      genes_common = intersect(rownames(fc.i), rownames(fc.j))
      fc.i = fc.i[genes_common, "logFC"]
      fc.j = fc.j[genes_common, "logFC"]
      label = paste0("Number of common genes: ", length(genes_common))
      pcc = cor(fc.i, fc.j)
      scc = cor(fc.i, fc.j, method = "spearman")
      pcc = as.character(round(pcc, 3))
      scc = as.character(round(scc, 3))
      pcc = paste0("PCC: ", pcc)
      scc = paste0("SCC: ", scc)
      label = paste0(label, "\n", pcc, "\n", scc)
      res[[k]] = ggplot()+
        annotate("text", x = 5, y = 5, size = 4, label = label)+
        theme_void()
    }
  }
}

library(ggpubr)
ggarrange(plotlist = res, ncol = 4, nrow = 4)
ggsave("F.DE_FC.pdf", width = 12, height = 12)

rm(list=ls())
setwd("c:/data/CVD_gene_DC/")

#
library(WGCNA);library(data.table)
library(flashClust)

#
load("CoExpression_network_220203.rdata")
#DS=4,MMS=80,DCOR=0.1
tree = cutreeHybrid(
  dendro = geneTreeA1, 
  minClusterSize = 80, ############
  pamStage = FALSE, 
  cutHeight = 0.999, 
  deepSplit = 4, ################
  distM=dissTOMA1)
#
merged = mergeCloseModules(exprData= t(dat), 
                           colors = tree$labels, 
                           cutHeight=0.1   ############
)

#
colors = labels2colors(merged$colors)
entrez = rownames(dat)
save(file = "final_modules", colors, entrez, merged, tree)

#
GSE60993 = read.table("GSE60993_deg.txt")
GSE20681 = read.table("GSE20681_deg.txt")
GSE59867 = read.table("GSE59867_deg.txt")

#
GSE60993 = GSE60993[entrez, ]
GSE20681 = GSE20681[entrez, ]
GSE59867 = GSE59867[entrez, ]

#
GSE60993_c = rep("white", length(entrez))
GSE20681_c = rep("white", length(entrez))
GSE59867_c = rep("white", length(entrez))

#
idx.up1 = which(GSE60993$logFC > 0)
idx.up2 = which(GSE20681$logFC > 0)
idx.up3 = which(GSE59867$logFC > 0)
GSE60993_c[idx.up1] = "#ffcccc"
GSE20681_c[idx.up2] = "#ffcccc"
GSE59867_c[idx.up3] = "#ffcccc"

#
GSE60993_c[-idx.up1] = "#ccccff"
GSE20681_c[-idx.up2] = "#ccccff"
GSE59867_c[-idx.up3] = "#ccccff"

#
idx.up1 = which(GSE60993$logFC > 0 & GSE60993$adj.P.Val < 0.05)
idx.up2 = which(GSE20681$logFC > 0 & GSE20681$adj.P.Val < 0.05)
idx.up3 = which(GSE59867$logFC > 0 & GSE59867$adj.P.Val < 0.05)
GSE60993_c[idx.up1] = "#ff6666"
GSE20681_c[idx.up2] = "#ff6666"
GSE59867_c[idx.up3] = "#ff6666"

#
idx.down1 = which(GSE60993$logFC < 0 & GSE60993$adj.P.Val < 0.05)
idx.down2 = which(GSE20681$logFC < 0 & GSE20681$adj.P.Val < 0.05)
idx.down3 = which(GSE59867$logFC < 0 & GSE59867$adj.P.Val < 0.05)
GSE60993_c[idx.down1] = "#6666ff"
GSE20681_c[idx.down2] = "#6666ff"
GSE59867_c[idx.down3] = "#6666ff"

#fc
idx.up1 = which(GSE60993$logFC > 0.5 & GSE60993$adj.P.Val < 0.05)
idx.up2 = which(GSE20681$logFC > 0.5 & GSE20681$adj.P.Val < 0.05)
idx.up3 = which(GSE59867$logFC > 0.3 & GSE59867$adj.P.Val < 0.05)
GSE60993_c[idx.up1] = "red"
GSE20681_c[idx.up2] = "red"
GSE59867_c[idx.up3] = "red"

#
idx.down1 = which(GSE60993$logFC < 0 & GSE60993$adj.P.Val < 0.05)
idx.down2 = which(GSE20681$logFC < 0 & GSE20681$adj.P.Val < 0.05)
idx.down3 = which(GSE59867$logFC < -0.3 & GSE59867$adj.P.Val < 0.05)
GSE60993_c[idx.down1] = "blue"
GSE20681_c[idx.down2] = "blue"
GSE59867_c[idx.down3] = "blue"

#
pdf("F.module.pdf",width=8,height=5)
plotDendroAndColors(geneTreeA1,
                    cbind(colors, GSE60993_c, GSE20681_c, GSE59867_c),
                    groupLabels = c("Module", 
                                    "GSE60993\n(ACS/CN)", 
                                    "GSE20681\n(CAD/CN)", 
                                    "GSE59867\n(ACS/CAD)"),
                    cex.colorLabels = 1,
                    addGuide=T,
                    dendroLabels=F)
dev.off()

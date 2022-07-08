rm(list=ls())
setwd("c:/data/CVD_gene_DC/")

#
library(WGCNA);library(data.table)
library(flashClust)

#
load("CoExpression_network_220203.rdata")
#Iterate WGCNA parameters for robustness -- this takes a while
colors = vector(mode="list")
labels = vector(mode="list")

#
for(pam in c(FALSE))
  for (minModSize in c(50, 60, 70, 80, 90, 100)) {
    for (dthresh in c(0.1, 0.2)) {
      for(ds in c(0:4)) { 
        print(paste("DS=", ds, ",MMS=", minModSize, ",DCOR=",dthresh,",PAM=",pam,sep=""))
        
        tree = cutreeHybrid(dendro = geneTreeA1, 
                            minClusterSize= minModSize, 
                            pamStage = pam, 
                            cutHeight = 0.999, 
                            deepSplit = ds, 
                            distM = dissTOMA1)
        merged = mergeCloseModules(exprData= t(dat), colors = tree$labels, cutHeight=dthresh)
        colors = cbind(colors, labels2colors(merged$colors))
        
        labels = c(labels, paste("DS=", ds, ",MMS=", minModSize, ",DCOR=",dthresh,sep=""))
      }
    }
  }

#par(mar=c(do,lt,up,rt))
pdf("FS.WGCNA_diffParams_1.pdf",width=13,height=8)
plotDendroAndColors(geneTreeA1, colors[, 1:15], groupLabels = labels[1:15],
                    addGuide=T, dendroLabels=F, cex.colorLabels=0.6)
dev.off()

#
pdf("FS.WGCNA_diffParams_2.pdf",width=13,height=8)
plotDendroAndColors(geneTreeA1, colors[, 16:30], groupLabels = labels[16:30],
                    addGuide=T, dendroLabels=F, cex.colorLabels=0.6)
dev.off()

#
pdf("FS.WGCNA_diffParams_3.pdf",width=13,height=8)
plotDendroAndColors(geneTreeA1, colors[, 31:45], groupLabels = labels[31:45],
                    addGuide=T, dendroLabels=F, cex.colorLabels=0.6)
dev.off()

#
pdf("FS.WGCNA_diffParams_4.pdf",width=13,height=8)
plotDendroAndColors(geneTreeA1, colors[, 46:60], groupLabels = labels[46:60],
                    addGuide=T, dendroLabels=F, cex.colorLabels=0.6)
dev.off()

width = 1000
png("FS.WGCNA_diffParams_all.png",width=width,height=width*2)
plotDendroAndColors(geneTreeA1, colors, groupLabels = labels[1:15],
                    addGuide=T, dendroLabels=F, cex.colorLabels=0.6)
dev.off()

#DS=4,MMS=80,DCOR=0.1

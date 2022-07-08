rm(list=ls())
setwd("c:/data/CVD_gene_DC/")

#
library(WGCNA);library(data.table)
library(ggplot2)

#
load("final_modules.rdata")
dataset = c("GSE60993", "GSE20681", "GSE59867")
CTRL.list = c(7, 99, 46)
table(colors)
color.type = unique(colors)
color.type = color.type[color.type != "grey"]
moduleTraitP = matrix(1, nrow=length(color.type), ncol=length(dataset))
moduleTraitB = matrix(0, nrow=length(color.type), ncol=length(dataset))
colnames(moduleTraitP) = colnames(moduleTraitB) = dataset
rownames(moduleTraitP) = rownames(moduleTraitB) = color.type
moduleTraitSE = moduleTraitB
table(colors)
for(i in 1:length(dataset)){
  dat = fread(paste0(dataset[i], "_maxmean.txt"))
  dat = data.frame(dat, stringsAsFactors = F)
  rownames(dat) = dat[, 1]
  dat = dat[entrez, -1]
  MEs = moduleEigengenes(expr = t(dat), colors, softPower = 15)
  group = c(rep(0, CTRL.list[i]), rep(1, ncol(dat) - CTRL.list[i]))
  subject = colnames(dat)
  for(j in 1:length(color.type)){
    me_name = paste0("ME", color.type[j])
    me = MEs$eigengenes[[me_name]]
    dat.t = data.frame(me = me, group = group)
    mixedmodel = lm(me ~ group, data = dat.t)
    mixedmodel = summary(mixedmodel)$coefficients
    moduleTraitP[j, i] = mixedmodel["group", "Pr(>|t|)"]
    moduleTraitB[j, i] = mixedmodel["group", "Estimate"]
    moduleTraitSE[j, i] = mixedmodel["group", "Std. Error"]
  }
}

#
bpdata = melt(moduleTraitB)
semdata = melt(moduleTraitSE)
pdata = melt(moduleTraitP)

#
bpdata$sem = semdata$value
bpdata$p = pdata$value
bpdata$p.symbol = ""
bpdata$p.symbol[bpdata$p<0.1] = "#"
bpdata$p.symbol[bpdata$p<0.05] = "*"
bpdata$p.symbol[bpdata$p<0.01] = "**"
bpdata$p.symbol[bpdata$p<0.001] = "***"

colnames(bpdata)[1:2] = c("module", "dataset")
write.table(bpdata, "DE.txt", sep = "\t", row.names = F)

#
bpdata$module = factor(bpdata$module, levels = color.type)
bpdata$dataset = factor(bpdata$dataset, levels = dataset)

ggplot(bpdata, aes(x=dataset, y=value, fill=module, group=module, label=p.symbol)) + 
  facet_wrap(~ module, ncol=3, scales="free") + 
  geom_bar(stat="identity", position=position_dodge(), color="black") +   
  geom_errorbar(aes(ymin=(value - sem), ymax=(value + sem)), 
                position=position_dodge(.9), size=0.3,width=.3) +
  theme_minimal() + 
  scale_fill_manual(name="Group",values=levels(bpdata$module)) + 
  labs(y="beta", x="") + 
  geom_text(color="red",size=6,
            aes(y=value+ sign(value)*sem + sign(value)*.001), 
            position=position_dodge(.9))  + 
  scale_x_discrete() + 
  theme(
    legend.position = "none", 
    axis.text.y = element_text(size=10),
    axis.text.x = element_text(angle=45, hjust=1, size=10))

ggsave("FS.DE.PDF", height = 20, width = 9)

#
temp = bpdata[bpdata$p.symbol != "", ]
module_de = data.frame(table(temp$module))
module_de$Var1 = as.character(module_de$Var1)
module_de = module_de$Var1[module_de$Freq >= 2]
module_de = setdiff(module_de, c("lightgreen", "red", "royalblue", "black", "purple"))

#
bpdata = bpdata[bpdata$module %in% module_de, ]
bpdata$module = factor(bpdata$module, levels = module_de)
ggplot(bpdata, aes(x=dataset, y=value, fill=module, group=module, label=p.symbol)) + 
  facet_wrap(~ module, ncol=3, scales="free") + 
  geom_bar(stat="identity", position=position_dodge(), color="black") +   
  geom_errorbar(aes(ymin=(value - sem), ymax=(value + sem)), 
                position=position_dodge(.9), size=0.3,width=.3) +
  theme_minimal() + 
  scale_fill_manual(name="Group",values=levels(bpdata$module)) + 
  labs(y="beta", x="") + 
  geom_text(color="red",size=6,
            aes(y=value+ sign(value)*sem + sign(value)*.001), 
            position=position_dodge(.9))  + 
  scale_x_discrete() + 
  theme(
    legend.position = "none", 
    axis.text.y = element_text(size=10),
    axis.text.x = element_text(angle=45, hjust=1, size=10))

ggsave("F.DE.PDF", width = 9, height = 6)

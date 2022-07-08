rm(list=ls())
setwd("c:/data/CVD_gene_DC/")

#
library(WGCNA);library(data.table)
library(ggplot2)
library(riskyr)

#
load("final_modules.rdata")
dataset = c("GSE60993", "GSE20681", "GSE59867")
CTRL.list = c(7, 99, 46)
all_color = c("turquoise", "magenta", "yellow", "greenyellow", "tan")

df.all = NULL
i = 1
j = 1
for(i in 1:length(dataset)){
  dat = fread(paste0(dataset[i], "_maxmean.txt"))
  dat = data.frame(dat, stringsAsFactors = F)
  rownames(dat) = dat[, 1]
  dat = dat[entrez, -1]
  CTRL = CTRL.list[i]
  
  #
  dat_ctrl = dat[, 1:CTRL]
  dat_dx = dat[, (CTRL+1):ncol(dat)]
  
  #12
  adjMat_ctrl = adjacency(t(dat_ctrl), power=15, type="signed", corFnc = "bicor")
  adjMat_dx = adjacency(t(dat_dx), power=15, type="signed", corFnc = "bicor")
  
  #
  res_ctrl = intramodularConnectivity(adjMat_ctrl, colors, scaleByMax = FALSE)$kWithin
  res_dx = intramodularConnectivity(adjMat_dx, colors, scaleByMax = FALSE)$kWithin
  
  df = NULL
  for(j in 1:length(all_color)){
    con_ctrl = res_ctrl[colors == all_color[j]]
    con_dx = res_dx[colors == all_color[j]]
    a = t.test(con_ctrl, con_dx, paired = T)
    t_value = a$statistic
    p_value = a$p.value
    df.t = data.frame(module = all_color[j], t = t_value, p = p_value, 
                      m = mean(con_dx - con_ctrl), se = sd(con_dx - con_ctrl), n = length(con_ctrl))
    df = rbind(df, df.t)
  }
  df$data = dataset[i]
  df.all = rbind(df.all, df)
  print(i)
}


#
df = df.all
df$label = ""
df$label[df$p*15 < 0.05] = "*"
df$label[df$p*15 < 0.01] = "**"
df$label[df$p*15 < 0.001] = "***"
write.table(df, "DcoE.txt", sep = "\t", row.names = F)

#
df = read.table("DcoE.txt", header = T)
df$module = factor(df$module, levels = all_color)
df$data = factor(df$data, levels = dataset)
ggplot(df, aes(x= data, 
               ymin= (m - 2*se),
               lower= (m - se),
               middle = m,
               upper = (m + se), 
               ymax = (m + 2*se),
               fill=module, 
               label = label))+
  theme_bw()+
  geom_boxplot(position=position_dodge(width=0), stat="identity")+
  facet_grid(~module, scales="free_x") +
  geom_hline(yintercept=0, color = "red", size=0.5) +
  ylab("Z-score") + 
  geom_text(color = "red", size = 10, aes(y = (m + 2*se) + (m + 2*se)*.1), position = position_dodge(0.9))+
  scale_fill_manual(name="Group",values=levels(df$module))
  
ggsave("FS.DcoM.pdf", width =  12, height = 6)





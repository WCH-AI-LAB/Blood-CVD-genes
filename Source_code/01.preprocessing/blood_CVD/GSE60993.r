rm(list = ls())
setwd("c:/data/CVD_gene_DC/")
library(data.table); library(limma)
library(plyr); library(WGCNA)

#GSE60993
dat1 <- fread("c:/download/back_up_220414/GSE60993_series_matrix.txt", fill = T)
dat1 <- dat1[-nrow(dat1), ]
dat1 <- data.frame(dat1, stringsAsFactors = F)

sample  <- unlist(dat1[37, -1])
table(sample)
dat1 <- dat1[-(1:60), ]
#48803

probe_id <- dat1[, 1]
dat1 <- dat1[, -1]
dat1[] <- lapply(dat1, function(x) as.numeric(as.character(x)))
dat1 <- data.frame(dat1, stringsAsFactors = F)

#7#26
idx1 <- which(sample[] == "disease: Non-disease control")
idx2 <- which(sample[] == "disease: Non-ST-elevation MI patients"|
                sample[] == "disease: ST-elevation myocardial infarction patients"|
                sample[] == "disease: unstable angina")

dat1 <- dat1[, c(idx1, idx2)]
boxplot(dat1[sample(1:nrow(dat1), 5000), ])

dat1 <- cbind(probe_id, dat1)
write.table(dat1, "temp.txt", col.names = T, row.names = F, sep = "\t")

#48803#34
dat1 <- read.table("temp.txt", header = T, stringsAsFactors = F)
source("c:/src/BI/GPL6884.r")

#remove_probe_id
dat1 <- dat1[, -1]

#high_variance
library(matrixStats)
sd = rowSds(as.matrix(dat1[, -1]))
cut_off = quantile(sd, 0.4)
dat1 = dat1[sd > cut_off, ]

#selecting a probe with maxmean value
rownames(dat1) = paste0("s", (1:nrow(dat1)))
dat2 = collapseRows(dat1[, -1], rowID = rownames(dat1), rowGroup = dat1[, 1])
datExpr = dat2$datETcollapsed
cn = paste0("CN_", c(1:length(idx1)))
dx = paste0("ACS_", c(1:length(idx2)))
colnames(datExpr) = c(cn, dx)
#13670
#write.table(datExpr, "GSE60993_maxmean.txt", sep = "\t")


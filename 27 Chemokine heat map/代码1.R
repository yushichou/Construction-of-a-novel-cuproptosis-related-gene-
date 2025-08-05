rm(list = ls())
setwd("D:\\科研\\王晓 宫颈癌 铜死亡分型 结果+组图\\27 趋化因子热图")

library(RColorBrewer)
library(circlize)
library(gplots)
library(viridis)
library(oompaBase)
library(tidyverse)
library(pheatmap)

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)

# 数据标准化函数
standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {
  outdata = t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata > halfwidth] = halfwidth
    outdata[outdata < -halfwidth] = -halfwidth
  }
  return(outdata)
}

# 读取数据
load("merge.RDATA")
dataa = outTab
colnames(dataa) = str_replace_all(colnames(dataa), "TCGA_", "")
colnames(dataa) = str_replace_all(colnames(dataa), "GSE63514_", "")

aaa = read.csv("Cytokine.csv", header = T, sep = ",")
genename1 = intersect(aaa$ID, rownames(dataa))
dataa = dataa[genename1,]
aaa = dplyr::filter(aaa, ID %in% genename1)

risk = read.table("score.group.txt", header = T, sep = "\t", row.names = 1)
sameid = intersect(colnames(dataa), rownames(risk))
dataa = dataa[, sameid]
risk = risk[sameid,]

hmdat = as.data.frame(t(dataa))

# 注释
type = data.frame(Type = aaa$Type)
row.names(type) = aaa$ID
Typeid <- type$Type

annCol <- data.frame(score = scale(risk$score),
                     RiskType = risk$group,
                     row.names = rownames(risk),
                     stringsAsFactors = F)

# 基因差异分析（仅保留 q 值）
bbb = annCol[rownames(hmdat),]
sampleType = bbb$RiskType
pvals = c()

for(i in colnames(hmdat)){
  test = wilcox.test(hmdat[,i] ~ sampleType)
  pvals = c(pvals, test$p.value)
}

qvals = p.adjust(pvals, method = "fdr")

# 修改列名，仅加上 q 值
for(i in 1:length(colnames(hmdat))){
  genename = colnames(hmdat)[i]
  qval = qvals[i]
  qval_sci = formatC(qval, format = "e", digits = 2)  # 科学计数法，保留2位有效数字
  colnames(hmdat)[i] = paste0(genename, " (q=", qval_sci, ")")
}

# 行注释
annRow <- data.frame(Type = factor(Typeid, levels = unique(Typeid)),
                     row.names = colnames(hmdat),
                     stringsAsFactors = F)

# 注释颜色
Type.col <- brewer.pal(n = length(unique(Typeid)), name = "Paired")
annColors <- list(Type = c("Chemokines and receptors" = Type.col[1],
                           "Interieukines and receptors" = Type.col[2],
                           "Interferons and receptors" = Type.col[3],
                           "Other cytokines" = Type.col[4]),
                  "score" = greenred(64),
                  "RiskType" = c("High" = "red", "Low" = "blue"))

# 样本排序
samorder <- rownames(risk[order(risk$score),])

# 数据标准化
indata <- t(hmdat)
indata <- indata[, colSums(indata) > 0]
plotdata <- standarize.fun(indata, halfwidth = 2)

# 热图1
pheatmap::pheatmap(mat = as.matrix(plotdata[, samorder]),
                   border_color = NA,
                   color = bluered(64),
                   cluster_rows = F,
                   cluster_cols = F,
                   show_rownames = T,
                   show_colnames = F,
                   annotation_col = annCol[samorder,, drop = F],
                   annotation_row = annRow,
                   annotation_colors = annColors,
                   gaps_col = table(annCol$RiskType)[2],
                   gaps_row = cumsum(table(annRow$Type)),
                   cellwidth = 0.8,
                   cellheight = 10,
                   filename = "immune heatmap by pheatmap_score1_qval_only.pdf")

# 平均表达热图
ccc = hmdat[rownames(annCol)[annCol$RiskType=="High"],]
ddd = hmdat[rownames(annCol)[annCol$RiskType=="Low"],]
ccc1 = colMeans(ccc)
ddd1 = colMeans(ddd)
mewdata = data.frame(High = ccc1, Low = ddd1)

annCol1 = data.frame(RiskType = c("High", "Low"))
row.names(annCol1) = c("High", "Low")

pheatmap::pheatmap(mat = as.matrix(mewdata[, row.names(annCol1)]),
                   border_color = "grey",
                   color = bluered(64),
                   cluster_rows = F,
                   cluster_cols = F,
                   show_rownames = T,
                   show_colnames = F,
                   annotation_col = annCol1,
                   annotation_row = annRow,
                   annotation_colors = annColors,
                   gaps_row = cumsum(table(annRow$Type)),
                   cellwidth = 20,
                   cellheight = 10,
                   filename = "immune heatmap by pheatmap_score_mean1_qval_only.pdf")


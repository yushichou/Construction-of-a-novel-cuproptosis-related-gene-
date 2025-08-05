rm(list = ls())

#引用包
library(corrplot)
scoreFile="score.txt"         #打分文件
immFile="ssGSEA.result.txt"      #ssGSEA结果文件

#读取m6A打分文件
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
colnames(score)="score"
#读取免疫细胞文件
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
library(tidyverse)
rownames(immune)=str_replace_all(rownames(immune),"na","")
immune=t(immune)

#数据合并
sameSample=intersect(row.names(score), row.names(immune))
data=cbind(score[sameSample,,drop=F], immune[sameSample,,drop=F])

#相关性矩阵
M=cor(data)
res1=cor.mtest(data, conf.level = 0.95)

#绘制相关性图形
pdf(file="cor.pdf", width=8, height=8)
corrplot(M,
         order="original",
         method = "circle",
         type = "upper",
         tl.cex=0.8, pch=T,
         p.mat = res1$p,
         insig = "label_sig",
         pch.cex = 1.6,
         sig.level=0.05,
         number.cex = 1,
         col=colorRampPalette(c("blue", "white", "red"))(50),
         tl.col="black")
dev.off()


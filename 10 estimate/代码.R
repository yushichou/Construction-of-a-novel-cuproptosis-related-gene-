#引用包
rm(list = ls())
library(limma)
library(estimate)
library(tidyverse)
inputFile="merge.txt"                                          #输入文件名字

#读取文件,并对输入文件整理
load("merge.RDATA")
rt=outTab


#输出整理后的矩阵文件
rt=rbind(ID=colnames(rt),rt)
write.table(rt,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)

#运行estimate包
filterCommonGenes(input.f="uniq.symbol.txt", 
                  output.f="commonGenes.gct", 
                  id="GeneSymbol")

estimateScore(input.ds = "commonGenes.gct",
              output.ds="estimateScore.gct", 
              platform="illumina")

#输出每个样品的打分
scores=read.table("estimateScore.gct",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
rownames(scores)=gsub("\\.","\\-",rownames(scores))
out=rbind(ID=colnames(scores),scores)
write.table(out,file="scores.txt",sep="\t",quote=F,col.names=F)

#读取分型文件
clusterFile="Cluster.txt" 
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
row.names(scores)<-str_replace_all(row.names(scores),"-",".")

#数据合并
sameSample=intersect(row.names(scores), row.names(cluster))
scores=scores[sameSample,,drop=F]
cluster=cluster[sameSample,,drop=F]
scoreCluster=cbind(scores, cluster)
library(reshape2)
#把数据转换成ggplot2输入文件
data=melt(scoreCluster, id.vars=c("cluster"))
colnames(data)=c("cluster", "Immune", "Fraction")

#绘制箱线图
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data[,"cluster"])))]
library(tidyverse)
library(ggpubr)
p=ggboxplot(data, x="Immune", y="Fraction", color="cluster", 
            ylab="Immune infiltration",
            xlab="",
            legend.title="cluster",
            palette=bioCol)
p=p+rotate_x_text(50)
pdf(file="boxplot.pdf", width=4, height=4.5)                          #输出图片文件
p+stat_compare_means(aes(group=cluster),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
dev.off()














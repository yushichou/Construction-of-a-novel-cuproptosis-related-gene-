rm(list = ls())
#引用包
library(limma)
library(ggpubr)
CluFile="Cluster.txt"        #分型文件
geneCluFile="geneCluster.txt"      #基因分型文件
scoreFile="score.txt"           #打分文件

#读取输入文件
Clu=read.table(CluFile, header=T, sep="\t", check.names=F, row.names=1)
geneClu=read.table(geneCluFile, header=T, sep="\t", check.names=F, row.names=1)
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
twoCluster=cbind(Clu, geneClu)
sameSample=intersect(row.names(twoCluster), row.names(score))
data=cbind(score[sameSample,,drop=F], twoCluster[sameSample,,drop=F])

#######分型与打分相关性########
#设置比较组
data$cluster=factor(data$cluster, levels=levels(factor(data$cluster)))
group=levels(factor(data$cluster))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#定义颜色
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data$cluster)))]

#绘制boxplot
boxplot=ggboxplot(data, x="cluster", y="score", color="cluster",
                  xlab="cluster",
                  ylab="score",
                  legend.title="cluster",
                  palette=bioCol,
                  add = "jitter")+ 
  stat_compare_means(comparisons = my_comparisons)

#输出图片
pdf(file="cluster.pdf", width=5, height=4.5)
print(boxplot)
dev.off()
#######分型与打分相关性########


#######基因分型与打分相关性########
#设置比较组
data$geneCluster=factor(data$geneCluster, levels=levels(factor(data$geneCluster)))
group=levels(factor(data$geneCluster))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#定义颜色
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data$geneCluster)))]

#绘制boxplot
boxplot=ggboxplot(data, x="geneCluster", y="score", color="geneCluster",
                  xlab="geneCluster",
                  ylab="score",
                  legend.title="geneCluster",
                  palette=bioCol,
                  add = "jitter")+ 
  stat_compare_means(comparisons = my_comparisons)

#输出图片
pdf(file="geneCluster.pdf", width=5, height=4.5)
print(boxplot)
dev.off()
#######基因分型与打分相关性########


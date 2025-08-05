rm(list = ls())
#引用包
library(ggalluvial)
library(ggplot2)
library(dplyr)
CluFile="Cluster.txt"         #分型文件
geneCluFile="geneCluster.txt"       #基因分型文件
scoreFile="score.group.txt"      #打分的分组文件
cliFile="clinicaldata.txt"              #临床数据文件
trait="fustat"                      #临床性状

#读取输入文件
Clu=read.table(CluFile, header=T, sep="\t", check.names=F, row.names=1)
geneClu=read.table(geneCluFile, header=T, sep="\t", check.names=F, row.names=1)
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

cli$fustat[which(cli$fustat == 1)]="Dead"
cli$fustat[which(cli$fustat == 0)]="Alive"
#合并数据
twoCluster=cbind(Clu, geneClu)
library(tidyverse)
rownames(twoCluster)=str_replace_all(rownames(twoCluster),"TCGA_","")
rownames(twoCluster)=str_replace_all(rownames(twoCluster),"GSE63514_","")

sameSample=intersect(row.names(twoCluster), row.names(score))
scoreClu=cbind(score[sameSample,,drop=F], twoCluster[sameSample,,drop=F])
sameSample=intersect(row.names(scoreClu), row.names(cli))
rt=cbind(scoreClu[sameSample,], cli[sameSample,])
rt$fustat[which(rt$fustat == 1)]="Dead"
rt$fustat[which(rt$fustat == 0)]="Alive"

#准备桑基图输入文件
rt=rt[,c("cluster", "geneCluster", "group", "fustat")]
colnames(rt)=c("cluster", "geneCluster", "score", "fustat")
rt=na.omit(rt)
corLodes=to_lodes_form(rt, axes = 1:ncol(rt), id = "Cohort")

#得到输出文件
pdf(file="ggalluvial.pdf", width=6, height=6)
mycol=rep(c("#0066FF","#FF9900","#FF0000","#029149","#6E568C","#E0367A","#D8D155","#223D6C","#D20A13","#431A3D","#91612D","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#64495D","#7CC767"),15)
ggplot(corLodes, aes(x = x, stratum = stratum, alluvium = Cohort,fill = stratum, label = stratum)) +
  scale_x_discrete(expand = c(0, 0)) +  
  #用aes.flow控制线条颜色，forward说明颜色和前面的柱状图一致，backward说明和后面的柱状图一致。
  geom_flow(width = 2/10,aes.flow = "forward") + 
  geom_stratum(alpha = .9,width = 2/10) +
  scale_fill_manual(values = mycol) +
  #size=3代表字体大小
  geom_text(stat = "stratum", size = 3,color="black") +
  xlab("") + ylab("") + theme_bw() + 
  theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text.y = element_blank()) + #去掉坐标轴
  theme(panel.grid =element_blank()) + 
  theme(panel.border = element_blank()) + 
  ggtitle("") + guides(fill = FALSE)                            
dev.off()


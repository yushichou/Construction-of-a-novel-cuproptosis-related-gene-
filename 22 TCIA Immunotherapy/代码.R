rm(list = ls())
library(ggpubr)                    #引用包
tciaFile="TCIA.txt"                #免疫治疗打分文件
scoreFile="score.group.txt"     #打分分组文件
library(tidyverse)  

#读取免疫治疗打分文件
ips=read.table(tciaFile, header=T, sep="\t", check.names=F, row.names=1)
rownames(ips)=str_replace_all(rownames(ips),"-",".")
#读取打分分组文件
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
library(tidyverse)

#合并数据
row.names(ips)=str_replace_all(row.names(ips),"-",".")
aaa=score
aaa$ID=row.names(score)
aaa$ID=substr(aaa$ID,1,12)
aaa<-distinct(aaa,ID,.keep_all = T) #去重
rownames(aaa)=aaa$ID
score=aaa
sameSample=intersect(row.names(ips), row.names(score))

ips=ips[sameSample, , drop=F]
score=score[sameSample, "group", drop=F]
data=cbind(ips, score)

#设置比较组
data$group=factor(data$group, levels=c("Low", "High"))
group=levels(factor(data$group))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#对免疫治疗打分进行循环,分别绘制小提琴图
for(i in colnames(data)[1:(ncol(data)-1)]){
  rt=data[,c(i, "group")]
  colnames(rt)=c("IPS", "group")
  gg1=ggviolin(rt, x="group", y="IPS", fill = "group", 
               xlab="score", ylab=i,
               legend.title="score",
               palette=c("#00AFBB", "#E7B800"),
               add = "boxplot", add.params = list(fill="white"))+ 
    stat_compare_means(comparisons = my_comparisons)
  #stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
  
  pdf(file=paste0(i, ".pdf"), width=6, height=5)
  print(gg1)
  dev.off()
}



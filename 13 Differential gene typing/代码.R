
rm(list = ls())
#引用包
library(limma)
library(ConsensusClusterPlus)
expFile="uniSigGeneExp.txt"        #表达输入文件

#读取输入文件，并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#聚类
maxK=9
results=ConsensusClusterPlus(data,
                             maxK=maxK,
                             reps=50,
                             pItem=0.8,
                             pFeature=1,
                             clusterAlg="km",
                             distance="euclidean",
                             seed=123456,
                             plot="png")

#输出分型结果
clusterNum=2        #分几类，根据判断标准判断
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("geneCluster")
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(cluster$geneCluster))
cluster$geneCluster=letter[match(cluster$geneCluster, uniqClu)]
clusterOut=rbind(ID=colnames(cluster), cluster)
write.table(clusterOut, file="geneCluster.txt", sep="\t", quote=F, col.names=F)



#引用包
library(survival)
library(survminer)
clusterFile="geneCluster.txt"     #基因分型文件
cliFile="clinicaldata.txt"                #生存数据文件

#读取输入文件
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
library(tidyverse)
rownames(cluster)=str_replace_all(rownames(cluster),"TCGA_","")
rownames(cluster)=str_replace_all(rownames(cluster),"GSE63514_","")

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[,1:2]
cli$futime=cli$futime/365

#数据合并
sameSample=intersect(row.names(cluster), row.names(cli))
rt=cbind(cli[sameSample,,drop=F], cluster[sameSample,,drop=F])
rt$futime=as.numeric(rt$futime)
rt=dplyr::filter(rt,futime < 8)
#生存差异统计
length=length(levels(factor(rt$geneCluster)))
diff=survdiff(Surv(futime, fustat) ~ geneCluster, data = rt)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ geneCluster, data = rt)
#print(surv_median(fit))

#绘制生存曲线
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(rt[,"geneCluster"])))]
surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=F,
                   pval=pValue,
                   pval.size=6,
                   legend.labs=levels(factor(rt[,"geneCluster"])),
                   legend.title="geneCluster",
                   legend = c(0.8, 0.8),
                   font.legend=10,
                   xlab="Time(years)",
                   break.time.by = 1,
                   palette = bioCol,
                   surv.median.line = "hv",
                   risk.table=T,
                   cumevents=F,
                   risk.table.height=.25)
pdf(file="survival.pdf", onefile = FALSE, width=7, height=5.5)
print(surPlot)
dev.off()



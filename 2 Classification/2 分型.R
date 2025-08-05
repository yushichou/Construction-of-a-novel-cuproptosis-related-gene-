rm(list = ls())
load("merge.RDATA")
Merge<-outTab
celltype="Cuproptosis"

genename<-read.table(paste0(celltype,".txt"),header = T,sep = "\t")
AAA<-intersect(genename[,1],rownames(Merge))
geneexpr<-Merge[AAA,]

write.table(geneexpr,"geneexpr.txt",quote = F,row.names =T ,sep = "\t")
save(geneexpr,file = "geneexpr.RDATA")
library(ConsensusClusterPlus)


#读取输入文件
data=geneexpr
data=as.matrix(data)

#聚类
maxK=9
results=ConsensusClusterPlus(data,
                             maxK=maxK,
                             reps=50,
                             pItem=0.8,
                             pFeature=1,
                             clusterAlg="pam",
                             distance="euclidean",
                             seed=123456,
                             title = celltype,
                             plot="png")


#输出分型结果
clusterNum=2      #分几类，根据判断标准判断
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("cluster")
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(cluster$cluster))
cluster$cluster=letter[match(cluster$cluster, uniqClu)]
clusterOut=rbind(ID=colnames(cluster), cluster)
write.table(clusterOut, file="Cluster.txt", sep="\t", quote=F, col.names=F)


#引用包
library(survival)
library(survminer)
library(tidyverse)
clusterFile="Cluster.txt"     #分型文件
cliFile="clinicaldata.txt"               #生存数据文件


#读取输入文件
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
rownames(cluster)=str_replace_all(rownames(cluster),"TCGA_","")
rownames(cluster)=str_replace_all(rownames(cluster),"GSE54460_","")

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli<-cli[,1:2]
cli$futime=cli$futime/365

#数据合并
sameSample=intersect(row.names(cluster), row.names(cli))

rt=cbind(cli[sameSample,,drop=F], cluster[sameSample,,drop=F])
#rt=dplyr::filter(rt,futime < 8)
#生存差异统计
#rt=dplyr::filter(rt,cluster != "B")
length=length(levels(factor(rt$cluster)))
diff=survdiff(Surv(futime, fustat) ~ cluster, data = rt)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ cluster, data = rt)
#print(surv_median(fit))

#绘制生存曲线
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=F,
                   pval=pValue,
                   pval.size=6,
                   legend.title="cluster",
                   legend.labs=levels(factor(rt[,"cluster"])),
                   legend = c(0.8, 0.8),
                   font.legend=10,
                   xlab="Time(years)",
                   break.time.by = 1,
                   palette = bioCol,
                   surv.median.line = "hv",
                   risk.table=T,
                   cumevents=F,
                   risk.table.height=.25)
surPlot
pdf(file="survival.pdf",onefile = FALSE,width=7,height=5.5)
print(surPlot)
dev.off()

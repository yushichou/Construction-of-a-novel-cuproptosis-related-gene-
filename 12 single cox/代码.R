
rm(list = ls())
#引用包
library(limma)
library(survival)
library(tidyverse)
expFile="interGeneExp.txt"     #差异基因表达文件
cliFile="clinicaldata.txt"             #生存数据文件

#读取表达文件，并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data1=t(data)
library(tidyverse)
rownames(data1)=str_replace_all(rownames(data1),"TCGA_","")
rownames(data1)=str_replace_all(rownames(data1),"GSE63514_","")

#读取生存数据
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)     #读取临床文件
cli=cli[,1:2]
cli$futime=cli$futime/365

#数据合并
sameSample=intersect(row.names(data1), row.names(cli))
data1=data1[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
rt=cbind(cli, data1)

#对基因进行循环，找出预后相关的基因
outTab=data.frame()
sigGenes=c()
for(i in colnames(rt[,3:ncol(rt)])){
  #cox分析
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<100){
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}

#输出单因素的结果
write.table(outTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)
outTab$pvalue=as.numeric(outTab$pvalue)
outTab=arrange(outTab,pvalue)
outTab=dplyr::filter(outTab,pvalue < 0.05)
sigGenes=outTab$id
#保存单因素显著基因的表达量
sigGeneExp=data[sigGenes,]
sigGeneExp=rbind(id=colnames(sigGeneExp), sigGeneExp)
write.table(sigGeneExp, file="uniSigGeneExp.txt", sep="\t", quote=F, col.names=F)


library(survminer)
library(survival)
library(survivalROC)
library(forestplot)

rt=outTab
rt=dplyr::filter(rt,id %in% sigGenes)
row.names(rt)=rt$id
rt=rt[,-1]
rt<-na.omit(rt)
rt<-arrange(rt,pvalue)
options(forestplot_new_page = FALSE)
clrs <- fpColors(box="red",line="darkblue", summary="royalblue") 
data=as.matrix(rt)
HR=as.data.frame(data[,1:3])
HR[,1]=as.numeric(HR[,1])
HR[,2]=as.numeric(HR[,2])
HR[,3]=as.numeric(HR[,3])
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=as.numeric(pVal)
pVal=ifelse(pVal<0.05, "<0.05", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )   
pdf(file="OS_forest.pdf",
    width = 6,            
    height = 8,           
)
forestplot(tabletext, 
           zero = 1,
           lwd.zero = 2,
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=T,
           lwd.ci=2,
           boxsize=0.2,
           xlab="Hazard ratio"
)
dev.off()

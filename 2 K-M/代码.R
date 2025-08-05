rm(list = ls())
library(limma)           #引用包
expFile="merge.txt"      #表达输入文件
geneFile="Cuproptosis.txt"      #基因列表文件

#读取输入文件，并对数据进行处理
load("merge.RDATA")
rt=outTab
data=outTab

gene=read.table(geneFile, header=T, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(data))
geneExp=data[sameGene,]
setdiff(gene[,1],rownames(data))
#输出结果
out=rbind(ID=colnames(geneExp),geneExp)
write.table(out,file="GeneExp.txt",sep="\t",quote=F,col.names=F)


#引用包
library(limma)
library(survival)
library(survminer)
expFile="GeneExp.txt"     #表达数据文件
cliFile="clinicaldata.txt"           #生存数据文件

#读取表达文件，并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=t(data)
library(tidyverse)
rownames(data)=str_replace_all(rownames(data),"TCGA_","")
rownames(data)=str_replace_all(rownames(data),"GSE63514_","")

#读取生存数据
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)     #读取临床文件
cli=cli[,1:2]
cli$futime=cli$futime/365
#数据合并
sameSample=intersect(row.names(data), row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
rt=cbind(cli, data)

#对基因进行循环，找出预后相关的基因
outTab=data.frame()
km=c()
for(i in colnames(rt[,3:ncol(rt)])){
  #cox分析
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
  #km分析
  data=rt[,c("futime", "fustat", i)]
  colnames(data)=c("futime", "fustat", "gene")
  #获取最优cutoff
  res.cut=surv_cutpoint(data, time = "futime", event = "fustat", variables =c("gene"))
  res.cat=surv_categorize(res.cut)
  fit=survfit(Surv(futime, fustat) ~gene, data = res.cat)
  #print(paste0(i, " ", res.cut$cutpoint[1]))
  #比较高低表达生存差异
  diff=survdiff(Surv(futime, fustat) ~gene,data =res.cat)
  pValue=1-pchisq(diff$chisq, df=1)
  km=c(km, pValue)
  #对pvalue<0.05的基因绘制生存曲线
  if(pValue<0.05){
    if(pValue<0.001){
      pValue="p<0.001"
    }else{
      pValue=paste0("p=",sprintf("%.03f",pValue))
    }
    
    #绘制生存曲线
    surPlot=ggsurvplot(fit,
                       data=res.cat,
                       pval=pValue,
                       pval.size=6,
                       legend.title=i,
                       legend.labs=c("high","low"),
                       xlab="Time(years)",
                       ylab="Overall survival",
                       palette=c("red", "blue"),
                       break.time.by=1,
                       conf.int=T,
                       risk.table=TRUE,
                       risk.table.title="",
                       risk.table.height=.25)
    pdf(file=paste0("sur.", i, ".pdf"),onefile = FALSE,
        width = 8,         #图片的宽度
        height =5)         #图片的高度
    print(surPlot)
    dev.off()
  }
}

#输出单因素的结果
outTab=cbind(outTab, km)
write.table(outTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)



#引用包
library(igraph)
library(psych)
library(reshape2)
library("RColorBrewer")
GeneExpfile <- "GeneExp.txt"     #表达数据文件
Genefile <- "Cuproptosis.txt"              #基因列表文件
Coxfile <- "uniCox.txt"             #预后结果文件

#读取输入文件
gene.group <- read.table(Genefile,header=T,sep="\t")
colnames(gene.group) <- c('id')
gene.group$group<-"Cuproptosis"
gene.exp <- read.table(GeneExpfile,header=T,sep="\t",row.names=1)
gene.cox <- read.table(Coxfile,header=T,sep="\t")

#基因取交集

genelist <- intersect(gene.group$id, gene.cox$id)
genelist <- intersect(genelist, rownames(gene.exp))
gene.group <- gene.group[match(genelist,gene.group$id),]
gene.group <- gene.group[order(gene.group$group),]
gene.exp <- gene.exp[match(gene.group$id,rownames(gene.exp)),]
gene.cox <- gene.cox[match(gene.group$id,gene.cox$id),]

#准备网络文件
gene.cor <- corr.test(t(gene.exp))
gene.cor.cor <- gene.cor$r
gene.cor.pvalue <- gene.cor$p
gene.cor.cor[upper.tri(gene.cor.cor)] = NA
gene.cor.pvalue[upper.tri(gene.cor.pvalue)] = NA
gene.cor.cor.melt <- melt(gene.cor.cor)   #gene1 \t gene2 \t cor
gene.cor.pvalue.melt <- melt(gene.cor.pvalue)
gene.melt <- data.frame(from = gene.cor.cor.melt$Var2,to=gene.cor.cor.melt$Var1,cor=gene.cor.cor.melt$value,pvalue=gene.cor.pvalue.melt$value)
gene.melt <- gene.melt[gene.melt$from!=gene.melt$to&!is.na(gene.melt$pvalue),,drop=F]
#选取p小于0.0001的关系对
gene.edge <- gene.melt[gene.melt$pvalue<0.0001,,drop=F]
gene.edge$color <- ifelse(gene.edge$cor>0,'pink','#6495ED')
gene.edge$weight <- abs(gene.edge$cor)*6

#准备结果属性文件
gene.node <- gene.group
group.color <- colorRampPalette(brewer.pal(9, "Set1"))(length(unique(gene.node$group)))
gene.node$color <- group.color[as.numeric(as.factor(gene.node$group))]
gene.node$shape <- "circle"
gene.node$frame <- ifelse(gene.cox$HR>1,'purple',"green")
gene.node$pvalue <- gene.cox$pvalue
# pvalue size
pvalue.breaks <- c(0,0.0001,0.001,0.01,0.05,1)
pvalue.size <- c(16,14,12,10,8)
cutpvalue <- cut(gene.node$pvalue,breaks=pvalue.breaks)
gene.node$size <- pvalue.size[as.numeric(cutpvalue)]
nodefile <- "network.node.txt"
edgefile <- "network.edge.txt"
write.table(gene.node, nodefile, sep="\t", col.names=T, row.names=F, quote=F)
write.table(gene.edge, edgefile, sep="\t", col.names=T, row.names=F, quote=F)


#绘制网络图
node = read.table(nodefile, header=T, sep="\t", comment.char="")
edge = read.table(edgefile, header=T, sep="\t", comment.char="")

g = graph.data.frame(edge,directed = FALSE)
node = node[match(names(components(g)$membership),node$id),]

if(!is.na(match('color',colnames(node)))) V(g)$color = node$color
if(!is.na(match('size',colnames(node)))) V(g)$size = node$size
if(!is.na(match('shape',colnames(node)))) V(g)$shape = node$shape
if(!is.na(match('frame',colnames(node)))) V(g)$frame = node$frame

#输出文件
pdf(file="network.pdf", width=15, height=12)
par(mar=c(0,0,0,0))
layout(matrix(c(1,1,4,2,3,4),nc=2),height=c(4,4,2),width=c(8,3))

#节点坐标 
coord = layout_in_circle(g)
degree.x = acos(coord[,1])
degree.y = asin(coord[,2])
degree.alpha = c()
for(i in 1:length(degree.x)){
  if(degree.y[i]<0) degree.alpha=c(degree.alpha,2*pi-degree.x[i]) else degree.alpha=c(degree.alpha,degree.x[i])
}
degree.cut.group = (0:8)/4*pi
degree.cut.group[1] = -0.0001
degree.cut = cut(degree.alpha,degree.cut.group)
degree.degree = c(-pi/4,-pi/4,-pi/2,-pi/2,pi/2,pi/2,pi/2,pi/4)
degree = degree.degree[as.numeric(degree.cut)]

#定义饼图,左半圆颜色代表的类型，右半圆代表基因是高风险基因,还是低风险基因
values <- lapply(node$id,function(x)c(1,1))
V(g)$pie.color = lapply(1:nrow(node),function(x)c(node$color[x],node$frame[x]))
V(g)$frame = NA 

#绘制图形
plot(g,layout=layout_in_circle,vertex.shape="pie",vertex.pie=values,
     vertex.label.cex=V(g)$lable.cex,edge.width = E(g)$weight,edge.arrow.size=0,
     vertex.label.color=V(g)$color,vertex.frame.color=V(g)$frame,edge.color=E(g)$color,
     vertex.label.cex=2,vertex.label.font=2,vertex.size=V(g)$size,edge.curved=0.4,
     vertex.color=V(g)$color,vertex.label.dist=1,vertex.label.degree=degree)
# label.degree : zero means to the right; and pi means to the left; up is -pi/2 and down is pi/2;  The default value is -pi/4
# label.dist If it is 0 then the label is centered on the vertex; If it is 1 then the label is displayed beside the vertex.

#绘制节点属性图例(的类型)
par(mar=c(0,0,0,0))
plot(1,type="n",xlab="",ylab="",axes=F)
groupinfo = unique(data.frame(group=node$group,color=node$color))
legend("left",legend=groupinfo$group,col=groupinfo$color,pch=16,bty="n",cex=3)
#绘制风险图例(哪些基因是高风险的基因,哪些基因是低风险的基因)
par(mar=c(0,0,0,0))
plot(1,type="n",xlab="",ylab="",axes=F)
legend("left",legend=c('Risk factors','Favorable factors'),col=c('purple','green'),pch=16,bty="n",cex=2.5)
#绘制预后pvalue图例
par(mar=c(0,0,0,0))
plot(1,type="n",xlab="",axes=F,ylab="")
legend("top",legend=c('Postive correlation with P<0.0001','Negative correlation with P<0.0001'),lty=1,lwd=4,col=c('pink','#6495ED'),bty="n",cex=2.2)
legend('bottom',legend=c(0.0001,0.001,0.01,0.05,1),pch=16,pt.cex=c(1.6,1.4,1.2,1,0.8)*6,bty="n",ncol=5,cex=2.2,col="black",title="Cox test, pvalue")
dev.off()





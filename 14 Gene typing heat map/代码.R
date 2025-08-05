
rm(list = ls())
library(pheatmap)        #引用包
expFile="uniSigGeneExp.txt"        #表达数据文件
geneCluFile="geneCluster.txt"      #基因分型结果文件
CluFile="Cluster.txt"        #分型结果文件
cliFile="clinicaldata.txt"             #临床数据文件


#读取输入文件
exp=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
Clu=read.table(CluFile, header=T, sep="\t", check.names=F, row.names=1)
geneClu=read.table(geneCluFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
exp=as.data.frame(t(exp))
sameSample=intersect(row.names(exp), row.names(Clu))
exp=exp[sameSample,,drop=F]
expData=cbind(exp, geneCluster=geneClu[sameSample,], cluster=Clu[sameSample,])
Project=gsub("(.*?)\\_.*", "\\1", rownames(expData))
library(tidyverse)
rownames(expData)=str_replace_all(rownames(expData),"TCGA_","")
rownames(expData)=str_replace_all(rownames(expData),"GSE63514_","")

expData=cbind(expData, Project)

#合并临床数据
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli$fustat[cli$fustat==1]="Dead"
cli$fustat[cli$fustat==0]="Alive"
sameSample=intersect(row.names(expData), row.names(cli))
expData=expData[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
data=cbind(expData, cli)
colnames(data)[which(colnames(data)=="cluster")]="cluster"


#提取热图数据
data=data[order(data$geneCluster),]
Type=data[,((ncol(data)-2-ncol(cli)):ncol(data))]
data=t(data[,1:(ncol(expData)-3)])

#聚类颜色
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
ann_colors=list()
col=bioCol[1:length(levels(factor(Type$cluster)))]
names(col)=levels(factor(Type$cluster))
ann_colors[["cluster"]]=col
GENEcol=bioCol[1:length(levels(factor(Type$geneCluster)))]
names(GENEcol)=levels(factor(Type$geneCluster))
ann_colors[["geneCluster"]]=GENEcol

#热图可视化
pdf("heatmap.pdf", height=6, width=8)
pheatmap(data,
         annotation=Type,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("blue",5), "black", rep("yellow",5)))(50),
         cluster_cols =F,
         cluster_rows =F,
         scale="row",
         show_colnames=F,
         show_rownames=T,
         fontsize=6,
         fontsize_row=8,
         fontsize_col=6)
dev.off()




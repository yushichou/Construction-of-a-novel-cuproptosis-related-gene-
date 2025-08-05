rm(list = ls())
#引用包
library(limma) 
library(VennDiagram)
expFile="merge.txt"           #表达输入文件
cluFile="Cluster.txt"      #聚类结果文件


#读取输入文件，并对输入文件整理
load("merge.RDATA")
rt=outTab
exp=outTab
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)


#读取cluster文件
cluster=read.table(cluFile, header=T, sep="\t", check.names=F, row.names=1)

#提取交集文件
sameSample=intersect(colnames(data), row.names(cluster))
data=data[,sameSample]
cluster=cluster[sameSample,]

#差异分析
logFCfilter=1
geneList=list()
Type=as.vector(cluster)
design=model.matrix(~0+factor(Type))
colnames(design)=levels(factor(Type))
comp=combn(levels(factor(Type)), 2)
allDiffGenes=c()
for(i in 1:ncol(comp)){
  fit=lmFit(data, design)
  contrast=paste0(comp[2,i], "-", comp[1,i])
  #print(contrast)
  cont.matrix=makeContrasts(contrast, levels=design)
  fit2=contrasts.fit(fit, cont.matrix)
  fit2=eBayes(fit2)
  
  #输出所有基因差异情况
  allDiff=topTable(fit2,adjust='fdr',number=200000)
  allDiffOut=rbind(id=colnames(allDiff),allDiff)
  write.table(allDiffOut, file=paste0(contrast, ".all.txt"), sep="\t", quote=F, col.names=F)
  
  #输出差异结果
  diffSig=allDiff[with(allDiff, (abs(logFC)>logFCfilter & P.Value < 0.05 )), ]
  diffSigOut=rbind(id=colnames(diffSig),diffSig)
  write.table(diffSigOut, file=paste0(contrast, ".diff.txt"), sep="\t", quote=F, col.names=F)
  geneList[[contrast]]=row.names(diffSig)
}


#保存并集基因

B.A=geneList$`B-A`


interGenes=unique(c(B.A))
write.table(file="interGene.txt",interGenes,sep="\t",quote=F,col.names=F,row.names=F)

#保存交集基因的表达量
interGeneExp=data[interGenes,]
interGeneExp=rbind(id=colnames(interGeneExp), interGeneExp)
write.table(interGeneExp, file="interGeneExp.txt", sep="\t", quote=F, col.names=F)


volcano<-read.table("B-A.all.txt",header = T,sep = "\t",row.names = 1)
p.cut=0.05
logFC.cut=1
volcano$type[(volcano$P.Value > p.cut|volcano$P.Value=="NA")|(volcano$logFC < logFC.cut)& volcano$logFC > -logFC.cut] <- "none significant"
volcano$type[volcano$P.Value <= p.cut & volcano$logFC >= logFC.cut] <- "up-regulated"
volcano$type[volcano$P.Value <= p.cut & volcano$logFC <= -logFC.cut] <- "down-regulated"  
library(tidyverse)
p = ggplot(volcano,aes(logFC,-1*log10(P.Value),color=type))       
p + geom_point()

x_lim <- max(volcano$logFC,-allDiff$logFC) 
gg=p + geom_point( aes(size = abs(logFC)),alpha = 0.4) + xlim(-x_lim,x_lim) +labs(x="log2(FC)",y="-log10(P.Value)")+
  scale_color_manual(values =c("blue","grey","red"))+
  geom_hline(aes(yintercept=-1*log10(p.cut)),colour="black", linetype="dashed") +
  geom_vline(xintercept=c(-logFC.cut,logFC.cut),colour="black", linetype="dashed")
print(gg) 

pdf(file="B-A火山图.PDF", width=6, height=6)
gg
dev.off()


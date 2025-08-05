rm(list = ls())
#引用包
library(limma)
library(sva)
library(ggplot2)
library(factoextra)
library(FactoMineR)
files=c("TCGA.txt","GSE63514.txt")       

#获取交集基因
geneList=list()
for(i in 1:length(files)){
  inputFile=files[i]
  rt=read.table(inputFile, header=T, sep="\t",check.names=F)
  header=unlist(strsplit(inputFile, "\\.|\\-"))
  geneList[[header[1]]]=as.vector(rt[,1])
}
intersectGenes=Reduce(intersect, geneList)

#数据合并
allTab=data.frame()
batchType=c()
for(i in 1:length(files)){
  inputFile=files[i]
  header=unlist(strsplit(inputFile, "\\.|\\-"))
  #读取输入文件，并对输入文件进行整理
  rt=read.table(inputFile, header=T, sep="\t", check.names=F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  rt=avereps(data)
  colnames(rt)=paste0(header[1], "_", colnames(rt))

  #对数值大的数据取log2
  qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
  if(LogC){
    rt[rt<0]=0
    rt=log2(rt+1)}
  if(header[1] != "TCGA"){
    rt=normalizeBetweenArrays(rt)
  }
  #数据合并
  if(i==1){
    allTab=rt[intersectGenes,]
  }else{
    allTab=cbind(allTab, rt[intersectGenes,])
  }
  batchType=c(batchType, rep(i,ncol(rt)))
}

#对数据进行矫正，输出矫正后的结果
outTab=ComBat(allTab, batchType, par.prior=TRUE)

library(FactoMineR)
library(factoextra)
library(pca3d) # 加载R包
pca <- prcomp(t(allTab),scale. = TRUE) # 使用R自带的主成分分析函数
pheno<-data.frame(ID=colnames(allTab))
pheno$cancer<-pheno$ID
pheno[1:306,2]<-"TCGA"
pheno[307:334,2]<-"GSE31210"

gr <- factor(pheno$cancer) 

pca3d(pca, group = gr, show.centroids = T, legend = 'topleft') # 画图
snapshotPCA3d(file="去除批次效应前.png")
library(FactoMineR)
library(factoextra)
pca <- prcomp(t(outTab),scale. = TRUE) # 使用R自带的主成分分析函数
pca3d(pca, group = gr, show.centroids = T, legend = 'topleft') # 画图
snapshotPCA3d(file="去除批次效应后.png")

save(outTab,file = "merge.RDATA")
outTab=rbind(geneNames=colnames(outTab), outTab)
write.table(outTab, file="merge.txt", sep="\t", quote=F, col.names=F)




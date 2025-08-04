rm(list = ls())
#引用包
library(limma)
library(sva)
library(ggplot2)
library(factoextra)
library(FactoMineR)
# 直接指定文件完整路径
setwd("D:/科研/王晓 宫颈癌 铜死亡分型 结果+组图/1 融合数据")
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



# PCA 矫正前
pca_before <- prcomp(t(allTab), scale. = TRUE)
group <- factor(pheno$cancer)
p1 <- fviz_pca_ind(pca_before, col.ind = group,   geom = "point",
                   addEllipses = TRUE, ellipse.level = 0.95,
                   legend.title = "Dataset") +
  ggtitle("PCA Before Batch Effect Removal")
ggsave("PCA_before_ComBat.pdf", plot = p1)



library(FactoMineR)
library(factoextra)
pca_after <- prcomp(t(outTab), scale. = TRUE)
p2 <- fviz_pca_ind(pca_after, col.ind = group,   geom = "point",
                   addEllipses = TRUE, ellipse.level = 0.95,
                   legend.title = "Dataset") +
  ggtitle("PCA After Batch Effect Removal")
ggsave("PCA_after_ComBat.pdf", plot = p2)
#绘制UMAP 图 
library(umap)
set.seed(123)
umap_before <- umap(t(allTab))
umap_after <- umap(t(outTab))

umap_df_before <- data.frame(UMAP1 = umap_before$layout[,1],
                             UMAP2 = umap_before$layout[,2],
                             Group = group)
ggplot(umap_df_before, aes(UMAP1, UMAP2, color = Group)) +
  geom_point(size = 2) + ggtitle("UMAP Before ComBat") +
  theme_minimal()
ggsave("UMAP_before_ComBat.pdf")

umap_df_after <- data.frame(UMAP1 = umap_after$layout[,1],
                            UMAP2 = umap_after$layout[,2],
                            Group = group)
ggplot(umap_df_after, aes(UMAP1, UMAP2, color = Group)) +
  geom_point(size = 2) + ggtitle("UMAP After ComBat") +
  theme_minimal()
ggsave("UMAP_after_ComBat.pdf")


# 计算 PCA 方差比例
var_explained_before <- summary(pca_before)$importance[2, 1:2]
var_explained_after <- summary(pca_after)$importance[2, 1:2]

cat("Before ComBat: PC1 explains", round(var_explained_before[1]*100,2), "% variance\n")
cat("After ComBat: PC1 explains", round(var_explained_after[1]*100,2), "% variance\n")



save(outTab,file = "merge.RDATA")
outTab=rbind(geneNames=colnames(outTab), outTab)
write.table(outTab, file="merge 1.txt", sep="\t", quote=F, col.names=F)




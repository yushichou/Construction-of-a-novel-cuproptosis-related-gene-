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
# 读取文件获取基因列表
geneList <- list()
for(i in 1:length(files)){
  rt <- read.table(files[i], header=TRUE, sep="\t", check.names=FALSE)
  genes <- as.character(rt[,1])
  genes <- trimws(genes)  # 去除前后空格
  geneList[[i]] <- genes
}

# 取交集
intersectGenes <- Reduce(intersect, geneList)

# 输出交集基因数量和示例
cat("交集基因数:", length(intersectGenes), "\n")
print(head(intersectGenes))
#数据合并
allTab <- list()
batchType <- c()

for(i in 1:length(files)){
  inputFile <- files[i]
  header <- unlist(strsplit(inputFile, "\\.|\\-"))
  dataset_name <- header[1]   # TCGA or GSE63514
  
  rt <- read.table(inputFile, header=TRUE, sep="\t", check.names=FALSE)
  rt <- as.matrix(rt)
  rownames(rt) <- rt[,1]
  exp <- rt[, -1]
  data <- matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=list(rownames(exp), colnames(exp)))
  rt <- avereps(data)
  
  # 处理log转换
  qx <- as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=TRUE))
  LogC <- (qx[5] > 100) || ((qx[6] - qx[1]) > 50 && qx[2] > 0)
  if(LogC){
    rt[rt < 0] <- 0
    rt <- log2(rt + 1)
  }
  
  # GSE数据归一化
  if(dataset_name != "TCGA"){
    rt <- normalizeBetweenArrays(rt)
  }
  
  # 修改列名加前缀
  colnames(rt) <- paste0(dataset_name, "_", colnames(rt))
  
  # 筛选交集基因
  rt <- rt[intersectGenes, ]
  
  # 添加到列表
  allTab[[i]] <- rt
  
  # 正确添加批次信息
  batchType <- c(batchType, rep(dataset_name, ncol(rt)))
}

# 合并数据矩阵
allTab <- do.call(cbind, allTab)
batchType <- as.factor(batchType)

outTab <- ComBat(dat = as.matrix(allTab), batch = batchType, mod = NULL, par.prior = TRUE)


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
p1 <- fviz_pca_ind(pca_before, geom = "point",
                   habillage = group, # 用于根据 group 着色
                   addEllipses = TRUE, ellipse.level = 0.95,
                   palette = c("blue", "#FF9900"),  # 设置颜色
                   legend.title = "Dataset") +
  ggtitle("PCA Before Batch Effect Removal")
ggsave("PCA_before_ComBat.pdf", plot = p1)


library(FactoMineR)
library(factoextra)
pca_after <- prcomp(t(outTab), scale. = TRUE)
p2 <- fviz_pca_ind(pca_after, geom = "point",
                   habillage = group, # 用于根据 group 着色
                   addEllipses = TRUE, ellipse.level = 0.95,
                   palette = c("blue", "#FF9900"),  # 设置颜色
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




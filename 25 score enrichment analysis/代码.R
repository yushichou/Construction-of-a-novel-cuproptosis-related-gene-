
rm(list = ls())
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggsci)
library(scater)
library(SingleR)
library(hdf5r)
library(tidyverse)
library(RColorBrewer)
library(dplyr)
library(Seurat)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(limma)
library(ggplot2)
load("merge.RDATA")

score=read.table("score.txt",header = T,sep = "\t")



# 导入gmt文件，这里以MsigDB中的Hallmark为例
# MsigDB链接(https://www.gsea-msigdb.org/gsea/msigdb/)
# 下载文件时，注意选择对应的gene格式吗，一般10X的结果，都是gene symbol格式的。
genesets <- getGmt("h.all.v7.5.1.symbols.gmt")
str(genesets)

gsvascore <- gsva(outTab, 
                  genesets, 
                  min.sz=10, 
                  max.sz=500, 
                  verbose=TRUE,
                  parallel.sz=1)


aaa=as.data.frame(t(gsvascore))
aaa$id=rownames(aaa)
aaa=inner_join(score,aaa,by="id")
#不分细胞，全部进行相关性分析
gene="score"
y <- as.numeric(aaa[,gene])#开始相关性分析
bbb=aaa[,3:52]
cor_data_df <- data.frame(colnames(bbb))
for (a in 1:length(colnames(bbb))){
  test <- cor.test(as.numeric(bbb[,a]),y,method = "pearson")
  cor_data_df[a,2] <- test$estimate
  cor_data_df[a,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")
write.csv (cor_data_df, file =paste(gene,'GSVA相关性分析.csv',sep=' '), row.names =FALSE)#将文件导出
#相关性分析结果排序选top15的做相关性图
cor_data_df<-na.omit(cor_data_df)
posneg<-cor_data_df %>%
  filter(symbol != gene)
posneg$group<-NA
posneg$group[ posneg$pvalue<0.05& posneg$correlation<0]<-"negative"
posneg$group[ posneg$pvalue>=0.05]<-"pvalue>0.05"
posneg$group[ posneg$pvalue<0.05& posneg$correlation>0]<-"positive"
posneg$symbol<-str_replace_all(posneg$symbol,'_',' ')
posneg$symbol<-str_replace_all(posneg$symbol,'HALLMARK ','')
library(ggpubr)
p<-ggbarplot(posneg, x = "symbol", y = "correlation",
             fill = "group",           # change fill color by mpg_level
             color = "white",            # Set bar border colors to white
             palette = "jco",            # jco journal color palett. see ?ggpar
             sort.val = "asc",           # Sort the value in ascending order
             sort.by.groups = FALSE,     # Don't sort inside each group
             #x.text.angle = 90,          # Rotate vertically x axis texts
             ylab = "GSVA correlation",
             xlab = "Term",
             legend.title = "group",
             lab.col = "black",
             rotate = TRUE,
             ggtheme = theme_minimal()
)
ggsave(p,filename = paste("ALL cell", gene,'GSVA相关性图.pdf',sep=' '),width = 10,height = 10)




















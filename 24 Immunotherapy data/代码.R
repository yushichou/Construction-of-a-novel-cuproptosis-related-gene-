#打分  #第3430行看突变的基因名是否需要变更  #3673开始免疫治疗
rm(list = ls())
library(GSEABase)
library(GSVAdata)
library(Biobase)
library(genefilter)
library(limma)
library(RColorBrewer)
library(GSVA)
library(clusterProfiler)
library(pheatmap)
library(tidyverse)
library(survminer)
library(survival)
library(limma)
library(ggpubr)
library(stringr)
library(RColorBrewer)
library(export)
library(survivalROC)
library(pheatmap)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
library(ReactomePA)
library(igraph)
library(ggraph)
library(forestplot)
library(ggradar)
library(fmsb)
library(circlize)
library(ggsci)
library(parallel)
library(maftools)
library(patchwork)
library(ComplexHeatmap)
options(stringsAsFactors = FALSE) #禁止chr转成factor
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
gene="Score"  #需要替换
dir.create(paste0("Y:\\lcd代码\\TCGA随意做\\TCGA任意做\\进阶版\\进阶/",gene))
dir.create(paste0("Y:\\lcd代码\\TCGA随意做\\TCGA任意做\\进阶版\\进阶/",gene,"/5免疫浸润分析"))
dir.create(paste0("Y:\\lcd代码\\TCGA随意做\\TCGA任意做\\进阶版\\进阶/",gene,"/5免疫浸润分析/6免疫治疗"))

setwd("Y:\\lcd代码\\TCGA随意做\\TCGA任意做\\进阶版\\进阶\\原始数据\\免疫治疗数据")
###################################################################################
genelist_input <- fread(file="GSE32894.csv", header = T, sep=',', data.table = F)
genename <- as.character(genelist_input[,1]) #提取第一列基因名
#x - 基因组注释R包
#keys - 需要转换的基因列表
#keytype - 基因名类型
#columns - 希望返回的基因名类型数据，如返回NCBI基因ID使用 ENTREZID
gene_map <- biomaRt::select(org.Hs.eg.db, keys=genename, keytype="ENTREZID", columns=c("SYMBOL"))
colnames(gene_map)[1]<-"ID"
colnames(genelist_input)[1]<-"ID"
genelist_input$ID<-as.character(genelist_input$ID)
immudata<-inner_join(gene_map,genelist_input,by="ID")
immudata<-immudata[,-1]
immudata<-na.omit(immudata)
rownames(immudata)<-immudata$SYMBOL
immudata<-immudata[,-1]
immudata<-as.data.frame(t(immudata))

CCC=read.table("uniSigGeneExp.txt",header = T,sep = "\t")
CCC=CCC$id
CCC=intersect(colnames(immudata),CCC)
AAA=immudata[,CCC]
#PCA分析
pca=prcomp(AAA, scale=TRUE)
value=predict(pca)
score=value[,1]+value[,2]
score=as.data.frame(score)
score$ID=rownames(score)
immudata$ID=row.names(immudata)
immudata=inner_join(immudata,score,by="ID")
colnames(immudata)[which(colnames(immudata)=="score")]=gene
immudata<-immudata %>%
  dplyr::select(ID,gene,everything())


clindata<-read.csv("GSE32894_clin.csv",header = T,sep = ",")
colnames(clindata)[1]<-"ID"
AAA<-intersect(clindata$ID,immudata$ID)
immudata<-dplyr::filter(immudata,ID %in% AAA)
clindata<-dplyr::filter(clindata,ID %in% AAA)
rt<-inner_join(clindata,immudata,by="ID")

rt$time<-rt$time/12
rt$status<-as.numeric(rt$status)

#生存分析
setwd(paste0("Y:\\lcd代码\\TCGA随意做\\TCGA任意做\\进阶版\\进阶/",gene,"/5免疫浸润分析/6免疫治疗"))
test<-rt
test<-na.omit(test)
res.cut <- surv_cutpoint(test, time = "time", event = "status",
                         variables = gene)
summary(res.cut)
#Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)
#Fit survival curves and visualize
fit <- survfit(Surv(time,status ) ~get(gene), data = res.cat)


pdf(file = "GSE32894最佳cutoff-DFS生存分析.pdf",width = 5.6,height = 5.34,onefile = F)

ggsurvplot(fit,
           data=res.cat,
           risk.table = T,
           risk.table.col = 'strata',
           risk.table.y.text = F,
           risk.table.title = 'Number at risk',
           pval = TRUE,
           #conf.int = TRUE,
           xlab = 'Time (years)',
           ggtheme = theme_light(),
           surv.median.line = 'hv',
           title=paste( gene,'Survival',sep = ' ' ))

dev.off()

#百分比图
res.cat$ID<-test$ID
res.cat$Progression<-test$Progression
res.cat$Group<-res.cat[,gene]
library(dplyr)
library(plyr)
library(ggplot2)
library(ggpubr)
dat = dplyr::count(res.cat,Group,Progression)
dat = dat %>% group_by(Group) %>% 
  dplyr::summarise(Progression = Progression,n = n/sum(n))
dat$Progression = factor(dat$Progression,levels = c("no","yes"))


#计算百分率
dat=ddply(dat, .(Group), transform, percent = n/sum(n) * 100)
#百分比位置
dat=ddply(dat, .(Group), transform, pos = (cumsum(n) - 0.5 * n))
dat$label=paste0(sprintf("%.0f", dat$percent), "%")

bioCol=c("#f87669","#2fa1dd")
p=ggplot(dat, aes(x = factor(Group), y = percent, fill = Progression)) +
  geom_bar(position = position_stack(), stat = "identity", width = .7) +
  scale_fill_manual(values=bioCol)+
  xlab("Group")+ ylab("Percent weight")+  guides(fill=guide_legend(title="Progression"))+
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  #coord_flip()+
  theme_bw()
pdf(file="GSE32894_barplot.pdf", width=4, height=5)
print(p)
dev.off()



test$Progression<-factor(test$Progression,levels = c('no','yes'))
my_comparisons <- list( c("no", "yes") )
#耐受非耐受差异表达
p <- ggboxplot(test, x = "Progression", y = gene,
               fill = "Progression",legend=F,palette =c("#E7B800", "#00AFBB"),bxp.errorbar=T)+
  theme(legend.position='none')+
  ylab(label = gene)+
  stat_compare_means(comparisons = my_comparisons,method = 't.test',aes(label = ..p.signif..))
ggsave(p,filename = paste("GSE32894",gene,'差异表达.pdf',sep=' '),width = 5.6,height = 4.22)



#预测反应性ROC曲线
library("pROC")
df<-dplyr::select(rt,Progression,gene,PDCD1,CTLA4,CD274)
#定义足够多的颜色，后面画线时从这里选颜色
mycol <- c("slateblue","seagreen3","dodgerblue","firebrick1","lightgoldenrod","magenta","orange2")

if(TRUE){
  pdf("GSE32894_ROC.pdf",height=6,width=6)
  auc.out <- c()
  
  #先画第一条线，此处是miRNA1
  x <- plot.roc(df[,1],df[,2],ylim=c(0,1),xlim=c(1,0),
                smooth=F, #绘制平滑曲线
                ci=TRUE, 
                main="",
                #print.thres="best", #把阈值写在图上，其sensitivity+ specificity之和最大
                col=mycol[2],#线的颜色
                lwd=2, #线的粗细
                legacy.axes=T)#采用大多数paper的画法，横坐标是“1-specificity”，从0到1
  
  ci.lower <- round(as.numeric(x$ci[1]),3) #置信区间下限
  ci.upper <- round(as.numeric(x$ci[3]),3) #置信区间上限
  
  auc.ci <- c(colnames(df)[2],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
  auc.out <- rbind(auc.out,auc.ci)
  
  
  #再用循环画第二条和后面更多条曲线
  for (i in 3:ncol(df)){
    x <- plot.roc(df[,1],df[,i],
                  add=T, #向前面画的图里添加
                  smooth=F,
                  ci=TRUE,
                  col=mycol[i],
                  lwd=2,
                  legacy.axes=T)
    
    ci.lower <- round(as.numeric(x$ci[1]),3)
    ci.upper <- round(as.numeric(x$ci[3]),3)
    
    auc.ci <- c(colnames(df)[i],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
    auc.out <- rbind(auc.out,auc.ci)
  }
  
  
  #对比多条曲线
  #在参数`method=`后面，有三种方法可选“delong”, “bootstrap”或“venkatraman”，计算p值
  
  auc.out <- as.data.frame(auc.out)
  colnames(auc.out) <- c("Name","AUC","AUC CI")
  
  #绘制图例
  legend.name <- paste(colnames(df)[2:length(df)],"AUC",auc.out$AUC,sep=" ")
  legend("bottomright", 
         legend=legend.name,
         col = mycol[2:length(df)],
         lwd = 2,
         bty="n")
  dev.off()
}


###################################################################################
rm(list = ls())
gene="Score"
setwd("Y:\\lcd代码\\TCGA随意做\\TCGA任意做\\进阶版\\进阶\\原始数据\\免疫治疗数据")
genelist_input <- fread(file="GSE13507.csv", header = T, sep=',', data.table = F)
genename <- as.character(genelist_input[,1]) #提取第一列基因名
#x - 基因组注释R包
#keys - 需要转换的基因列表
#keytype - 基因名类型
#columns - 希望返回的基因名类型数据，如返回NCBI基因ID使用 ENTREZID
gene_map <- biomaRt::select(org.Hs.eg.db, keys=genename, keytype="ENTREZID", columns=c("SYMBOL"))
colnames(gene_map)[1]<-"ID"
colnames(genelist_input)[1]<-"ID"
genelist_input$ID<-as.character(genelist_input$ID)
immudata<-inner_join(gene_map,genelist_input,by="ID")
immudata<-immudata[,-1]
immudata<-na.omit(immudata)
rownames(immudata)<-immudata$SYMBOL
immudata<-immudata[,-1]
immudata<-as.data.frame(t(immudata))

CCC=read.table("uniSigGeneExp.txt",header = T,sep = "\t")
CCC=CCC$id
CCC=intersect(colnames(immudata),CCC)
AAA=immudata[,CCC]
#PCA分析
pca=prcomp(AAA, scale=TRUE)
value=predict(pca)
score=value[,1]+value[,2]
score=as.data.frame(score)
score$ID=rownames(score)
immudata$ID=row.names(immudata)
immudata=inner_join(immudata,score,by="ID")
colnames(immudata)[which(colnames(immudata)=="score")]=gene
immudata<-immudata %>%
  dplyr::select(ID,gene,everything())


clindata<-read.csv("GSE13507_clin.csv",header = T,sep = ",")
colnames(clindata)[1]<-"ID"
AAA<-intersect(clindata$ID,immudata$ID)
immudata<-dplyr::filter(immudata,ID %in% AAA)
clindata<-dplyr::filter(clindata,ID %in% AAA)
rt<-inner_join(clindata,immudata,by="ID")

rt$OS.time<-rt$OS.time/12
rt$OS<-as.numeric(rt$OS)

#生存分析
setwd(paste0("Y:\\lcd代码\\TCGA随意做\\TCGA任意做\\进阶版\\进阶/",gene,"/5免疫浸润分析/6免疫治疗"))
test<-rt
test<-na.omit(test)
res.cut <- surv_cutpoint(test, time = "OS.time", event = "OS",
                         variables = gene)
summary(res.cut)
#Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)
#Fit survival curves and visualize
fit <- survfit(Surv(OS.time,OS ) ~get(gene), data = res.cat)


pdf(file = "GSE13507最佳cutoff-OS生存分析.pdf",width = 5.6,height = 5.34,onefile = F)

ggsurvplot(fit,
           data=res.cat,
           risk.table = T,
           risk.table.col = 'strata',
           risk.table.y.text = F,
           risk.table.title = 'Number at risk',
           pval = TRUE,
           #conf.int = TRUE,
           xlab = 'Time (years)',
           ggtheme = theme_light(),
           surv.median.line = 'hv',
           title=paste( gene,'Survival',sep = ' ' ))

dev.off()



#百分比图
res.cat$ID<-test$ID
res.cat$Progression<-test$Progression
res.cat$Group<-res.cat[,gene]
library(dplyr)
library(plyr)
library(ggplot2)
library(ggpubr)
dat = dplyr::count(res.cat,Group,Progression)
dat = dat %>% group_by(Group) %>% 
  dplyr::summarise(Progression = Progression,n = n/sum(n))
dat$Progression = factor(dat$Progression,levels = c("No","Yes"))


#计算百分率
dat=ddply(dat, .(Group), transform, percent = n/sum(n) * 100)
#百分比位置
dat=ddply(dat, .(Group), transform, pos = (cumsum(n) - 0.5 * n))
dat$label=paste0(sprintf("%.0f", dat$percent), "%")

bioCol=c("#f87669","#2fa1dd")
p=ggplot(dat, aes(x = factor(Group), y = percent, fill = Progression)) +
  geom_bar(position = position_stack(), stat = "identity", width = .7) +
  scale_fill_manual(values=bioCol)+
  xlab("Group")+ ylab("Percent weight")+  guides(fill=guide_legend(title="Progression"))+
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  #coord_flip()+
  theme_bw()
pdf(file="GSE13507_barplot.pdf", width=4, height=5)
print(p)
dev.off()

test$Progression<-factor(test$Progression,levels = c('No','Yes'))
my_comparisons <- list( c("No", "Yes") )
#耐受非耐受差异表达
p <- ggboxplot(test, x = "Progression", y = gene,
               fill = "Progression",legend=F,palette =c("#E7B800", "#00AFBB"),bxp.errorbar=T)+
  theme(legend.position='none')+
  ylab(label = gene)+
  stat_compare_means(comparisons = my_comparisons,method = 't.test',aes(label = ..p.signif..))
ggsave(p,filename = paste("GSE13507",gene,'差异表达.pdf',sep=' '),width = 5.6,height = 4.22)



#预测反应性ROC曲线
library("pROC")
df<-dplyr::select(rt,Progression,gene,PDCD1,CTLA4,CD274)
df$Progression[df$Progression=="Yes"]=1
df$Progression[df$Progression=="No"]=0
#定义足够多的颜色，后面画线时从这里选颜色
mycol <- c("slateblue","seagreen3","dodgerblue","firebrick1","lightgoldenrod","magenta","orange2")

if(TRUE){
  pdf("GSE13507_ROC.pdf",height=6,width=6)
  auc.out <- c()
  
  #先画第一条线，此处是miRNA1
  x <- plot.roc(df[,1],df[,2],ylim=c(0,1),xlim=c(1,0),
                smooth=F, #绘制平滑曲线
                ci=TRUE, 
                main="",
                #print.thres="best", #把阈值写在图上，其sensitivity+ specificity之和最大
                col=mycol[2],#线的颜色
                lwd=2, #线的粗细
                legacy.axes=T)#采用大多数paper的画法，横坐标是“1-specificity”，从0到1
  
  ci.lower <- round(as.numeric(x$ci[1]),3) #置信区间下限
  ci.upper <- round(as.numeric(x$ci[3]),3) #置信区间上限
  
  auc.ci <- c(colnames(df)[2],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
  auc.out <- rbind(auc.out,auc.ci)
  
  
  #再用循环画第二条和后面更多条曲线
  for (i in 3:ncol(df)){
    x <- plot.roc(df[,1],df[,i],
                  add=T, #向前面画的图里添加
                  smooth=F,
                  ci=TRUE,
                  col=mycol[i],
                  lwd=2,
                  legacy.axes=T)
    
    ci.lower <- round(as.numeric(x$ci[1]),3)
    ci.upper <- round(as.numeric(x$ci[3]),3)
    
    auc.ci <- c(colnames(df)[i],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
    auc.out <- rbind(auc.out,auc.ci)
  }
  
  
  auc.out <- as.data.frame(auc.out)
  colnames(auc.out) <- c("Name","AUC","AUC CI")
  #输出p value到文件
  
  #绘制图例
  legend.name <- paste(colnames(df)[2:length(df)],"AUC",auc.out$AUC,sep=" ")
  legend("bottomright", 
         legend=legend.name,
         col = mycol[2:length(df)],
         lwd = 2,
         bty="n")
  dev.off()
}

###################################################################################
rm(list = ls())
gene="Score"
setwd("Y:\\lcd代码\\TCGA随意做\\TCGA任意做\\进阶版\\进阶\\原始数据\\免疫治疗数据")
genelist_input <- fread(file="GSE35640.csv", header = T, sep=',', data.table = F)
genename <- as.character(genelist_input[,1]) #提取第一列基因名
#x - 基因组注释R包
#keys - 需要转换的基因列表
#keytype - 基因名类型
#columns - 希望返回的基因名类型数据，如返回NCBI基因ID使用 ENTREZID
gene_map <- biomaRt::select(org.Hs.eg.db, keys=genename, keytype="ENTREZID", columns=c("SYMBOL"))
colnames(gene_map)[1]<-"ID"
colnames(genelist_input)[1]<-"ID"
genelist_input$ID<-as.character(genelist_input$ID)
immudata<-inner_join(gene_map,genelist_input,by="ID")
immudata<-immudata[,-1]
immudata<-na.omit(immudata)
rownames(immudata)<-immudata$SYMBOL
immudata<-immudata[,-1]
immudata<-as.data.frame(t(immudata))

CCC=read.table("uniSigGeneExp.txt",header = T,sep = "\t")
CCC=CCC$id
CCC=intersect(colnames(immudata),CCC)
AAA=immudata[,CCC]
#PCA分析
pca=prcomp(AAA, scale=TRUE)
value=predict(pca)
score=value[,1]+value[,2]
score=as.data.frame(score)
score$ID=rownames(score)
immudata$ID=row.names(immudata)
immudata=inner_join(immudata,score,by="ID")
colnames(immudata)[which(colnames(immudata)=="score")]=gene
immudata<-immudata %>%
  dplyr::select(ID,gene,everything())

clindata<-read.csv("GSE35640_clin.csv",header = T,sep = ",")
colnames(clindata)[1]<-"ID"
AAA<-intersect(clindata$ID,immudata$ID)
immudata<-dplyr::filter(immudata,ID %in% AAA)
clindata<-dplyr::filter(clindata,ID %in% AAA)
rt<-inner_join(clindata,immudata,by="ID")

test<-rt
test$Response<-factor(test$Response,levels = c('non-responder','responder'))
my_comparisons <- list( c("non-responder", "responder") )
#耐受非耐受差异表达
setwd(paste0("Y:\\lcd代码\\TCGA随意做\\TCGA任意做\\进阶版\\进阶/",gene,"/5免疫浸润分析/6免疫治疗"))
p <- ggboxplot(test, x = "Response", y = gene,
               fill = "Response",legend=F,palette =c("#E7B800", "#00AFBB"),bxp.errorbar=T)+
  theme(legend.position='none')+
  ylab(label = gene)+
  stat_compare_means(comparisons = my_comparisons,method = 't.test',aes(label = ..p.signif..))
ggsave(p,filename = paste("GSE35640",gene,'差异表达.pdf',sep=' '),width = 5.6,height = 4.22)



#预测反应性ROC曲线
library("pROC")
df<-dplyr::select(rt,Response,gene,PDCD1,CTLA4,CD274)
#定义足够多的颜色，后面画线时从这里选颜色
mycol <- c("slateblue","seagreen3","dodgerblue","firebrick1","lightgoldenrod","magenta","orange2")

if(TRUE){
  pdf("GSE35640_ROC.pdf",height=6,width=6)
  auc.out <- c()
  
  #先画第一条线，此处是miRNA1
  x <- plot.roc(df[,1],df[,2],ylim=c(0,1),xlim=c(1,0),
                smooth=F, #绘制平滑曲线
                ci=TRUE, 
                main="",
                #print.thres="best", #把阈值写在图上，其sensitivity+ specificity之和最大
                col=mycol[2],#线的颜色
                lwd=2, #线的粗细
                legacy.axes=T)#采用大多数paper的画法，横坐标是“1-specificity”，从0到1
  
  ci.lower <- round(as.numeric(x$ci[1]),3) #置信区间下限
  ci.upper <- round(as.numeric(x$ci[3]),3) #置信区间上限
  
  auc.ci <- c(colnames(df)[2],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
  auc.out <- rbind(auc.out,auc.ci)
  
  
  #再用循环画第二条和后面更多条曲线
  for (i in 3:ncol(df)){
    x <- plot.roc(df[,1],df[,i],
                  add=T, #向前面画的图里添加
                  smooth=F,
                  ci=TRUE,
                  col=mycol[i],
                  lwd=2,
                  legacy.axes=T)
    
    ci.lower <- round(as.numeric(x$ci[1]),3)
    ci.upper <- round(as.numeric(x$ci[3]),3)
    
    auc.ci <- c(colnames(df)[i],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
    auc.out <- rbind(auc.out,auc.ci)
  }
  
  
  #对比多条曲线
  #在参数`method=`后面，有三种方法可选“delong”, “bootstrap”或“venkatraman”，计算p值
  auc.out <- as.data.frame(auc.out)
  colnames(auc.out) <- c("Name","AUC","AUC CI")
  #输出p value到文件
  
  #绘制图例
  legend.name <- paste(colnames(df)[2:length(df)],"AUC",auc.out$AUC,sep=" ")
  legend("bottomright", 
         legend=legend.name,
         col = mycol[2:length(df)],
         lwd = 2,
         bty="n")
  dev.off()
}

###################################################################################
rm(list = ls())
gene="Score"
setwd("Y:\\lcd代码\\TCGA随意做\\TCGA任意做\\进阶版\\进阶\\原始数据\\免疫治疗数据")
genelist_input <- fread(file="GSE61676.csv", header = T, sep=',', data.table = F)
immudata<-genelist_input
rownames(immudata)<-immudata$ID
immudata<-immudata[,-1]
immudata<-as.matrix(immudata)
immudata<-as.data.frame(t(immudata))

CCC=read.table("uniSigGeneExp.txt",header = T,sep = "\t")
CCC=CCC$id
CCC=intersect(colnames(immudata),CCC)
AAA=immudata[,CCC]
#PCA分析
pca=prcomp(AAA, scale=TRUE)
value=predict(pca)
score=value[,1]+value[,2]
score=as.data.frame(score)
score$ID=rownames(score)
immudata$ID=row.names(immudata)
immudata=inner_join(immudata,score,by="ID")
colnames(immudata)[which(colnames(immudata)=="score")]=gene
immudata<-immudata %>%
  dplyr::select(ID,gene,everything())


clindata<-read.csv("GSE61676_clin.csv",header = T,sep = ",")
colnames(clindata)[1]<-"ID"
AAA<-intersect(clindata$ID,immudata$ID)
immudata<-dplyr::filter(immudata,ID %in% AAA)
clindata<-dplyr::filter(clindata,ID %in% AAA)
rt<-inner_join(clindata,immudata,by="ID")

rt$time<-rt$time/12
rt$status<-as.numeric(rt$status)

#生存分析
setwd(paste0("Y:\\lcd代码\\TCGA随意做\\TCGA任意做\\进阶版\\进阶/",gene,"/5免疫浸润分析/6免疫治疗"))
test<-rt
res.cut <- surv_cutpoint(test, time = "time", event = "status",
                         variables = gene)
summary(res.cut)
#Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)
#Fit survival curves and visualize
fit <- survfit(Surv(time,status ) ~get(gene), data = res.cat)


pdf(file = "GSE61676最佳cutoff生存分析.pdf",width = 5.6,height = 5.34,onefile = F)

ggsurvplot(fit,
           data=res.cat,
           risk.table = T,
           risk.table.col = 'strata',
           risk.table.y.text = F,
           risk.table.title = 'Number at risk',
           pval = TRUE,
           #conf.int = TRUE,
           xlab = 'Time (years)',
           ggtheme = theme_light(),
           surv.median.line = 'hv',
           title=paste( gene,'Survival',sep = ' ' ))

dev.off()


#百分比图
res.cat$ID<-test$ID
res.cat$Response<-test$Response
res.cat$Group<-res.cat[,gene]
library(dplyr)
library(plyr)
library(ggplot2)
library(ggpubr)
dat = dplyr::count(res.cat,Group,Response)
dat = dat %>% group_by(Group) %>% 
  dplyr::summarise(Response = Response,n = n/sum(n))
dat$Response = factor(dat$Response,levels = c("Non-responder","Responder"))


#计算百分率
dat=ddply(dat, .(Group), transform, percent = n/sum(n) * 100)
#百分比位置
dat=ddply(dat, .(Group), transform, pos = (cumsum(n) - 0.5 * n))
dat$label=paste0(sprintf("%.0f", dat$percent), "%")

bioCol=c("#f87669","#2fa1dd")
p=ggplot(dat, aes(x = factor(Group), y = percent, fill = Response)) +
  geom_bar(position = position_stack(), stat = "identity", width = .7) +
  scale_fill_manual(values=bioCol)+
  xlab("Group")+ ylab("Percent weight")+  guides(fill=guide_legend(title="Response"))+
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  #coord_flip()+
  theme_bw()
pdf(file="GSE61676_barplot.pdf", width=4, height=5)
print(p)
dev.off()

test$Response<-factor(test$Response,levels = c('Non-responder','Responder'))
my_comparisons <- list( c("Non-responder", "Responder") )
#耐受非耐受差异表达
p <- ggboxplot(test, x = "Response", y = gene,
               fill = "Response",legend=F,palette =c("#E7B800", "#00AFBB"),bxp.errorbar=T)+
  theme(legend.position='none')+
  ylab(label = gene)+
  stat_compare_means(comparisons = my_comparisons,method = 't.test',aes(label = ..p.signif..))
ggsave(p,filename = paste("GSE61676",gene,'差异表达.pdf',sep=' '),width = 5.6,height = 4.22)



#预测反应性ROC曲线
library("pROC")
df<-dplyr::select(rt,Response,gene,CTLA4,CD274)
df$Response[df$Response=="Responder"]=1
df$Response[df$Response=="Non-responder"]=0
#定义足够多的颜色，后面画线时从这里选颜色
mycol <- c("slateblue","seagreen3","dodgerblue","firebrick1","lightgoldenrod","magenta","orange2")

if(TRUE){
  pdf("GSE61676_ROC.pdf",height=6,width=6)
  auc.out <- c()
  
  #先画第一条线，此处是miRNA1
  x <- plot.roc(df[,1],df[,2],ylim=c(0,1),xlim=c(1,0),
                smooth=F, #绘制平滑曲线
                ci=TRUE, 
                main="",
                #print.thres="best", #把阈值写在图上，其sensitivity+ specificity之和最大
                col=mycol[2],#线的颜色
                lwd=2, #线的粗细
                legacy.axes=T)#采用大多数paper的画法，横坐标是“1-specificity”，从0到1
  
  ci.lower <- round(as.numeric(x$ci[1]),3) #置信区间下限
  ci.upper <- round(as.numeric(x$ci[3]),3) #置信区间上限
  
  auc.ci <- c(colnames(df)[2],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
  auc.out <- rbind(auc.out,auc.ci)
  
  
  #再用循环画第二条和后面更多条曲线
  for (i in 3:ncol(df)){
    x <- plot.roc(df[,1],df[,i],
                  add=T, #向前面画的图里添加
                  smooth=F,
                  ci=TRUE,
                  col=mycol[i],
                  lwd=2,
                  legacy.axes=T)
    
    ci.lower <- round(as.numeric(x$ci[1]),3)
    ci.upper <- round(as.numeric(x$ci[3]),3)
    
    auc.ci <- c(colnames(df)[i],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
    auc.out <- rbind(auc.out,auc.ci)
  }
  
  auc.out <- as.data.frame(auc.out)
  colnames(auc.out) <- c("Name","AUC","AUC CI")
  #输出p value到文件
  
  #绘制图例
  legend.name <- paste(colnames(df)[2:length(df)],"AUC",auc.out$AUC,sep=" ")
  legend("bottomright", 
         legend=legend.name,
         col = mycol[2:length(df)],
         lwd = 2,
         bty="n")
  dev.off()
}


###################################################################################
rm(list = ls())
gene="Score"
setwd("Y:\\lcd代码\\TCGA随意做\\TCGA任意做\\进阶版\\进阶\\原始数据\\免疫治疗数据")
genelist_input <- fread(file="IMvigor210CoreBiologies.csv", header = T, sep=',', data.table = F)
immudata<-genelist_input


#构建一个布尔向量，索引
index<-duplicated(immudata$symbol)
index
#筛选数据
immudata<-immudata[!index,]  #选中了非重复的数据
immudata<-na.omit(immudata)

rownames(immudata)<-immudata$symbol
immudata<-immudata[,-1]
immudata<-log2(immudata+1)
immudata<-as.matrix(immudata)
immudata<-as.data.frame(t(immudata))

CCC=read.table("uniSigGeneExp.txt",header = T,sep = "\t")
CCC=CCC$id
CCC=intersect(colnames(immudata),CCC)
AAA=immudata[,CCC]
#PCA分析
pca=prcomp(AAA, scale=TRUE)
value=predict(pca)
score=value[,1]+value[,2]
score=as.data.frame(score)
score$ID=rownames(score)
immudata$ID=row.names(immudata)
immudata=inner_join(immudata,score,by="ID")
colnames(immudata)[which(colnames(immudata)=="score")]=gene
immudata<-immudata %>%
  dplyr::select(ID,gene,everything())


clindata<-read.csv("IMvigor210CoreBiologies_clin.csv",header = T,sep = ",")
colnames(clindata)[1]<-"ID"
AAA<-intersect(clindata$ID,immudata$ID)
immudata<-dplyr::filter(immudata,ID %in% AAA)
clindata<-dplyr::filter(clindata,ID %in% AAA)
rt<-inner_join(clindata,immudata,by="ID")

rt$futime<-rt$futime/12
rt$fustat<-as.numeric(rt$fustat)

#生存分析
setwd(paste0("Y:\\lcd代码\\TCGA随意做\\TCGA任意做\\进阶版\\进阶/",gene,"/5免疫浸润分析/6免疫治疗"))
test<-rt

res.cut <- surv_cutpoint(test, time = "futime", event = "fustat",
                         variables = gene)
summary(res.cut)
#Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)
#Fit survival curves and visualize
fit <- survfit(Surv(futime,fustat ) ~get(gene), data = res.cat)


pdf(file = "IMvigor210CoreBiologies最佳cutoff生存分析.pdf",width = 5.6,height = 5.34,onefile = F)

ggsurvplot(fit,
           data=res.cat,
           risk.table = T,
           risk.table.col = 'strata',
           risk.table.y.text = F,
           risk.table.title = 'Number at risk',
           pval = TRUE,
           #conf.int = TRUE,
           xlab = 'Time (years)',
           ggtheme = theme_light(),
           surv.median.line = 'hv',
           title=paste( gene,'Survival',sep = ' ' ))

dev.off()

#百分比图
res.cat$ID<-test$ID
res.cat$Response<-test$Response
res.cat$Group<-res.cat[,gene]
res.cat<-na.omit(res.cat)
library(dplyr)
library(plyr)
library(ggplot2)
library(ggpubr)
dat = dplyr::count(res.cat,Group,Response)
dat = dat %>% group_by(Group) %>% 
  dplyr::summarise(Response = Response,n = n/sum(n))
dat$Response = factor(dat$Response,levels = c("CR/PR","SD/PD"))


#计算百分率
dat=ddply(dat, .(Group), transform, percent = n/sum(n) * 100)
#百分比位置
dat=ddply(dat, .(Group), transform, pos = (cumsum(n) - 0.5 * n))
dat$label=paste0(sprintf("%.0f", dat$percent), "%")

bioCol=c("#f87669","#2fa1dd")
p=ggplot(dat, aes(x = factor(Group), y = percent, fill = Response)) +
  geom_bar(position = position_stack(), stat = "identity", width = .7) +
  scale_fill_manual(values=bioCol)+
  xlab("Group")+ ylab("Percent weight")+  guides(fill=guide_legend(title="Response"))+
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  #coord_flip()+
  theme_bw()
pdf(file="IMvigor210CoreBiologies_barplot.pdf", width=4, height=5)
print(p)
dev.off()

test$Response<-factor(test$Response,levels = c('CR/PR','SD/PD'))
my_comparisons <- list( c("CR/PR", "SD/PD") )
#耐受非耐受差异表达
test<-na.omit(test)
p <- ggboxplot(test, x = "Response", y = gene,
               fill = "Response",legend=F,palette =c("#E7B800", "#00AFBB"),bxp.errorbar=T)+
  theme(legend.position='none')+
  ylab(label = gene)+
  stat_compare_means(comparisons = my_comparisons,method = 't.test',aes(label = ..p.signif..))
ggsave(p,filename = paste("IMvigor210CoreBiologies",gene,'差异表达.pdf',sep=' '),width = 5.6,height = 4.22)



#预测反应性ROC曲线
library("pROC")
df<-dplyr::select(rt,Response,gene,CTLA4,CD274,PDCD1)
df$Response[df$Response=="SD/PD"]=1
df$Response[df$Response=="CR/PR"]=0
#定义足够多的颜色，后面画线时从这里选颜色
mycol <- c("slateblue","seagreen3","dodgerblue","firebrick1","lightgoldenrod","magenta","orange2")

if(TRUE){
  pdf("IMvigor210CoreBiologies_ROC.pdf",height=6,width=6)
  auc.out <- c()
  
  #先画第一条线，此处是miRNA1
  x <- plot.roc(df[,1],df[,2],ylim=c(0,1),xlim=c(1,0),
                smooth=F, #绘制平滑曲线
                ci=TRUE, 
                main="",
                #print.thres="best", #把阈值写在图上，其sensitivity+ specificity之和最大
                col=mycol[2],#线的颜色
                lwd=2, #线的粗细
                legacy.axes=T)#采用大多数paper的画法，横坐标是“1-specificity”，从0到1
  
  ci.lower <- round(as.numeric(x$ci[1]),3) #置信区间下限
  ci.upper <- round(as.numeric(x$ci[3]),3) #置信区间上限
  
  auc.ci <- c(colnames(df)[2],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
  auc.out <- rbind(auc.out,auc.ci)
  
  
  #再用循环画第二条和后面更多条曲线
  for (i in 3:ncol(df)){
    x <- plot.roc(df[,1],df[,i],
                  add=T, #向前面画的图里添加
                  smooth=F,
                  ci=TRUE,
                  col=mycol[i],
                  lwd=2,
                  legacy.axes=T)
    
    ci.lower <- round(as.numeric(x$ci[1]),3)
    ci.upper <- round(as.numeric(x$ci[3]),3)
    
    auc.ci <- c(colnames(df)[i],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
    auc.out <- rbind(auc.out,auc.ci)
  }
  
  
  #对比多条曲线
  #在参数`method=`后面，有三种方法可选“delong”, “bootstrap”或“venkatraman”，计算p值
  # p.out <- c()
  # for (i in 2:(ncol(df)-1)){
  #   for (j in (i+1):ncol(df)){
  #     p <- roc.test(df[,1],df[,i],df[,j], method="bootstrap")
  #     p.tmp <- c(colnames(df)[i],colnames(df)[j],p$p.value)
  #     p.out <- rbind(p.out,p.tmp)
  #   }
  # }
  # p.out <- as.data.frame(p.out)
  # colnames(p.out) <- c("ROC1","ROC2","p.value")
  auc.out <- as.data.frame(auc.out)
  colnames(auc.out) <- c("Name","AUC","AUC CI")
  #输出p value到文件
  
  #绘制图例
  legend.name <- paste(colnames(df)[2:length(df)],"AUC",auc.out$AUC,sep=" ")
  legend("bottomright", 
         legend=legend.name,
         col = mycol[2:length(df)],
         lwd = 2,
         bty="n")
  dev.off()
}

###################################################################################
rm(list = ls())
gene="Score"
setwd("Y:\\lcd代码\\TCGA随意做\\TCGA任意做\\进阶版\\进阶\\原始数据\\免疫治疗数据")
genelist_input <- fread(file="GSE78220.csv", header = T, sep=',', data.table = F)
immudata<-genelist_input


#构建一个布尔向量，索引
index<-duplicated(immudata$Gene)
index
#筛选数据
immudata<-immudata[!index,]  #选中了非重复的数据
immudata<-na.omit(immudata)

rownames(immudata)<-immudata$Gene
immudata<-immudata[,-1]
immudata<-as.matrix(immudata)
immudata<-as.data.frame(t(immudata))

CCC=read.table("uniSigGeneExp.txt",header = T,sep = "\t")
CCC=CCC$id
CCC=intersect(colnames(immudata),CCC)
AAA=immudata[,CCC]
#PCA分析
pca=prcomp(AAA, scale=TRUE)
value=predict(pca)
score=value[,1]+value[,2]
score=as.data.frame(score)
score$ID=rownames(score)
immudata$ID=row.names(immudata)
immudata=inner_join(immudata,score,by="ID")
colnames(immudata)[which(colnames(immudata)=="score")]=gene
immudata<-immudata %>%
  dplyr::select(ID,gene,everything())


clindata<-read.csv("GSE78220_clin.csv",header = T,sep = ",")
colnames(clindata)[1]<-"ID"
AAA<-intersect(clindata$ID,immudata$ID)
immudata<-dplyr::filter(immudata,ID %in% AAA)
clindata<-dplyr::filter(clindata,ID %in% AAA)
rt<-inner_join(clindata,immudata,by="ID")

rt$OS<-rt$OS/365
rt$Status<-as.numeric(rt$Status)

#生存分析
setwd(paste0("Y:\\lcd代码\\TCGA随意做\\TCGA任意做\\进阶版\\进阶/",gene,"/5免疫浸润分析/6免疫治疗"))
test<-rt
test<-na.omit(test)
res.cut <- surv_cutpoint(test, time = "OS", event = "Status",
                         variables = gene)
summary(res.cut)
#Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)
#Fit survival curves and visualize
fit <- survfit(Surv(OS,Status ) ~get(gene), data = res.cat)


pdf(file = "GSE78220最佳cutoff生存分析.pdf",width = 5.6,height = 5.34,onefile = F)

ggsurvplot(fit,
           data=res.cat,
           risk.table = T,
           risk.table.col = 'strata',
           risk.table.y.text = F,
           risk.table.title = 'Number at risk',
           pval = TRUE,
           #conf.int = TRUE,
           xlab = 'Time (years)',
           ggtheme = theme_light(),
           surv.median.line = 'hv',
           title=paste( gene,'Survival',sep = ' ' ))

dev.off()
#百分比图
res.cat$ID<-test$ID
res.cat$Response<-test$Response
res.cat$Group<-res.cat[,gene]
library(dplyr)
library(plyr)
library(ggplot2)
library(ggpubr)
dat = dplyr::count(res.cat,Group,Response)
dat = dat %>% group_by(Group) %>% 
  dplyr::summarise(Response = Response,n = n/sum(n))
dat$Response = factor(dat$Response,levels = c("PD/SD","PR/CR"))


#计算百分率
dat=ddply(dat, .(Group), transform, percent = n/sum(n) * 100)
#百分比位置
dat=ddply(dat, .(Group), transform, pos = (cumsum(n) - 0.5 * n))
dat$label=paste0(sprintf("%.0f", dat$percent), "%")

bioCol=c("#f87669","#2fa1dd")
p=ggplot(dat, aes(x = factor(Group), y = percent, fill = Response)) +
  geom_bar(position = position_stack(), stat = "identity", width = .7) +
  scale_fill_manual(values=bioCol)+
  xlab("Group")+ ylab("Percent weight")+  guides(fill=guide_legend(title="Response"))+
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  #coord_flip()+
  theme_bw()
pdf(file="GSE78220_barplot.pdf", width=4, height=5)
print(p)
dev.off()


test$Response<-factor(test$Response,levels = c('PR/CR','PD/SD'))
my_comparisons <- list( c("PR/CR", "PD/SD") )
#耐受非耐受差异表达
test<-na.omit(test)
p <- ggboxplot(test, x = "Response", y = gene,
               fill = "Response",legend=F,palette =c("#E7B800", "#00AFBB"),bxp.errorbar=T)+
  theme(legend.position='none')+
  ylab(label = gene)+
  stat_compare_means(comparisons = my_comparisons,method = 't.test',aes(label = ..p.signif..))
ggsave(p,filename = paste("GSE78220",gene,'差异表达.pdf',sep=' '),width = 5.6,height = 4.22)


#预测反应性ROC曲线
library("pROC")
df<-dplyr::select(rt,Response,gene,CTLA4,CD274,PDCD1)
df$Response[df$Response=="SD/PD"]=1
df$Response[df$Response=="CR/PR"]=0
#定义足够多的颜色，后面画线时从这里选颜色
mycol <- c("slateblue","seagreen3","dodgerblue","firebrick1","lightgoldenrod","magenta","orange2")

if(TRUE){
  pdf("GSE78220_ROC.pdf",height=6,width=6)
  auc.out <- c()
  
  #先画第一条线，此处是miRNA1
  x <- plot.roc(df[,1],df[,2],ylim=c(0,1),xlim=c(1,0),
                smooth=F, #绘制平滑曲线
                ci=TRUE, 
                main="",
                #print.thres="best", #把阈值写在图上，其sensitivity+ specificity之和最大
                col=mycol[2],#线的颜色
                lwd=2, #线的粗细
                legacy.axes=T)#采用大多数paper的画法，横坐标是“1-specificity”，从0到1
  
  ci.lower <- round(as.numeric(x$ci[1]),3) #置信区间下限
  ci.upper <- round(as.numeric(x$ci[3]),3) #置信区间上限
  
  auc.ci <- c(colnames(df)[2],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
  auc.out <- rbind(auc.out,auc.ci)
  
  
  #再用循环画第二条和后面更多条曲线
  for (i in 3:ncol(df)){
    x <- plot.roc(df[,1],df[,i],
                  add=T, #向前面画的图里添加
                  smooth=F,
                  ci=TRUE,
                  col=mycol[i],
                  lwd=2,
                  legacy.axes=T)
    
    ci.lower <- round(as.numeric(x$ci[1]),3)
    ci.upper <- round(as.numeric(x$ci[3]),3)
    
    auc.ci <- c(colnames(df)[i],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
    auc.out <- rbind(auc.out,auc.ci)
  }
  
  
  # #对比多条曲线
  # #在参数`method=`后面，有三种方法可选“delong”, “bootstrap”或“venkatraman”，计算p值
  # p.out <- c()
  # for (i in 2:(ncol(df)-1)){
  #   for (j in (i+1):ncol(df)){
  #     p <- roc.test(df[,1],df[,i],df[,j], method="bootstrap")
  #     p.tmp <- c(colnames(df)[i],colnames(df)[j],p$p.value)
  #     p.out <- rbind(p.out,p.tmp)
  #   }
  # }
  # p.out <- as.data.frame(p.out)
  # colnames(p.out) <- c("ROC1","ROC2","p.value")
  auc.out <- as.data.frame(auc.out)
  colnames(auc.out) <- c("Name","AUC","AUC CI")
  #输出p value到文件
  
  #绘制图例
  legend.name <- paste(colnames(df)[2:length(df)],"AUC",auc.out$AUC,sep=" ")
  legend("bottomright", 
         legend=legend.name,
         col = mycol[2:length(df)],
         lwd = 2,
         bty="n")
  dev.off()
}

###################################################################################################
###################################################################################
rm(list = ls())
gene="Score"

setwd("Y:\\lcd代码\\TCGA随意做\\TCGA任意做\\进阶版\\进阶\\原始数据\\免疫治疗数据")
genelist_input <- fread(file="ICB.Riaz2017_Nivolumab_Melanoma_Naive.csv", header = T, sep=',', data.table = F)
genename <- as.character(genelist_input[,1]) #提取第一列基因名
#x - 基因组注释R包
#keys - 需要转换的基因列表
#keytype - 基因名类型
#columns - 希望返回的基因名类型数据，如返回NCBI基因ID使用 ENTREZID
gene_map <- biomaRt::select(org.Hs.eg.db, keys=genename, keytype="ENTREZID", columns=c("SYMBOL"))
colnames(gene_map)[1]<-"ID"
colnames(genelist_input)[1]<-"ID"
genelist_input$ID<-as.character(genelist_input$ID)
immudata<-inner_join(gene_map,genelist_input,by="ID")
immudata<-immudata[,-1]
immudata<-na.omit(immudata)
rownames(immudata)<-immudata$SYMBOL
immudata<-immudata[,-1]
immudata<-as.matrix(immudata)
immudata<-as.data.frame(t(immudata))

CCC=read.table("uniSigGeneExp.txt",header = T,sep = "\t")
CCC=CCC$id
CCC=intersect(colnames(immudata),CCC)
AAA=immudata[,CCC]
#PCA分析
pca=prcomp(AAA, scale=TRUE)
value=predict(pca)
score=value[,1]+value[,2]
score=as.data.frame(score)
score$ID=rownames(score)
immudata$ID=row.names(immudata)
immudata=inner_join(immudata,score,by="ID")
colnames(immudata)[which(colnames(immudata)=="score")]=gene
immudata<-immudata %>%
  dplyr::select(ID,gene,everything())
clindata<-read.csv("ICB.Riaz2017_Nivolumab_Melanoma_Naive_clin.csv",header = T,sep = ",")
colnames(clindata)[1]<-"ID"
AAA<-intersect(clindata$ID,immudata$ID)
immudata<-dplyr::filter(immudata,ID %in% AAA)
clindata<-dplyr::filter(clindata,ID %in% AAA)
rt<-inner_join(clindata,immudata,by="ID")

rt$OS<-rt$OS/365
rt$OS.Event<-as.numeric(rt$OS.Event)

#生存分析
setwd(paste0("Y:\\lcd代码\\TCGA随意做\\TCGA任意做\\进阶版\\进阶/",gene,"/5免疫浸润分析/6免疫治疗"))
test<-rt
test<-na.omit(test)
res.cut <- surv_cutpoint(test, time = "OS", event = "OS.Event",
                         variables = gene)
summary(res.cut)
#Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)
#Fit survival curves and visualize
fit <- survfit(Surv(OS,OS.Event ) ~get(gene), data = res.cat)


pdf(file = "ICB.Riaz2017_Nivolumab_Melanoma_Naive_OS生存分析.pdf",width = 5.6,height = 5.34,onefile = F)

ggsurvplot(fit,
           data=res.cat,
           risk.table = T,
           risk.table.col = 'strata',
           risk.table.y.text = F,
           risk.table.title = 'Number at risk',
           pval = TRUE,
           #conf.int = TRUE,
           xlab = 'Time (years)',
           ggtheme = theme_light(),
           surv.median.line = 'hv',
           title=paste( gene,'Survival',sep = ' ' ))

dev.off()

#百分比图
res.cat$ID<-test$ID
res.cat$Response<-test$Response
res.cat$Response[res.cat$Response==1]="Yes"
res.cat$Response[res.cat$Response==0]="No"
res.cat$Group<-res.cat[,gene]
library(dplyr)
library(plyr)
library(ggplot2)
library(ggpubr)
dat = dplyr::count(res.cat,Group,Response)
dat = dat %>% group_by(Group) %>% 
  dplyr::summarise(Response = Response,n = n/sum(n))
dat$Response = factor(dat$Response,levels = c("No","Yes"))


#计算百分率
dat=ddply(dat, .(Group), transform, percent = n/sum(n) * 100)
#百分比位置
dat=ddply(dat, .(Group), transform, pos = (cumsum(n) - 0.5 * n))
dat$label=paste0(sprintf("%.0f", dat$percent), "%")

bioCol=c("#f87669","#2fa1dd")
p=ggplot(dat, aes(x = factor(Group), y = percent, fill = Response)) +
  geom_bar(position = position_stack(), stat = "identity", width = .7) +
  scale_fill_manual(values=bioCol)+
  xlab("Group")+ ylab("Percent weight")+  guides(fill=guide_legend(title="Response"))+
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  #coord_flip()+
  theme_bw()
pdf(file="ICB.Riaz2017_Nivolumab_Melanoma_Naive_barplot.pdf", width=4, height=5)
print(p)
dev.off()

test<-rt
test<-na.omit(test)
test$Response[test$Response==1]="Yes"
test$Response[test$Response==0]="No"
test$Response<-factor(test$Response,levels = c('No','Yes'))
my_comparisons <- list( c("No", "Yes") )
#耐受非耐受差异表达
p <- ggboxplot(test, x = "Response", y = gene,
               fill = "Response",legend=F,palette =c("#E7B800", "#00AFBB"),bxp.errorbar=T)+
  theme(legend.position='none')+
  ylab(label = gene)+
  stat_compare_means(comparisons = my_comparisons,method = 't.test',aes(label = ..p.signif..))
ggsave(p,filename = paste("ICB.Riaz2017_Nivolumab_Melanoma_Naive",gene,'差异表达.pdf',sep=' '),width = 5.6,height = 4.22)



#预测反应性ROC曲线
library("pROC")
df<-dplyr::select(rt,Response,gene,PDCD1,CTLA4,CD274)
#定义足够多的颜色，后面画线时从这里选颜色
mycol <- c("slateblue","seagreen3","dodgerblue","firebrick1","lightgoldenrod","magenta","orange2")

if(TRUE){
  pdf("ICB.Riaz2017_Nivolumab_Melanoma_Naive_ROC.pdf",height=6,width=6)
  auc.out <- c()
  
  #先画第一条线，此处是miRNA1
  x <- plot.roc(df[,1],df[,2],ylim=c(0,1),xlim=c(1,0),
                smooth=F, #绘制平滑曲线
                ci=TRUE, 
                main="",
                #print.thres="best", #把阈值写在图上，其sensitivity+ specificity之和最大
                col=mycol[2],#线的颜色
                lwd=2, #线的粗细
                legacy.axes=T)#采用大多数paper的画法，横坐标是“1-specificity”，从0到1
  
  ci.lower <- round(as.numeric(x$ci[1]),3) #置信区间下限
  ci.upper <- round(as.numeric(x$ci[3]),3) #置信区间上限
  
  auc.ci <- c(colnames(df)[2],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
  auc.out <- rbind(auc.out,auc.ci)
  
  
  #再用循环画第二条和后面更多条曲线
  for (i in 3:ncol(df)){
    x <- plot.roc(df[,1],df[,i],
                  add=T, #向前面画的图里添加
                  smooth=F,
                  ci=TRUE,
                  col=mycol[i],
                  lwd=2,
                  legacy.axes=T)
    
    ci.lower <- round(as.numeric(x$ci[1]),3)
    ci.upper <- round(as.numeric(x$ci[3]),3)
    
    auc.ci <- c(colnames(df)[i],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
    auc.out <- rbind(auc.out,auc.ci)
  }
  
  
  # #对比多条曲线
  # #在参数`method=`后面，有三种方法可选“delong”, “bootstrap”或“venkatraman”，计算p值
  # p.out <- c()
  # for (i in 2:(ncol(df)-1)){
  #   for (j in (i+1):ncol(df)){
  #     p <- roc.test(df[,1],df[,i],df[,j], method="bootstrap")
  #     p.tmp <- c(colnames(df)[i],colnames(df)[j],p$p.value)
  #     p.out <- rbind(p.out,p.tmp)
  #   }
  # }
  # p.out <- as.data.frame(p.out)
  # colnames(p.out) <- c("ROC1","ROC2","p.value")
  auc.out <- as.data.frame(auc.out)
  colnames(auc.out) <- c("Name","AUC","AUC CI")
  #输出p value到文件
  
  #绘制图例
  legend.name <- paste(colnames(df)[2:length(df)],"AUC",auc.out$AUC,sep=" ")
  legend("bottomright", 
         legend=legend.name,
         col = mycol[2:length(df)],
         lwd = 2,
         bty="n")
  dev.off()
}


#生存分析
setwd(paste0("Y:\\lcd代码\\TCGA随意做\\TCGA任意做\\进阶版\\进阶/",gene,"/5免疫浸润分析/6免疫治疗"))
test<-rt
test<-na.omit(test)
res.cut <- surv_cutpoint(test, time = "PFS", event = "PFS.Event",
                         variables = gene)
summary(res.cut)
#Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)
#Fit survival curves and visualize
fit <- survfit(Surv(PFS,PFS.Event ) ~get(gene), data = res.cat)


pdf(file = "ICB.Riaz2017_Nivolumab_Melanoma_Naive_PFS生存分析.pdf",width = 5.6,height = 5.34,onefile = F)

ggsurvplot(fit,
           data=res.cat,
           risk.table = T,
           risk.table.col = 'strata',
           risk.table.y.text = F,
           risk.table.title = 'Number at risk',
           pval = TRUE,
           #conf.int = TRUE,
           xlab = 'Time (years)',
           ggtheme = theme_light(),
           surv.median.line = 'hv',
           title=paste( gene,'Survival',sep = ' ' ))

dev.off()

###################################################################################
rm(list = ls())
gene="Score"
setwd("Y:\\lcd代码\\TCGA随意做\\TCGA任意做\\进阶版\\进阶\\原始数据\\免疫治疗数据")
genelist_input <- fread(file="ICB.Riaz2017_Nivolumab_Melanoma_Prog.csv", header = T, sep=',', data.table = F)
genename <- as.character(genelist_input[,1]) #提取第一列基因名
#x - 基因组注释R包
#keys - 需要转换的基因列表
#keytype - 基因名类型
#columns - 希望返回的基因名类型数据，如返回NCBI基因ID使用 ENTREZID
gene_map <- biomaRt::select(org.Hs.eg.db, keys=genename, keytype="ENTREZID", columns=c("SYMBOL"))
colnames(gene_map)[1]<-"ID"
colnames(genelist_input)[1]<-"ID"
genelist_input$ID<-as.character(genelist_input$ID)
immudata<-inner_join(gene_map,genelist_input,by="ID")
immudata<-immudata[,-1]
immudata<-na.omit(immudata)
rownames(immudata)<-immudata$SYMBOL
immudata<-immudata[,-1]
immudata<-as.matrix(immudata)
immudata<-as.data.frame(t(immudata))

CCC=read.table("uniSigGeneExp.txt",header = T,sep = "\t")
CCC=CCC$id
CCC=intersect(colnames(immudata),CCC)
AAA=immudata[,CCC]
#PCA分析
pca=prcomp(AAA, scale=TRUE)
value=predict(pca)
score=value[,1]+value[,2]
score=as.data.frame(score)
score$ID=rownames(score)
immudata$ID=row.names(immudata)
immudata=inner_join(immudata,score,by="ID")
colnames(immudata)[which(colnames(immudata)=="score")]=gene
immudata<-immudata %>%
  dplyr::select(ID,gene,everything())


clindata<-read.csv("ICB.Riaz2017_Nivolumab_Melanoma_Prog_clin.csv",header = T,sep = ",")
colnames(clindata)[1]<-"ID"
AAA<-intersect(clindata$ID,immudata$ID)
immudata<-dplyr::filter(immudata,ID %in% AAA)
clindata<-dplyr::filter(clindata,ID %in% AAA)
rt<-inner_join(clindata,immudata,by="ID")

rt$OS<-rt$OS/365
rt$OS.Event<-as.numeric(rt$OS.Event)

#生存分析
setwd(paste0("Y:\\lcd代码\\TCGA随意做\\TCGA任意做\\进阶版\\进阶/",gene,"/5免疫浸润分析/6免疫治疗"))
test<-rt
test<-na.omit(test)
res.cut <- surv_cutpoint(test, time = "OS", event = "OS.Event",
                         variables = gene)
summary(res.cut)
#Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)
#Fit survival curves and visualize
fit <- survfit(Surv(OS,OS.Event ) ~get(gene), data = res.cat)


pdf(file = "ICB.Riaz2017_Nivolumab_Melanoma_Prog_OS生存分析.pdf",width = 5.6,height = 5.34,onefile = F)

ggsurvplot(fit,
           data=res.cat,
           risk.table = T,
           risk.table.col = 'strata',
           risk.table.y.text = F,
           risk.table.title = 'Number at risk',
           pval = TRUE,
           #conf.int = TRUE,
           xlab = 'Time (years)',
           ggtheme = theme_light(),
           surv.median.line = 'hv',
           title=paste( gene,'Survival',sep = ' ' ))

dev.off()

#百分比图
res.cat$ID<-test$ID
res.cat$Response<-test$Response
res.cat$Response[res.cat$Response==1]="Yes"
res.cat$Response[res.cat$Response==0]="No"
res.cat$Group<-res.cat[,gene]
library(dplyr)
library(plyr)
library(ggplot2)
library(ggpubr)
dat = dplyr::count(res.cat,Group,Response)
dat = dat %>% group_by(Group) %>% 
  dplyr::summarise(Response = Response,n = n/sum(n))
dat$Response = factor(dat$Response,levels = c("No","Yes"))


#计算百分率
dat=ddply(dat, .(Group), transform, percent = n/sum(n) * 100)
#百分比位置
dat=ddply(dat, .(Group), transform, pos = (cumsum(n) - 0.5 * n))
dat$label=paste0(sprintf("%.0f", dat$percent), "%")

bioCol=c("#f87669","#2fa1dd")
p=ggplot(dat, aes(x = factor(Group), y = percent, fill = Response)) +
  geom_bar(position = position_stack(), stat = "identity", width = .7) +
  scale_fill_manual(values=bioCol)+
  xlab("Group")+ ylab("Percent weight")+  guides(fill=guide_legend(title="Response"))+
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  #coord_flip()+
  theme_bw()
pdf(file="ICB.Riaz2017_Nivolumab_Melanoma_Prog_barplot.pdf", width=4, height=5)
print(p)
dev.off()

test<-rt
test<-na.omit(test)
test$Response[test$Response==1]="Yes"
test$Response[test$Response==0]="No"
test$Response<-factor(test$Response,levels = c('No','Yes'))
my_comparisons <- list( c("No", "Yes") )
#耐受非耐受差异表达
p <- ggboxplot(test, x = "Response", y = gene,
               fill = "Response",legend=F,palette =c("#E7B800", "#00AFBB"),bxp.errorbar=T)+
  theme(legend.position='none')+
  ylab(label = gene)+
  stat_compare_means(comparisons = my_comparisons,method = 't.test',aes(label = ..p.signif..))
ggsave(p,filename = paste("ICB.Riaz2017_Nivolumab_Melanoma_Prog",gene,'差异表达.pdf',sep=' '),width = 5.6,height = 4.22)



#预测反应性ROC曲线
library("pROC")
df<-dplyr::select(rt,Response,gene,PDCD1,CTLA4,CD274)
#定义足够多的颜色，后面画线时从这里选颜色
mycol <- c("slateblue","seagreen3","dodgerblue","firebrick1","lightgoldenrod","magenta","orange2")

if(TRUE){
  pdf("ICB.Riaz2017_Nivolumab_Melanoma_Prog_ROC.pdf",height=6,width=6)
  auc.out <- c()
  
  #先画第一条线，此处是miRNA1
  x <- plot.roc(df[,1],df[,2],ylim=c(0,1),xlim=c(1,0),
                smooth=F, #绘制平滑曲线
                ci=TRUE, 
                main="",
                #print.thres="best", #把阈值写在图上，其sensitivity+ specificity之和最大
                col=mycol[2],#线的颜色
                lwd=2, #线的粗细
                legacy.axes=T)#采用大多数paper的画法，横坐标是“1-specificity”，从0到1
  
  ci.lower <- round(as.numeric(x$ci[1]),3) #置信区间下限
  ci.upper <- round(as.numeric(x$ci[3]),3) #置信区间上限
  
  auc.ci <- c(colnames(df)[2],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
  auc.out <- rbind(auc.out,auc.ci)
  
  
  #再用循环画第二条和后面更多条曲线
  for (i in 3:ncol(df)){
    x <- plot.roc(df[,1],df[,i],
                  add=T, #向前面画的图里添加
                  smooth=F,
                  ci=TRUE,
                  col=mycol[i],
                  lwd=2,
                  legacy.axes=T)
    
    ci.lower <- round(as.numeric(x$ci[1]),3)
    ci.upper <- round(as.numeric(x$ci[3]),3)
    
    auc.ci <- c(colnames(df)[i],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
    auc.out <- rbind(auc.out,auc.ci)
  }
  
  
  # #对比多条曲线
  # #在参数`method=`后面，有三种方法可选“delong”, “bootstrap”或“venkatraman”，计算p值
  # p.out <- c()
  # for (i in 2:(ncol(df)-1)){
  #   for (j in (i+1):ncol(df)){
  #     p <- roc.test(df[,1],df[,i],df[,j], method="bootstrap")
  #     p.tmp <- c(colnames(df)[i],colnames(df)[j],p$p.value)
  #     p.out <- rbind(p.out,p.tmp)
  #   }
  # }
  # p.out <- as.data.frame(p.out)
  # colnames(p.out) <- c("ROC1","ROC2","p.value")
  auc.out <- as.data.frame(auc.out)
  colnames(auc.out) <- c("Name","AUC","AUC CI")
  #输出p value到文件
  
  #绘制图例
  legend.name <- paste(colnames(df)[2:length(df)],"AUC",auc.out$AUC,sep=" ")
  legend("bottomright", 
         legend=legend.name,
         col = mycol[2:length(df)],
         lwd = 2,
         bty="n")
  dev.off()
}

#生存分析
setwd(paste0("Y:\\lcd代码\\TCGA随意做\\TCGA任意做\\进阶版\\进阶/",gene,"/5免疫浸润分析/6免疫治疗"))
rt$PFS<-rt$PFS/365
rt$PFS.Event<-as.numeric(rt$PFS.Event)
test<-rt
test<-na.omit(test)
res.cut <- surv_cutpoint(test, time = "PFS", event = "PFS.Event",
                         variables = gene)
summary(res.cut)
#Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)
#Fit survival curves and visualize
fit <- survfit(Surv(PFS,PFS.Event ) ~get(gene), data = res.cat)


pdf(file = "ICB.Riaz2017_Nivolumab_Melanoma_Prog_PFS生存分析.pdf",width = 5.6,height = 5.34,onefile = F)

ggsurvplot(fit,
           data=res.cat,
           risk.table = T,
           risk.table.col = 'strata',
           risk.table.y.text = F,
           risk.table.title = 'Number at risk',
           pval = TRUE,
           #conf.int = TRUE,
           xlab = 'Time (years)',
           ggtheme = theme_light(),
           surv.median.line = 'hv',
           title=paste( gene,'Survival',sep = ' ' ))

dev.off()

###################################################################################
rm(list = ls())
gene="Score"
setwd("Y:\\lcd代码\\TCGA随意做\\TCGA任意做\\进阶版\\进阶\\原始数据\\免疫治疗数据")
genelist_input <- fread(file="GSE135222.csv", header = T, sep=',', data.table = F)
genelist_input<-separate(genelist_input,gene_id, c("ID", NA),sep = "\\.")
genename <- as.character(genelist_input[,1]) #提取第一列基因名
#x - 基因组注释R包
#keys - 需要转换的基因列表
#keytype - 基因名类型
#columns - 希望返回的基因名类型数据，如返回NCBI基因ID使用 ENTREZID
gene_map <- biomaRt::select(org.Hs.eg.db, keys=genename, keytype="ENSEMBL", columns=c("SYMBOL"))
colnames(gene_map)[1]<-"ID"
colnames(genelist_input)[1]<-"ID"
genelist_input$ID<-as.character(genelist_input$ID)
immudata<-inner_join(gene_map,genelist_input,by="ID")
immudata<-immudata[,-1]
immudata<-na.omit(immudata)

#构建一个布尔向量，索引
index<-duplicated(immudata$SYMBOL)
index
#筛选数据
immudata<-immudata[!index,]  #选中了非重复的数据
immudata<-na.omit(immudata)
rownames(immudata)<-immudata$SYMBOL
immudata<-immudata[,-1]
immudata<-as.matrix(immudata)
immudata<-as.data.frame(t(immudata))

CCC=read.table("uniSigGeneExp.txt",header = T,sep = "\t")
CCC=CCC$id
CCC=intersect(colnames(immudata),CCC)

AAA=immudata[,CCC]
#PCA分析
pca=prcomp(AAA, scale=TRUE)
value=predict(pca)
score=value[,1]+value[,2]
score=as.data.frame(score)
score$ID=rownames(score)
immudata$ID=row.names(immudata)
immudata=inner_join(immudata,score,by="ID")
colnames(immudata)[which(colnames(immudata)=="score")]=gene
immudata<-immudata %>%
  dplyr::select(ID,gene,everything())


clindata<-read.csv("GSE135222_clin.csv",header = T,sep = ",")
colnames(clindata)[1]<-"ID"
AAA<-intersect(clindata$ID,immudata$ID)
immudata<-dplyr::filter(immudata,ID %in% AAA)
clindata<-dplyr::filter(clindata,ID %in% AAA)
rt<-inner_join(clindata,immudata,by="ID")

rt$PFS<-rt$PFS/365
rt$PFS.Event<-as.numeric(rt$PFS.Event)

#生存分析
setwd(paste0("Y:\\lcd代码\\TCGA随意做\\TCGA任意做\\进阶版\\进阶/",gene,"/5免疫浸润分析/6免疫治疗"))
test<-rt
test<-na.omit(test)
res.cut <- surv_cutpoint(test, time = "PFS", event = "PFS.Event",
                         variables = gene)
summary(res.cut)
#Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)
#Fit survival curves and visualize
fit <- survfit(Surv(PFS,PFS.Event ) ~get(gene), data = res.cat)


pdf(file = "GSE135222_PFS生存分析.pdf",width = 5.6,height = 5.34,onefile = F)

ggsurvplot(fit,
           data=res.cat,
           risk.table = T,
           risk.table.col = 'strata',
           risk.table.y.text = F,
           risk.table.title = 'Number at risk',
           pval = TRUE,
           #conf.int = TRUE,
           xlab = 'Time (years)',
           ggtheme = theme_light(),
           surv.median.line = 'hv',
           title=paste( gene,'Survival',sep = ' ' ))

dev.off()

#百分比图
res.cat$ID<-test$ID
res.cat$Progression<-test$PFS.Event
res.cat$Progression[res.cat$Progression==1]="Yes"
res.cat$Progression[res.cat$Progression==0]="No"
res.cat$Group<-res.cat[,gene]
library(dplyr)
library(plyr)
library(ggplot2)
library(ggpubr)
dat = dplyr::count(res.cat,Group,Progression)
dat = dat %>% group_by(Group) %>% 
  dplyr::summarise(Progression = Progression,n = n/sum(n))
dat$Progression = factor(dat$Progression,levels = c("No","Yes"))


#计算百分率
dat=ddply(dat, .(Group), transform, percent = n/sum(n) * 100)
#百分比位置
dat=ddply(dat, .(Group), transform, pos = (cumsum(n) - 0.5 * n))
dat$label=paste0(sprintf("%.0f", dat$percent), "%")

bioCol=c("#f87669","#2fa1dd")
p=ggplot(dat, aes(x = factor(Group), y = percent, fill = Progression)) +
  geom_bar(position = position_stack(), stat = "identity", width = .7) +
  scale_fill_manual(values=bioCol)+
  xlab("Group")+ ylab("Percent weight")+  guides(fill=guide_legend(title="Progression"))+
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  #coord_flip()+
  theme_bw()
pdf(file="GSE135222_barplot.pdf", width=4, height=5)
print(p)
dev.off()

test<-rt
test<-na.omit(test)
test$PFS.Event[test$PFS.Event==1]="Yes"
test$PFS.Event[test$PFS.Event==0]="No"
test$PFS.Event<-factor(test$PFS.Event,levels = c('No','Yes'))
my_comparisons <- list( c("No", "Yes") )
#耐受非耐受差异表达
p <- ggboxplot(test, x = "PFS.Event", y = gene,
               fill = "PFS.Event",legend=F,palette =c("#E7B800", "#00AFBB"),bxp.errorbar=T)+
  theme(legend.position='none')+
  ylab(label = gene)+
  xlab(label = "Progression")+
  stat_compare_means(comparisons = my_comparisons,method = 't.test',aes(label = ..p.signif..))
ggsave(p,filename = paste("GSE135222",gene,'差异表达.pdf',sep=' '),width = 5.6,height = 4.22)


#预测反应性ROC曲线
library("pROC")
df<-dplyr::select(rt,PFS.Event,gene,PDCD1,CTLA4,CD274)
#定义足够多的颜色，后面画线时从这里选颜色
mycol <- c("slateblue","seagreen3","dodgerblue","firebrick1","lightgoldenrod","magenta","orange2")

if(TRUE){
  pdf("GSE135222_ROC.pdf",height=6,width=6)
  auc.out <- c()
  
  #先画第一条线，此处是miRNA1
  x <- plot.roc(df[,1],df[,2],ylim=c(0,1),xlim=c(1,0),
                smooth=F, #绘制平滑曲线
                ci=TRUE, 
                main="",
                #print.thres="best", #把阈值写在图上，其sensitivity+ specificity之和最大
                col=mycol[2],#线的颜色
                lwd=2, #线的粗细
                legacy.axes=T)#采用大多数paper的画法，横坐标是“1-specificity”，从0到1
  
  ci.lower <- round(as.numeric(x$ci[1]),3) #置信区间下限
  ci.upper <- round(as.numeric(x$ci[3]),3) #置信区间上限
  
  auc.ci <- c(colnames(df)[2],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
  auc.out <- rbind(auc.out,auc.ci)
  
  
  #再用循环画第二条和后面更多条曲线
  for (i in 3:ncol(df)){
    x <- plot.roc(df[,1],df[,i],
                  add=T, #向前面画的图里添加
                  smooth=F,
                  ci=TRUE,
                  col=mycol[i],
                  lwd=2,
                  legacy.axes=T)
    
    ci.lower <- round(as.numeric(x$ci[1]),3)
    ci.upper <- round(as.numeric(x$ci[3]),3)
    
    auc.ci <- c(colnames(df)[i],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
    auc.out <- rbind(auc.out,auc.ci)
  }
  
  
  #对比多条曲线
  #在参数`method=`后面，有三种方法可选“delong”, “bootstrap”或“venkatraman”，计算p值
  # p.out <- c()
  # for (i in 2:(ncol(df)-1)){
  #   for (j in (i+1):ncol(df)){
  #     p <- roc.test(df[,1],df[,i],df[,j], method="bootstrap")
  #     p.tmp <- c(colnames(df)[i],colnames(df)[j],p$p.value)
  #     p.out <- rbind(p.out,p.tmp)
  #   }
  # }
  # p.out <- as.data.frame(p.out)
  # colnames(p.out) <- c("ROC1","ROC2","p.value")
  auc.out <- as.data.frame(auc.out)
  colnames(auc.out) <- c("Name","AUC","AUC CI")
  #输出p value到文件
  
  #绘制图例
  legend.name <- paste(colnames(df)[2:length(df)],"AUC",auc.out$AUC,sep=" ")
  legend("bottomright", 
         legend=legend.name,
         col = mycol[2:length(df)],
         lwd = 2,
         bty="n")
  dev.off()
}

###################################################################################################
###################################################################################
rm(list = ls())
gene="Score"
setwd("Y:\\lcd代码\\TCGA随意做\\TCGA任意做\\进阶版\\进阶\\原始数据\\免疫治疗数据")
genelist_input <- fread(file="GSE176307.csv", header = T, sep=',', data.table = F)
rownames(genelist_input)<-genelist_input$ID_REF
genelist_input<-genelist_input[,-1]
immudata<-na.omit(genelist_input)
immudata<-log2(immudata+1)
immudata<-as.matrix(immudata)
immudata<-as.data.frame(t(immudata))

CCC=read.table("uniSigGeneExp.txt",header = T,sep = "\t")
CCC=CCC$id
CCC=intersect(colnames(immudata),CCC)

AAA=immudata[,CCC]
#PCA分析
pca=prcomp(AAA, scale=TRUE)
value=predict(pca)
score=value[,1]+value[,2]
score=as.data.frame(score)
score$ID=rownames(score)
immudata$ID=row.names(immudata)
immudata=inner_join(immudata,score,by="ID")
colnames(immudata)[which(colnames(immudata)=="score")]=gene
immudata<-immudata %>%
  dplyr::select(ID,gene,everything())


clindata<-read.csv("GSE176307_clin.csv",header = T,sep = ",")
colnames(clindata)[1]<-"ID"
AAA<-intersect(clindata$ID,immudata$ID)
immudata<-dplyr::filter(immudata,ID %in% AAA)
clindata<-dplyr::filter(clindata,ID %in% AAA)
rt<-inner_join(clindata,immudata,by="ID")

rt$OS<-rt$OS/365
rt$PFS<-rt$PFS/365
rt$status<-as.numeric(rt$status)
rt$Progression<-as.numeric(rt$Progression)
#生存分析
setwd(paste0("Y:\\lcd代码\\TCGA随意做\\TCGA任意做\\进阶版\\进阶/",gene,"/5免疫浸润分析/6免疫治疗"))
test<-rt
test<-dplyr::select(test,OS,status,gene,Response)
res.cut <- surv_cutpoint(test, time = "OS", event = "status",
                         variables = gene)
summary(res.cut)
#Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)
#Fit survival curves and visualize
fit <- survfit(Surv(OS,status ) ~get(gene), data = res.cat)


pdf(file = "GSE176307_OS生存分析.pdf",width = 5.6,height = 5.34,onefile = F)

ggsurvplot(fit,
           data=res.cat,
           risk.table = T,
           risk.table.col = 'strata',
           risk.table.y.text = F,
           risk.table.title = 'Number at risk',
           pval = TRUE,
           #conf.int = TRUE,
           xlab = 'Time (years)',
           ggtheme = theme_light(),
           surv.median.line = 'hv',
           title=paste( gene,'Survival',sep = ' ' ))

dev.off()

#百分比图
res.cat$ID<-test$ID
res.cat$Response<-test$Response
res.cat$Group<-res.cat[,gene]
res.cat<-na.omit(res.cat)
library(dplyr)
library(plyr)
library(ggplot2)
library(ggpubr)
dat = dplyr::count(res.cat,Group,Response)
dat = dat %>% group_by(Group) %>% 
  dplyr::summarise(Response = Response,n = n/sum(n))
dat$Response = factor(dat$Response,levels = c("CR/PR","PD/SD"))


#计算百分率
dat=ddply(dat, .(Group), transform, percent = n/sum(n) * 100)
#百分比位置
dat=ddply(dat, .(Group), transform, pos = (cumsum(n) - 0.5 * n))
dat$label=paste0(sprintf("%.0f", dat$percent), "%")

bioCol=c("#f87669","#2fa1dd")
p=ggplot(dat, aes(x = factor(Group), y = percent, fill = Response)) +
  geom_bar(position = position_stack(), stat = "identity", width = .7) +
  scale_fill_manual(values=bioCol)+
  xlab("Group")+ ylab("Percent weight")+  guides(fill=guide_legend(title="Response"))+
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  #coord_flip()+
  theme_bw()
pdf(file="GSE176307_barplot.pdf", width=4, height=5)
print(p)
dev.off()

test<-dplyr::select(rt,ID,Response,gene)
test<-na.omit(test)
test$Response<-factor(test$Response,levels = c('CR/PR','PD/SD'))
my_comparisons <- list( c("CR/PR", "PD/SD") )
#耐受非耐受差异表达
p <- ggboxplot(test, x = "Response", y = gene,
               fill = "Response",legend=F,palette =c("#E7B800", "#00AFBB"),bxp.errorbar=T)+
  theme(legend.position='none')+
  ylab(label = gene)+
  stat_compare_means(comparisons = my_comparisons,method = 't.test',aes(label = ..p.signif..))
ggsave(p,filename = paste("GSE176307",gene,'差异表达.pdf',sep=' '),width = 5.6,height = 4.22)



test<-dplyr::select(rt,PFS,Progression,gene)
res.cut <- surv_cutpoint(test, time = "PFS", event = "Progression",
                         variables = gene)
summary(res.cut)
#Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)
#Fit survival curves and visualize
fit <- survfit(Surv(PFS,Progression ) ~get(gene), data = res.cat)


pdf(file = "GSE176307_PFS生存分析.pdf",width = 5.6,height = 5.34,onefile = F)

ggsurvplot(fit,
           data=res.cat,
           risk.table = T,
           risk.table.col = 'strata',
           risk.table.y.text = F,
           risk.table.title = 'Number at risk',
           pval = TRUE,
           #conf.int = TRUE,
           xlab = 'Time (years)',
           ggtheme = theme_light(),
           surv.median.line = 'hv',
           title=paste( gene,'Survival',sep = ' ' ))

dev.off()

#预测反应性ROC曲线
library("pROC")
df<-dplyr::select(rt,Response,gene,PDCD1,CTLA4,CD274)
df<-na.omit(df)
df$Response[df$Response=="PD/SD"]=1
df$Response[df$Response=="CR/PR"]=0
#定义足够多的颜色，后面画线时从这里选颜色
mycol <- c("slateblue","seagreen3","dodgerblue","firebrick1","lightgoldenrod","magenta","orange2")

if(TRUE){
  pdf("GSE176307_ROC.pdf",height=6,width=6)
  auc.out <- c()
  
  #先画第一条线，此处是miRNA1
  x <- plot.roc(df[,1],df[,2],ylim=c(0,1),xlim=c(1,0),
                smooth=F, #绘制平滑曲线
                ci=TRUE, 
                main="",
                #print.thres="best", #把阈值写在图上，其sensitivity+ specificity之和最大
                col=mycol[2],#线的颜色
                lwd=2, #线的粗细
                legacy.axes=T)#采用大多数paper的画法，横坐标是“1-specificity”，从0到1
  
  ci.lower <- round(as.numeric(x$ci[1]),3) #置信区间下限
  ci.upper <- round(as.numeric(x$ci[3]),3) #置信区间上限
  
  auc.ci <- c(colnames(df)[2],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
  auc.out <- rbind(auc.out,auc.ci)
  
  
  #再用循环画第二条和后面更多条曲线
  for (i in 3:ncol(df)){
    x <- plot.roc(df[,1],df[,i],
                  add=T, #向前面画的图里添加
                  smooth=F,
                  ci=TRUE,
                  col=mycol[i],
                  lwd=2,
                  legacy.axes=T)
    
    ci.lower <- round(as.numeric(x$ci[1]),3)
    ci.upper <- round(as.numeric(x$ci[3]),3)
    
    auc.ci <- c(colnames(df)[i],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
    auc.out <- rbind(auc.out,auc.ci)
  }
  
  
  # #对比多条曲线
  # #在参数`method=`后面，有三种方法可选“delong”, “bootstrap”或“venkatraman”，计算p值
  # p.out <- c()
  # for (i in 2:(ncol(df)-1)){
  #   for (j in (i+1):ncol(df)){
  #     p <- roc.test(df[,1],df[,i],df[,j], method="bootstrap")
  #     p.tmp <- c(colnames(df)[i],colnames(df)[j],p$p.value)
  #     p.out <- rbind(p.out,p.tmp)
  #   }
  # }
  # p.out <- as.data.frame(p.out)
  # colnames(p.out) <- c("ROC1","ROC2","p.value")
  auc.out <- as.data.frame(auc.out)
  colnames(auc.out) <- c("Name","AUC","AUC CI")
  #输出p value到文件
  
  #绘制图例
  legend.name <- paste(colnames(df)[2:length(df)],"AUC",auc.out$AUC,sep=" ")
  legend("bottomright", 
         legend=legend.name,
         col = mycol[2:length(df)],
         lwd = 2,
         bty="n")
  dev.off()
}

###################################################################################################
###################################################################################
rm(list = ls())
gene="Score"
setwd("Y:\\lcd代码\\TCGA随意做\\TCGA任意做\\进阶版\\进阶\\原始数据\\免疫治疗数据")
genelist_input <- fread(file="NCT02684006.csv", header = T, sep=',', data.table = F)
rownames(genelist_input)<-genelist_input$HUGO
genelist_input<-genelist_input[,-1]
immudata<-na.omit(genelist_input)
range(immudata)
#immudata<-log2(immudata+1)
immudata<-as.matrix(immudata)
immudata<-as.data.frame(t(immudata))

CCC=read.table("uniSigGeneExp.txt",header = T,sep = "\t")
CCC=CCC$id
CCC=intersect(colnames(immudata),CCC)

AAA=immudata[,CCC]
#PCA分析
pca=prcomp(AAA, scale=TRUE)
value=predict(pca)
score=value[,1]+value[,2]
score=as.data.frame(score)
score$ID=rownames(score)
immudata$ID=row.names(immudata)
immudata=inner_join(immudata,score,by="ID")
colnames(immudata)[which(colnames(immudata)=="score")]=gene
immudata<-immudata %>%
  dplyr::select(ID,gene,everything())


clindata<-read.csv("NCT02684006_clin.csv",header = T,sep = ",")
colnames(clindata)[1]<-"ID"
AAA<-intersect(clindata$ID,immudata$ID)
immudata<-dplyr::filter(immudata,ID %in% AAA)
clindata<-dplyr::filter(clindata,ID %in% AAA)
rt<-inner_join(clindata,immudata,by="ID")

rt$PFS.time<-rt$PFS.time/12
rt$PFS<-as.numeric(rt$PFS)

#生存分析
setwd(paste0("Y:\\lcd代码\\TCGA随意做\\TCGA任意做\\进阶版\\进阶/",gene,"/5免疫浸润分析/6免疫治疗"))
test<-rt
test<-dplyr::select(test,ID,PFS.time,PFS,gene,Treatment)
test1<-dplyr::filter(test,Treatment=="Avelumab+Axitinib")
test2<-dplyr::filter(test,Treatment=="Sunitinib")

res.cut <- surv_cutpoint(test1, time = "PFS.time", event = "PFS",
                         variables = gene)
summary(res.cut)
#Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)
#Fit survival curves and visualize
fit <- survfit(Surv(PFS.time,PFS ) ~get(gene), data = res.cat)


pdf(file = "NCT02684006_免疫治疗_PFS生存分析.pdf",width = 5.6,height = 5.34,onefile = F)

ggsurvplot(fit,
           data=res.cat,
           risk.table = T,
           risk.table.col = 'strata',
           risk.table.y.text = F,
           risk.table.title = 'Number at risk',
           pval = TRUE,
           #conf.int = TRUE,
           xlab = 'Time (years)',
           ggtheme = theme_light(),
           surv.median.line = 'hv',
           title=paste( gene,'Survival',sep = ' ' ))

dev.off()

#百分比图
res.cat$ID<-test1$ID
res.cat$Progression<-test1$PFS
res.cat$Group<-res.cat[,gene]
res.cat<-na.omit(res.cat)
res.cat$Progression[res.cat$PFS==1]="Yes"
res.cat$Progression[res.cat$PFS==0]="No"
library(dplyr)
library(plyr)
library(ggplot2)
library(ggpubr)
dat = dplyr::count(res.cat,Group,Progression)
dat = dat %>% group_by(Group) %>% 
  dplyr::summarise(Progression = Progression,n = n/sum(n))
dat$Progression = factor(dat$Progression,levels = c("Yes","No"))


#计算百分率
dat=ddply(dat, .(Group), transform, percent = n/sum(n) * 100)
#百分比位置
dat=ddply(dat, .(Group), transform, pos = (cumsum(n) - 0.5 * n))
dat$label=paste0(sprintf("%.0f", dat$percent), "%")

bioCol=c("#f87669","#2fa1dd")
p=ggplot(dat, aes(x = factor(Group), y = percent, fill = Progression)) +
  geom_bar(position = position_stack(), stat = "identity", width = .7) +
  scale_fill_manual(values=bioCol)+
  xlab("Group")+ ylab("Percent weight")+  guides(fill=guide_legend(title="Progression"))+
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  #coord_flip()+
  theme_bw()
pdf(file="NCT02684006_免疫治疗_barplot.pdf", width=4, height=5)
print(p)
dev.off()

DAT<-dplyr::select(test1,ID,PFS,gene)
DAT<-na.omit(DAT)
colnames(DAT)[2]<-"Progression"
DAT$Progression[DAT$Progression==1]="Yes"
DAT$Progression[DAT$Progression==0]="No"


DAT$Progression<-factor(DAT$Progression,levels = c('Yes','No'))
my_comparisons <- list( c("Yes", "No") )
#耐受非耐受差异表达
p <- ggboxplot(DAT, x = "Progression", y = gene,
               fill = "Progression",legend=F,palette =c("#E7B800", "#00AFBB"),bxp.errorbar=T)+
  theme(legend.position='none')+
  ylab(label = gene)+
  stat_compare_means(comparisons = my_comparisons,method = 't.test',aes(label = ..p.signif..))
ggsave(p,filename = paste("NCT02684006_免疫治疗",gene,'差异表达.pdf',sep=' '),width = 5.6,height = 4.22)


#预测反应性ROC曲线
library("pROC")
df<-rt %>%
  dplyr::filter(Treatment=="Avelumab+Axitinib") %>%
  dplyr::select(ID,PFS,gene,CTLA4,CD274)


df<-na.omit(df)
colnames(df)[2]<-"Progression"
df<-df[,-1]
#定义足够多的颜色，后面画线时从这里选颜色
mycol <- c("slateblue","seagreen3","dodgerblue","firebrick1","lightgoldenrod","magenta","orange2")

if(TRUE){
  pdf("NCT02684006_免疫治疗_ROC.pdf",height=6,width=6)
  auc.out <- c()
  
  #先画第一条线，此处是miRNA1
  x <- plot.roc(df[,1],df[,2],ylim=c(0,1),xlim=c(1,0),
                smooth=F, #绘制平滑曲线
                ci=TRUE, 
                main="",
                #print.thres="best", #把阈值写在图上，其sensitivity+ specificity之和最大
                col=mycol[2],#线的颜色
                lwd=2, #线的粗细
                legacy.axes=T)#采用大多数paper的画法，横坐标是“1-specificity”，从0到1
  
  ci.lower <- round(as.numeric(x$ci[1]),3) #置信区间下限
  ci.upper <- round(as.numeric(x$ci[3]),3) #置信区间上限
  
  auc.ci <- c(colnames(df)[2],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
  auc.out <- rbind(auc.out,auc.ci)
  
  
  #再用循环画第二条和后面更多条曲线
  for (i in 3:ncol(df)){
    x <- plot.roc(df[,1],df[,i],
                  add=T, #向前面画的图里添加
                  smooth=F,
                  ci=TRUE,
                  col=mycol[i],
                  lwd=2,
                  legacy.axes=T)
    
    ci.lower <- round(as.numeric(x$ci[1]),3)
    ci.upper <- round(as.numeric(x$ci[3]),3)
    
    auc.ci <- c(colnames(df)[i],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
    auc.out <- rbind(auc.out,auc.ci)
  }
  
  
  # #对比多条曲线
  # #在参数`method=`后面，有三种方法可选“delong”, “bootstrap”或“venkatraman”，计算p值
  # p.out <- c()
  # for (i in 2:(ncol(df)-1)){
  #   for (j in (i+1):ncol(df)){
  #     p <- roc.test(df[,1],df[,i],df[,j], method="bootstrap")
  #     p.tmp <- c(colnames(df)[i],colnames(df)[j],p$p.value)
  #     p.out <- rbind(p.out,p.tmp)
  #   }
  # }
  # p.out <- as.data.frame(p.out)
  # colnames(p.out) <- c("ROC1","ROC2","p.value")
  auc.out <- as.data.frame(auc.out)
  colnames(auc.out) <- c("Name","AUC","AUC CI")
  #输出p value到文件
  
  #绘制图例
  legend.name <- paste(colnames(df)[2:length(df)],"AUC",auc.out$AUC,sep=" ")
  legend("bottomright", 
         legend=legend.name,
         col = mycol[2:length(df)],
         lwd = 2,
         bty="n")
  dev.off()
}

res.cut <- surv_cutpoint(test2, time = "PFS.time", event = "PFS",
                         variables = gene)
summary(res.cut)
#Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)
#Fit survival curves and visualize
fit <- survfit(Surv(PFS.time,PFS ) ~get(gene), data = res.cat)


pdf(file = "NCT02684006_Sunitinib_PFS生存分析.pdf",width = 5.6,height = 5.34,onefile = F)

ggsurvplot(fit,
           data=res.cat,
           risk.table = T,
           risk.table.col = 'strata',
           risk.table.y.text = F,
           risk.table.title = 'Number at risk',
           pval = TRUE,
           #conf.int = TRUE,
           xlab = 'Time (years)',
           ggtheme = theme_light(),
           surv.median.line = 'hv',
           title=paste( gene,'Survival',sep = ' ' ))

dev.off()

#百分比图
res.cat$ID<-test2$ID
res.cat$Progression<-test2$PFS
res.cat$Group<-res.cat[,gene]
res.cat<-na.omit(res.cat)
res.cat$Progression[res.cat$PFS==1]="Yes"
res.cat$Progression[res.cat$PFS==0]="No"
library(dplyr)
library(plyr)
library(ggplot2)
library(ggpubr)
dat = dplyr::count(res.cat,Group,Progression)
dat = dat %>% group_by(Group) %>% 
  dplyr::summarise(Progression = Progression,n = n/sum(n))
dat$Progression = factor(dat$Progression,levels = c("Yes","No"))


#计算百分率
dat=ddply(dat, .(Group), transform, percent = n/sum(n) * 100)
#百分比位置
dat=ddply(dat, .(Group), transform, pos = (cumsum(n) - 0.5 * n))
dat$label=paste0(sprintf("%.0f", dat$percent), "%")

bioCol=c("#f87669","#2fa1dd")
p=ggplot(dat, aes(x = factor(Group), y = percent, fill = Progression)) +
  geom_bar(position = position_stack(), stat = "identity", width = .7) +
  scale_fill_manual(values=bioCol)+
  xlab("Group")+ ylab("Percent weight")+  guides(fill=guide_legend(title="Progression"))+
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  #coord_flip()+
  theme_bw()
pdf(file="NCT02684006_Sunitinib_barplot.pdf", width=4, height=5)
print(p)
dev.off()

DAT<-dplyr::select(test2,ID,PFS,gene)
DAT<-na.omit(DAT)
colnames(DAT)[2]<-"Progression"
DAT$Progression[DAT$Progression==1]="Yes"
DAT$Progression[DAT$Progression==0]="No"


DAT$Progression<-factor(DAT$Progression,levels = c('Yes','No'))
my_comparisons <- list( c("Yes", "No") )
#耐受非耐受差异表达
p <- ggboxplot(DAT, x = "Progression", y = gene,
               fill = "Progression",legend=F,palette =c("#E7B800", "#00AFBB"),bxp.errorbar=T)+
  theme(legend.position='none')+
  ylab(label = gene)+
  stat_compare_means(comparisons = my_comparisons,method = 't.test',aes(label = ..p.signif..))
ggsave(p,filename = paste("NCT02684006_Sunitinib",gene,'差异表达.pdf',sep=' '),width = 5.6,height = 4.22)


#预测反应性ROC曲线
library("pROC")
df<-rt %>%
  dplyr::filter(Treatment=="Sunitinib") %>%
  dplyr::select(ID,PFS,gene,CTLA4,CD274)


df<-na.omit(df)
colnames(df)[2]<-"Progression"
df<-df[,-1]
#定义足够多的颜色，后面画线时从这里选颜色
mycol <- c("slateblue","seagreen3","dodgerblue","firebrick1","lightgoldenrod","magenta","orange2")

if(TRUE){
  pdf("NCT02684006_Sunitinib_ROC.pdf",height=6,width=6)
  auc.out <- c()
  
  #先画第一条线，此处是miRNA1
  x <- plot.roc(df[,1],df[,2],ylim=c(0,1),xlim=c(1,0),
                smooth=F, #绘制平滑曲线
                ci=TRUE, 
                main="",
                #print.thres="best", #把阈值写在图上，其sensitivity+ specificity之和最大
                col=mycol[2],#线的颜色
                lwd=2, #线的粗细
                legacy.axes=T)#采用大多数paper的画法，横坐标是“1-specificity”，从0到1
  
  ci.lower <- round(as.numeric(x$ci[1]),3) #置信区间下限
  ci.upper <- round(as.numeric(x$ci[3]),3) #置信区间上限
  
  auc.ci <- c(colnames(df)[2],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
  auc.out <- rbind(auc.out,auc.ci)
  
  
  #再用循环画第二条和后面更多条曲线
  for (i in 3:ncol(df)){
    x <- plot.roc(df[,1],df[,i],
                  add=T, #向前面画的图里添加
                  smooth=F,
                  ci=TRUE,
                  col=mycol[i],
                  lwd=2,
                  legacy.axes=T)
    
    ci.lower <- round(as.numeric(x$ci[1]),3)
    ci.upper <- round(as.numeric(x$ci[3]),3)
    
    auc.ci <- c(colnames(df)[i],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
    auc.out <- rbind(auc.out,auc.ci)
  }
  
  
  # #对比多条曲线
  # #在参数`method=`后面，有三种方法可选“delong”, “bootstrap”或“venkatraman”，计算p值
  # p.out <- c()
  # for (i in 2:(ncol(df)-1)){
  #   for (j in (i+1):ncol(df)){
  #     p <- roc.test(df[,1],df[,i],df[,j], method="bootstrap")
  #     p.tmp <- c(colnames(df)[i],colnames(df)[j],p$p.value)
  #     p.out <- rbind(p.out,p.tmp)
  #   }
  # }
  # p.out <- as.data.frame(p.out)
  # colnames(p.out) <- c("ROC1","ROC2","p.value")
  auc.out <- as.data.frame(auc.out)
  colnames(auc.out) <- c("Name","AUC","AUC CI")
  #输出p value到文件
  
  #绘制图例
  legend.name <- paste(colnames(df)[2:length(df)],"AUC",auc.out$AUC,sep=" ")
  legend("bottomright", 
         legend=legend.name,
         col = mycol[2:length(df)],
         lwd = 2,
         bty="n")
  dev.off()
}

###################################################################################################
###################################################################################
rm(list = ls())
gene="Score"

setwd("Y:\\lcd代码\\TCGA随意做\\TCGA任意做\\进阶版\\进阶\\原始数据\\免疫治疗数据")
genelist_input <- fread(file="PMID32472114.csv", header = T, sep=',', data.table = F)
rownames(genelist_input)<-genelist_input$gene_name
genelist_input<-genelist_input[,-1]
immudata<-na.omit(genelist_input)
range(immudata)
#immudata<-log2(immudata+1)
immudata<-as.matrix(immudata)
immudata<-as.data.frame(t(immudata))

CCC=read.table("uniSigGeneExp.txt",header = T,sep = "\t")
CCC=CCC$id
CCC=intersect(colnames(immudata),CCC)

AAA=immudata[,CCC]
#PCA分析
pca=prcomp(AAA, scale=TRUE)
value=predict(pca)
score=value[,1]+value[,2]
score=as.data.frame(score)
score$ID=rownames(score)
immudata$ID=row.names(immudata)
immudata=inner_join(immudata,score,by="ID")
colnames(immudata)[which(colnames(immudata)=="score")]=gene
immudata<-immudata %>%
  dplyr::select(ID,gene,everything())


clindata<-read.csv("PMID32472114_clin.csv",header = T,sep = ",")
colnames(clindata)[1]<-"ID"
AAA<-intersect(clindata$ID,immudata$ID)
immudata<-dplyr::filter(immudata,ID %in% AAA)
clindata<-dplyr::filter(clindata,ID %in% AAA)
rt<-inner_join(clindata,immudata,by="ID")

rt$PFS.time<-rt$PFS.time/12
rt$PFS<-as.numeric(rt$PFS)
rt$OS.time<-rt$OS.time/12
rt$OS<-as.numeric(rt$OS)
#生存分析
setwd(paste0("Y:\\lcd代码\\TCGA随意做\\TCGA任意做\\进阶版\\进阶/",gene,"/5免疫浸润分析/6免疫治疗"))

test<-dplyr::select(rt,ID,OS.time,OS,gene,Treatment,ORR)
test1<-dplyr::filter(test,Treatment=="Nivolumab")
test2<-dplyr::filter(test,Treatment=="Everolimus")

res.cut <- surv_cutpoint(test1, time = "OS.time", event = "OS",
                         variables = gene)
summary(res.cut)
#Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)
#Fit survival curves and visualize
fit <- survfit(Surv(OS.time,OS ) ~get(gene), data = res.cat)


pdf(file = "PMID32472114_免疫治疗_OS生存分析.pdf",width = 5.6,height = 5.34,onefile = F)

ggsurvplot(fit,
           data=res.cat,
           risk.table = T,
           risk.table.col = 'strata',
           risk.table.y.text = F,
           risk.table.title = 'Number at risk',
           pval = TRUE,
           #conf.int = TRUE,
           xlab = 'Time (years)',
           ggtheme = theme_light(),
           surv.median.line = 'hv',
           title=paste( gene,'Survival',sep = ' ' ))

dev.off()

#百分比图
res.cat$ID<-test1$ID
res.cat$Response<-test1$ORR
res.cat$Group<-res.cat[,gene]
res.cat<-na.omit(res.cat)
library(dplyr)
library(plyr)
library(ggplot2)
library(ggpubr)
dat = dplyr::count(res.cat,Group,Response)
dat = dat %>% group_by(Group) %>% 
  dplyr::summarise(Response = Response,n = n/sum(n))
dat$Response = factor(dat$Response,levels = c("CR/PR","PD/SD"))


#计算百分率
dat=ddply(dat, .(Group), transform, percent = n/sum(n) * 100)
#百分比位置
dat=ddply(dat, .(Group), transform, pos = (cumsum(n) - 0.5 * n))
dat$label=paste0(sprintf("%.0f", dat$percent), "%")

bioCol=c("#f87669","#2fa1dd")
p=ggplot(dat, aes(x = factor(Group), y = percent, fill = Response)) +
  geom_bar(position = position_stack(), stat = "identity", width = .7) +
  scale_fill_manual(values=bioCol)+
  xlab("Group")+ ylab("Percent weight")+  guides(fill=guide_legend(title="Response"))+
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  #coord_flip()+
  theme_bw()
pdf(file="PMID32472114_免疫治疗_barplot.pdf", width=4, height=5)
print(p)
dev.off()

DAT<-dplyr::select(test1,ID,ORR,gene)
DAT<-na.omit(DAT)
colnames(DAT)[2]<-"Response"


DAT$Response<-factor(DAT$Response,levels = c('CR/PR','PD/SD'))
my_comparisons <- list( c("CR/PR", "PD/SD") )
#耐受非耐受差异表达
p <- ggboxplot(DAT, x = "Response", y = gene,
               fill = "Response",legend=F,palette =c("#E7B800", "#00AFBB"),bxp.errorbar=T)+
  theme(legend.position='none')+
  ylab(label = gene)+
  stat_compare_means(comparisons = my_comparisons,method = 't.test',aes(label = ..p.signif..))
ggsave(p,filename = paste("PMID32472114_免疫治疗",gene,'差异表达.pdf',sep=' '),width = 5.6,height = 4.22)


#预测反应性ROC曲线
library("pROC")
df<-rt %>%
  dplyr::filter(Treatment=="Nivolumab") %>%
  dplyr::select(ID,ORR,gene,PDCD1,CTLA4,CD274) %>%
  na.omit()

df$ORR[df$ORR=="PD/SD"]=1
df$ORR[df$ORR=="CR/PR"]=0
df<-df[,-1]
#定义足够多的颜色，后面画线时从这里选颜色
mycol <- c("slateblue","seagreen3","dodgerblue","firebrick1","lightgoldenrod","magenta","orange2")

if(TRUE){
  pdf("PMID32472114_免疫治疗_ROC.pdf",height=6,width=6)
  auc.out <- c()
  
  #先画第一条线，此处是miRNA1
  x <- plot.roc(df[,1],df[,2],ylim=c(0,1),xlim=c(1,0),
                smooth=F, #绘制平滑曲线
                ci=TRUE, 
                main="",
                #print.thres="best", #把阈值写在图上，其sensitivity+ specificity之和最大
                col=mycol[2],#线的颜色
                lwd=2, #线的粗细
                legacy.axes=T)#采用大多数paper的画法，横坐标是“1-specificity”，从0到1
  
  ci.lower <- round(as.numeric(x$ci[1]),3) #置信区间下限
  ci.upper <- round(as.numeric(x$ci[3]),3) #置信区间上限
  
  auc.ci <- c(colnames(df)[2],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
  auc.out <- rbind(auc.out,auc.ci)
  
  
  #再用循环画第二条和后面更多条曲线
  for (i in 3:ncol(df)){
    x <- plot.roc(df[,1],df[,i],
                  add=T, #向前面画的图里添加
                  smooth=F,
                  ci=TRUE,
                  col=mycol[i],
                  lwd=2,
                  legacy.axes=T)
    
    ci.lower <- round(as.numeric(x$ci[1]),3)
    ci.upper <- round(as.numeric(x$ci[3]),3)
    
    auc.ci <- c(colnames(df)[i],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
    auc.out <- rbind(auc.out,auc.ci)
  }
  
  
  # #对比多条曲线
  # #在参数`method=`后面，有三种方法可选“delong”, “bootstrap”或“venkatraman”，计算p值
  # p.out <- c()
  # for (i in 2:(ncol(df)-1)){
  #   for (j in (i+1):ncol(df)){
  #     p <- roc.test(df[,1],df[,i],df[,j], method="bootstrap")
  #     p.tmp <- c(colnames(df)[i],colnames(df)[j],p$p.value)
  #     p.out <- rbind(p.out,p.tmp)
  #   }
  # }
  # p.out <- as.data.frame(p.out)
  # colnames(p.out) <- c("ROC1","ROC2","p.value")
  auc.out <- as.data.frame(auc.out)
  colnames(auc.out) <- c("Name","AUC","AUC CI")
  #输出p value到文件
  
  #绘制图例
  legend.name <- paste(colnames(df)[2:length(df)],"AUC",auc.out$AUC,sep=" ")
  legend("bottomright", 
         legend=legend.name,
         col = mycol[2:length(df)],
         lwd = 2,
         bty="n")
  dev.off()
}


res.cut <- surv_cutpoint(test2, time = "OS.time", event = "OS",
                         variables = gene)
summary(res.cut)
#Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)
#Fit survival curves and visualize
fit <- survfit(Surv(OS.time,OS ) ~get(gene), data = res.cat)


pdf(file = "PMID32472114_Everolimus_OS生存分析.pdf",width = 5.6,height = 5.34,onefile = F)

ggsurvplot(fit,
           data=res.cat,
           risk.table = T,
           risk.table.col = 'strata',
           risk.table.y.text = F,
           risk.table.title = 'Number at risk',
           pval = TRUE,
           #conf.int = TRUE,
           xlab = 'Time (years)',
           ggtheme = theme_light(),
           surv.median.line = 'hv',
           title=paste( gene,'Survival',sep = ' ' ))

dev.off()

#百分比图
res.cat$ID<-test2$ID
res.cat$Progression<-test2$ORR
res.cat$Group<-res.cat[,gene]
res.cat<-na.omit(res.cat)

library(dplyr)
library(plyr)
library(ggplot2)
library(ggpubr)
dat = dplyr::count(res.cat,Group,Progression)
dat = dat %>% group_by(Group) %>% 
  dplyr::summarise(Progression = Progression,n = n/sum(n))
dat$Progression = factor(dat$Progression,levels = c("CR/PR","PD/SD"))


#计算百分率
dat=ddply(dat, .(Group), transform, percent = n/sum(n) * 100)
#百分比位置
dat=ddply(dat, .(Group), transform, pos = (cumsum(n) - 0.5 * n))
dat$label=paste0(sprintf("%.0f", dat$percent), "%")

bioCol=c("#f87669","#2fa1dd")
p=ggplot(dat, aes(x = factor(Group), y = percent, fill = Progression)) +
  geom_bar(position = position_stack(), stat = "identity", width = .7) +
  scale_fill_manual(values=bioCol)+
  xlab("Group")+ ylab("Percent weight")+  guides(fill=guide_legend(title="Progression"))+
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  #coord_flip()+
  theme_bw()
pdf(file="PMID32472114_Everolimus_barplot.pdf", width=4, height=5)
print(p)
dev.off()

DAT<-dplyr::select(test2,ID,ORR,gene)
DAT<-na.omit(DAT)
colnames(DAT)[2]<-"Progression"



DAT$Progression<-factor(DAT$Progression,levels = c('CR/PR','PD/SD'))
my_comparisons <- list( c("CR/PR", "PD/SD") )
#耐受非耐受差异表达
p <- ggboxplot(DAT, x = "Progression", y = gene,
               fill = "Progression",legend=F,palette =c("#E7B800", "#00AFBB"),bxp.errorbar=T)+
  theme(legend.position='none')+
  ylab(label = gene)+
  stat_compare_means(comparisons = my_comparisons,method = 't.test',aes(label = ..p.signif..))
ggsave(p,filename = paste("PMID32472114_Everolimus",gene,'差异表达.pdf',sep=' '),width = 5.6,height = 4.22)


#预测反应性ROC曲线
library("pROC")
df<-rt %>%
  dplyr::filter(Treatment=="Everolimus") %>%
  dplyr::select(ID,ORR,gene,CTLA4,CD274)


df<-na.omit(df)
colnames(df)[2]<-"Progression"
df$Progression[df$Progression=="PD/SD"]=1
df$Progression[df$Progression=="CR/PR"]=0
df<-df[,-1]
#定义足够多的颜色，后面画线时从这里选颜色
mycol <- c("slateblue","seagreen3","dodgerblue","firebrick1","lightgoldenrod","magenta","orange2")

if(TRUE){
  pdf("PMID32472114_Everolimus_ROC.pdf",height=6,width=6)
  auc.out <- c()
  
  #先画第一条线，此处是miRNA1
  x <- plot.roc(df[,1],df[,2],ylim=c(0,1),xlim=c(1,0),
                smooth=F, #绘制平滑曲线
                ci=TRUE, 
                main="",
                #print.thres="best", #把阈值写在图上，其sensitivity+ specificity之和最大
                col=mycol[2],#线的颜色
                lwd=2, #线的粗细
                legacy.axes=T)#采用大多数paper的画法，横坐标是“1-specificity”，从0到1
  
  ci.lower <- round(as.numeric(x$ci[1]),3) #置信区间下限
  ci.upper <- round(as.numeric(x$ci[3]),3) #置信区间上限
  
  auc.ci <- c(colnames(df)[2],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
  auc.out <- rbind(auc.out,auc.ci)
  
  
  #再用循环画第二条和后面更多条曲线
  for (i in 3:ncol(df)){
    x <- plot.roc(df[,1],df[,i],
                  add=T, #向前面画的图里添加
                  smooth=F,
                  ci=TRUE,
                  col=mycol[i],
                  lwd=2,
                  legacy.axes=T)
    
    ci.lower <- round(as.numeric(x$ci[1]),3)
    ci.upper <- round(as.numeric(x$ci[3]),3)
    
    auc.ci <- c(colnames(df)[i],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
    auc.out <- rbind(auc.out,auc.ci)
  }
  
  
  # #对比多条曲线
  # #在参数`method=`后面，有三种方法可选“delong”, “bootstrap”或“venkatraman”，计算p值
  # p.out <- c()
  # for (i in 2:(ncol(df)-1)){
  #   for (j in (i+1):ncol(df)){
  #     p <- roc.test(df[,1],df[,i],df[,j], method="bootstrap")
  #     p.tmp <- c(colnames(df)[i],colnames(df)[j],p$p.value)
  #     p.out <- rbind(p.out,p.tmp)
  #   }
  # }
  # p.out <- as.data.frame(p.out)
  # colnames(p.out) <- c("ROC1","ROC2","p.value")
  auc.out <- as.data.frame(auc.out)
  colnames(auc.out) <- c("Name","AUC","AUC CI")
  #输出p value到文件
  
  #绘制图例
  legend.name <- paste(colnames(df)[2:length(df)],"AUC",auc.out$AUC,sep=" ")
  legend("bottomright", 
         legend=legend.name,
         col = mycol[2:length(df)],
         lwd = 2,
         bty="n")
  dev.off()
}

test<-dplyr::select(rt,ID,PFS.time,PFS,gene,Treatment,ORR)
test1<-dplyr::filter(test,Treatment=="Nivolumab")
test2<-dplyr::filter(test,Treatment=="Everolimus")

res.cut <- surv_cutpoint(test1, time = "PFS.time", event = "PFS",
                         variables = gene)
summary(res.cut)
#Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)
#Fit survival curves and visualize
fit <- survfit(Surv(PFS.time,PFS ) ~get(gene), data = res.cat)


pdf(file = "PMID32472114_免疫治疗_PFS生存分析.pdf",width = 5.6,height = 5.34,onefile = F)

ggsurvplot(fit,
           data=res.cat,
           risk.table = T,
           risk.table.col = 'strata',
           risk.table.y.text = F,
           risk.table.title = 'Number at risk',
           pval = TRUE,
           #conf.int = TRUE,
           xlab = 'Time (years)',
           ggtheme = theme_light(),
           surv.median.line = 'hv',
           title=paste( gene,'Survival',sep = ' ' ))

dev.off()


res.cut <- surv_cutpoint(test2, time = "PFS.time", event = "PFS",
                         variables = gene)
summary(res.cut)
#Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)
#Fit survival curves and visualize
fit <- survfit(Surv(PFS.time,PFS ) ~get(gene), data = res.cat)


pdf(file = "PMID32472114_Everolimus_PFS生存分析.pdf",width = 5.6,height = 5.34,onefile = F)

ggsurvplot(fit,
           data=res.cat,
           risk.table = T,
           risk.table.col = 'strata',
           risk.table.y.text = F,
           risk.table.title = 'Number at risk',
           pval = TRUE,
           #conf.int = TRUE,
           xlab = 'Time (years)',
           ggtheme = theme_light(),
           surv.median.line = 'hv',
           title=paste( gene,'Survival',sep = ' ' ))

dev.off()



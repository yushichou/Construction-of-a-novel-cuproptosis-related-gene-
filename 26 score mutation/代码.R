rm(list = ls())
#加载包
library(tidyverse)
library(maftools)

drawdata<-read.table("score.group.txt",header = T,sep = "\t")
drawdata$id=substr(drawdata$id,1,12)
laml.maf <- data.table::fread(paste0("TCGA.","CESC",".varscan.gz"),data.table = F)

laml.maf$Tumor_Sample_Barcode<-substr(laml.maf$Tumor_Sample_Barcode,1,12)
laml.maf$Tumor_Sample_Barcode<-str_replace_all(laml.maf$Tumor_Sample_Barcode,"-",".")

newdata=drawdata
colnames(newdata)[colnames(newdata)=="id"]="ID"
value<-newdata$ID[which(newdata$ID %in% laml.maf$Tumor_Sample_Barcode)]
newdata<-newdata%>%
  dplyr::filter(ID %in%value )


Highexpr<-newdata$ID[which(newdata[,"group"]=="High")]
Lowexpr<-newdata$ID[which(newdata[,"group"]=="Low")]

High.laml.maf<-laml.maf%>%
  dplyr::filter(Tumor_Sample_Barcode %in% Highexpr )
Low.laml.maf<-laml.maf%>%
  dplyr::filter(Tumor_Sample_Barcode %in% Lowexpr )

High.laml = read.maf(maf = High.laml.maf)
Low.laml = read.maf(maf = Low.laml.maf)

vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)

pdf(file=paste0("CESC","高score组突变瀑布图.pdf"),width=8,height=6)
oncoplot(maf = High.laml, colors = vc_cols,top = 10)
dev.off()
pdf(file=paste0("CESC","低score组突变瀑布图.pdf"),width=8,height=6)
oncoplot(maf = Low.laml, colors = vc_cols, top = 10)
dev.off()
pt.vs.rt <- mafCompare(m1 = High.laml, m2 = Low.laml, m1Name = paste0("CESC","_High"), m2Name =  paste0("CESC","_Low"), minMut = 0)
write.csv(pt.vs.rt$results,file =paste0("CESC","高低score组差异突变基因对比.csv"))
#选取前12个基因展示
pt.vs.rt$results<-pt.vs.rt$results[1:12,]
pdf(file=paste0("CESC","高低score组差异突变基因对比.pdf"),width=8,height=8)
maftools::forestPlot(mafCompareRes = pt.vs.rt, pVal = 1, geneFontSize = 0.8)
dev.off()

  


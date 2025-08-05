aaa=read.table("TCGA-PRAD.MRNA.TXT",header = T,sep = "\t",row.names = 1)
aaa=as.data.frame(t(aaa))
rownames(aaa)
substr(rownames(aaa),14,15)
table(substr(rownames(aaa),14,15))
aaa$type=substr(rownames(aaa),14,15)
aaa$type=as.numeric(aaa$type)
aaa$tumor=NULL



aaa$type=as.character(aaa$type)
aaa$tumor=ifelse(aaa$type %in% c("11","06","1") ,"else","Tumor")





library(tidyverse)
aaa=dplyr::select(aaa,tumor,type,everything())
aaa=dplyr::filter(aaa,tumor == "Tumor")
aaa=aaa[,-1]
aaa=aaa[,-1]
aaa=as.data.frame(t(aaa))
write.table(aaa,file = "TCGA.PRAD.TUMOR.txt",sep = "\t",row.names = T)











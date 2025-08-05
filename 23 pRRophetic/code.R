
#------------------------------------------------------------
rm(list = ls())
library(limma)
library(ggpubr)
library(pRRophetic)
library(ggplot2)
set.seed(12345)
expFile="merge.txt" 
load("merge.RDATA")
rt = outTab
data=outTab
library(tidyverse)
colnames(data)=str_replace_all(colnames(data),"TCGA_","")
colnames(data)=str_replace_all(colnames(data),"GSE63514_","")

group=read.table("score.group.txt",header = T,sep = "\t")
group=dplyr::select(group,id,group)
risk=group
colnames(risk)[2]="risk"
rownames(risk)=risk$id
drugname<-c(  "A.443654", "A.770041", "ABT.263", "ABT.888", "AG.014699", "AICAR", "AKT.inhibitor.VIII", "AMG.706", "AP.24534",
              "AS601245", "ATRA", "AUY922", "Axitinib", "AZ628", "AZD.0530", "AZD.2281", "AZD6244", "AZD6482", "AZD7762","AZD8055",
              "BAY.61.3606", "Bexarotene", "BI.2536", "BIBW2992", "Bicalutamide", "BI.D1870", "BIRB.0796", "Bleomycin",
              "BMS.509744", "BMS.536924", "BMS.708163", "BMS.754807", "Bortezomib", "Bosutinib", "Bryostatin.1", "BX.795",
              "Camptothecin", "CCT007093", "CCT018159", "CEP.701", "CGP.082996", "CGP.60474", "CHIR.99021", "CI.1040", "Cisplatin",
              "CMK", "Cyclopamine", "Cytarabine", "Dasatinib", "DMOG", "Docetaxel", "Doxorubicin", "EHT.1864", "Elesclomol",
              "Embelin", "Epothilone.B", "Erlotinib", "Etoposide", "FH535", "FTI.277", "GDC.0449", "GDC0941", "Gefitinib",
              "Gemcitabine", "GNF.2", "GSK269962A", "GSK.650394", "GW.441756", "GW843682X", "Imatinib", "IPA.3",
              "JNJ.26854165", "JNK.9L", "JNK.Inhibitor.VIII", "JW.7.52.1", "KIN001.135", "KU.55933", "Lapatinib",
              "Lenalidomide", "LFM.A13", "Metformin", "Methotrexate", "MG.132", "Midostaurin", "Mitomycin.C", "MK.2206",
              "MS.275", "Nilotinib", "NSC.87877", "NU.7441", "Nutlin.3a", "NVP.BEZ235", "NVP.TAE684", "Obatoclax.Mesylate",
              "OSI.906", "PAC.1", "Paclitaxel", "Parthenolide", "Pazopanib", "PD.0325901", "PD.0332991", "PD.173074", "PF.02341066", 
              "PF.4708671", "PF.562271", "PHA.665752", "PLX4720", "Pyrimethamine", "QS11", "Rapamycin", "RDEA119", "RO.3306", 
              "Roscovitine", "Salubrinal", "SB.216763", "SB590885", "Shikonin", "SL.0101.1", "Sorafenib", "S.Trityl.L.cysteine", 
              "Sunitinib", "Temsirolimus", "Thapsigargin", "Tipifarnib", "TW.37", "Vinblastine", "Vinorelbine", "Vorinostat", 
              "VX.680", "VX.702", "WH.4.023", "WO2009093972", "WZ.1.84", "X17.AAG", "X681640", "XMD8.85", "Z.LLNle.CHO", "ZM.447439"          )


for (i in 1:138) {
  drug=drugname[[i]]
  senstivity=pRRopheticPredict(data, drug, selection=1)
  senstivity=senstivity[senstivity!="NaN"]
  sameSample=intersect(row.names(risk), names(senstivity))
  risk=risk[sameSample, "risk",drop=F]
  senstivity=senstivity[sameSample]
  rt=cbind(risk, senstivity)
  rt$risk=factor(rt$risk, levels=c("Low", "High"))
  type=levels(factor(rt[,"risk"]))
  comp=combn(type, 2)
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  boxplot=ggboxplot(rt, x="risk", y="senstivity", fill="risk",bxp.errorbar=T,
                    xlab="score",
                    ylab=paste0(drug, " senstivity (IC50)"),
                    legend.title="score",
                    palette =c("#00AFBB", "#E7B800")
  )+ 
    stat_compare_means(comparisons=my_comparisons)
  pdf(file=paste0(drug, ".pdf"), width=5, height=4.5)
  print(boxplot)
  dev.off()
}










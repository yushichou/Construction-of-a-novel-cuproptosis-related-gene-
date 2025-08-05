rm(list = ls())

#引用包
library(plyr)
library(ggplot2)
library(ggpubr)
scoreFile="score.group.txt"    #打分文件
cliFile="clinicaldata.txt"            #临床数据文件
trait="fustat"                    #临床性状


#读取输入文件
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(score), row.names(cli))
rt=cbind(score[sameSample,,drop=F], cli[sameSample,,drop=F])
rt$fustat[rt$fustat==1]="Dead"
rt$fustat[rt$fustat==0]="Alive"
#定义临床性状的颜色
bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(unique(rt[,trait]))]

#统计高低评分组病人数目
rt1=rt[,c(trait, "group")]
rt1<-na.omit(rt1)
colnames(rt1)=c("trait", "group")
df=as.data.frame(table(rt1))
#计算高低评分组的百分率
df=ddply(df, .(group), transform, percent = Freq/sum(Freq) * 100)
#百分比位置
df=ddply(df, .(group), transform, pos = (cumsum(Freq) - 0.5 * Freq))
df$label=paste0(sprintf("%.0f", df$percent), "%")
df$group=factor(df$group, levels=c("Low", "High"))

#绘制百分率图
p=ggplot(df, aes(x = factor(group), y = percent, fill = trait)) +
  geom_bar(position = position_stack(), stat = "identity", width = .7) +
  scale_fill_manual(values=bioCol)+
  xlab("score")+ ylab("Percent weight")+  guides(fill=guide_legend(title=trait))+
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  #coord_flip()+
  theme_bw()
pdf(file="fustat_barplot.pdf", width=4, height=5)
print(p)
dev.off()

#设置比较组
rt2=rt[,c(trait, "score")]
rt2<-na.omit(rt2)
colnames(rt2)=c("trait", "score")
type=levels(factor(rt2[,"trait"]))
rt2$trait=factor(rt2$trait,levels =type )
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
#绘制箱线图
boxplot=ggboxplot(rt2, x="trait", y="score", fill="trait",
                  xlab=trait,
                  ylab="score",
                  legend.title=trait,
                  palette=bioCol
)+ 
  stat_compare_means(comparisons=my_comparisons)
pdf(file="fustat_boxplot.pdf",width=5,height=4.5)
print(boxplot)
dev.off()

colnames(rt)


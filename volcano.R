######原始数据格式保存于 volcano格式.png

getwd()
setwd("/Users/liruiting")
#install.packages("ggplot2")
library(ggplot2)
data<-read.csv("HEME-DEG.csv",header=T,row.names = 1)
pdf("volcano_heme.pdf")
#定义“点”，修改’点‘的大小与透明度  参数:size=,alpha=
p <- ggplot()+geom_point(data,mapping=aes(x=log2(T_VS_C),y=-log10(pval),color=up.down), size=2, alpha = 0.6)
#修改’点’的颜色(参数：scale_colour_manual)
color=c("#00EEEE","grey","#FF6A6A")
p<-p+scale_colour_manual(values=color)
#增加图标题和修改x与y轴的名称
 p<-p+labs(title="Serum.vs.Con",x=expression(Log[2]*'Fold Change Serum.vs.Con'),y=expression(-Log[10]*'P value'))
#修改图的背景
 p<-p+theme(panel.background = element_blank(),axis.line = element_line(colour = "black"))
#增加辅助线 
 p<-p+geom_hline(aes(yintercept=-log10(0.05)),linetype=2)+    geom_vline(aes(xintercept=log2(1.5)),linetype=2)+
 geom_vline(aes(xintercept=-log2(1.5)),linetype=2)
 p
 dev.off()



all <- ggplot()+geom_point(data,mapping=aes(x=log2FoldChange,y=-log10(pval),color=marker), size=5, alpha = 0.6)+scale_colour_manual(values=color)+labs(title="Con.vs.Treat",x=expression(Log[2]*'Fold Change Con.vs.Treat'),y=expression(-Log[10]*'P value'))+theme(panel.background = element_blank(),axis.line = element_line(colour = "black"))+geom_hline(aes(yintercept=-log10(0.05)),linetype=2)+    geom_vline(aes(xintercept=log2(1.5)),linetype=2)+
     geom_vline(aes(xintercept=-log2(1.5)),linetype=2)
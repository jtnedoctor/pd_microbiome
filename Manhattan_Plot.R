setwd("~/project/16s/ndd/r_analysis")
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
#save(blautia_genes_cor_pd,blautia_genes_cor_pd_all,file = "blautia_genes_cor_pd.Rdata")
load("blautia_genes_cor_pd.Rdata")
#save(blautia_genes_cor_pd_all_biotype,file = "blautia_genes_cor_pd_all_biotype.Rdata")
load("blautia_genes_cor_pd_all_biotype.Rdata")
df <- blautia_genes_cor_pd_all_biotype
df$cor<- as.numeric(df$cor)
df$pvalue <- as.numeric(df$pvalue)
#挑选含量大于100的biotype分类
type100 <-as.data.frame(table(df$gene_biotype))$Var1[as.data.frame(table(df$gene_biotype))$Freq>100]
type100 <- as.character(type100)
df$Chromosome <- ifelse(df$gene_biotype %in% type100,df$gene_biotype,"others")
df$Chromosome <- factor(df$Chromosome)
pdf(file = 'manhattan_blautia_cor_allgenes_biotypes.pdf')
ggplot(df,aes(x=Chromosome,y=-log10(pvalue)))+  #trait1对应特征的p值，也可以是你感兴趣的参数
  geom_jitter(aes(color=Chromosome))+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle=60,hjust=1))+
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,20))+
  #scale_x_discrete(labels=paste0("Chr",c(1:18,"X")))+
  scale_x_discrete(labels=levels(df$Chromosome))+
  geom_text_repel(
    data=df[-log10(df$pvalue)>14,],  #差异最显著的结果基因
    aes(label=togene),  #不能用rownames,只能用数据框里面的某一列
    size=4.5,
    color="black",
    segment.color="black",show.legend=FALSE) +#添加关注的点的基因名
  labs(x=NULL,y="-log10(Pvalue)") +
  geom_hline(yintercept = 14,lty="dashed")
dev.off()
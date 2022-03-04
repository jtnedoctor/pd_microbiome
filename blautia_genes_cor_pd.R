setwd("~/project/16s/ndd/r_analysis")
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
#save(blautia_genes_cor_pd,blautia_genes_cor_pd_all,file = "blautia_genes_cor_pd.Rdata") #此处blautia_genes_cor_pd的r>0.3
load("blautia_genes_cor_pd.Rdata")
#save(res_rnaseq_metaallgenes,res_rnaseq_metaNAgenes,res_rnaseq_metasigdiffgenes.fixed,res_rnaseq_metasigdiffgenes.random,file = "res_rnaseq_metagenes.Rdata")
load("res_rnaseq_metagenes.Rdata")
allDEG <- res_rnaseq_metaallgenes[rownames(res_rnaseq_metaallgenes) %in% rownames(blautia_genes_cor_pd),]
allDEG$togene <- rownames(allDEG)
colnames(allDEG)[c(1,4)] <- c("logFC","adj.P.Val")
xMax=max(abs(allDEG$logFC))
yMax=max(-log10(allDEG$adj.P.Val))   
allDEG$change <- ifelse(allDEG$adj.P.Val < 0.05 & abs(allDEG$logFC) > 0,
                        ifelse(allDEG$logFC > 0,'UP','DOWN'),
                        'NOT')
table(allDEG$change)
pdf(file = 'volcano_DEGs_blautia_cor_genes_0.3_label0.5.pdf')
ggplot(data= allDEG, aes(x = logFC, y = -log10(adj.P.Val), color = change)) +
  geom_point(alpha=0.8, size = 1) +
  theme_bw(base_size = 15) +
  theme(plot.title=element_text(hjust=0.5),   #  标题居中
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) + # 网格线设置为空白
  geom_vline(xintercept= 0 ,linetype= 2 ) +
  scale_color_manual(name = "", 
                     values = c("red", "green", "black"),
                     limits = c("UP", "DOWN", "NOT")) +
  ylim(0,yMax) + 
  xlim(-xMax,xMax) +
  #geom_vline(xintercept=c(-0.3,0.3),lty=2,col="black",lwd=0.5) + #添加横线|logFoldChange|>0.25
  geom_hline(yintercept=-log10(0.05),lty=2,col="black",lwd=0.5) + #添加竖线padj<0.05
  geom_text_repel(
    data=allDEG[allDEG$adj.P.Val<0.05 & abs(allDEG$logFC)>0.5,],
    aes(label=togene),  #不能用rownames,只能用数据框里面的某一列
    size=4.5,
    color="black",
    segment.color="black",show.legend=FALSE)+#添加关注的点的基因名
  labs(title = 'Volcano', x = 'Effect Size', y = '-Log10(pvalue)')
dev.off()
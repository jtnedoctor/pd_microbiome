setwd("~/project/16s/ndd/r_analysis")
#save(blautia_genes_cor_pd,blautia_genes_cor_pd_all,file = "blautia_genes_cor_pd.Rdata") #此处blautia_genes_cor_pd的r>0.3
load("blautia_genes_cor_pd.Rdata")
#save(res_rnaseq_metaallgenes,res_rnaseq_metaNAgenes,res_rnaseq_metasigdiffgenes.fixed,res_rnaseq_metasigdiffgenes.random,file = "res_rnaseq_metagenes.Rdata")
load("res_rnaseq_metagenes.Rdata")
#save(meta_ana_pd_chip_all,meta_ana_pd_chip_sig,file = "meta_chip_pd_unrep_fixed.Rdata")
load("meta_chip_pd_unrep_fixed.Rdata")

input_all <- list(Blautia_cor_genes=rownames(blautia_genes_cor_pd),DEGs_RNASeq=rownames(res_rnaseq_metasigdiffgenes.fixed),DEGs_Microarray=rownames(meta_ana_pd_chip_sig))

input_genes_rnaseq_chip <- list(RNASeq_genes=rownames(res_rnaseq_metaallgenes),Microarray_genes=rownames(meta_ana_pd_chip_all))

input_blautia_genes_rnaseq_chip <- list(Blautia_cor_genes=rownames(blautia_genes_cor_pd),RNASeq_genes=rownames(res_rnaseq_metaallgenes),Microarray_genes=rownames(meta_ana_pd_chip_all))
#方法三
#install_github("js229/Vennerable")
library(devtools)
library('Vennerable')
data1 <- Vennerable::Venn(input_all)
data2 <- Vennerable::Venn(input_genes_rnaseq_chip)
data3 <- Vennerable::Venn(input_blautia_genes_rnaseq_chip)
Weights(data1)
#plot(data1,doWeight=F,type='ellipses')
pdf(file = 'Venndiagram_blautia_rnaseq_chip_DEGs.pdf')
plot(data1,doWeight=F)
dev.off()

pdf(file = 'Venndiagram_genes_rnaseq_chip.pdf')
plot(data2,doWeight=F)
dev.off()

pdf(file = 'Venndiagram_blautia_genes_rnaseq_chip.pdf')
plot(data3,doWeight=F)
dev.off()
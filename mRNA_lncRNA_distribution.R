setwd("~/project/16s/ndd/r_analysis")
#save(blautia_genes_cor_pd,blautia_genes_cor_pd_all,file = "blautia_genes_cor_pd.Rdata") #此处blautia_genes_cor_pd的r>0.3
load("blautia_genes_cor_pd.Rdata")
#save(res_rnaseq_metaallgenes,res_rnaseq_metaNAgenes,res_rnaseq_metasigdiffgenes.fixed,res_rnaseq_metasigdiffgenes.random,file = "res_rnaseq_metagenes.Rdata")
load("res_rnaseq_metagenes.Rdata")
#save(meta_ana_pd_chip_all,meta_ana_pd_chip_sig,file = "meta_chip_pd_unrep_fixed.Rdata")
#load("meta_chip_pd_unrep_fixed.Rdata")
library(biomaRt)
mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset('hsapiens_gene_ensembl', mart)
annot <- getBM(
  mart = mart,
  attributes = c(
    'ensembl_gene_id',
    'hgnc_symbol',
    'entrezgene_id',
    'gene_biotype'),
  uniqueRows = TRUE)
#与blautia的r在0.3以上的基因的分布
#res_rnaseq_metaallgenes,res_rnaseq_metasigdiffgenes.fixed
rnatype <- res_rnaseq_metasigdiffgenes.fixed[rownames(res_rnaseq_metasigdiffgenes.fixed) %in% rownames(blautia_genes_cor_pd),]
#与blautia所有基因的基因分布
rnatype <- blautia_genes_cor_pd_all
rnatype$hgnc_symbol <- rownames(rnatype)
uniannot <- annot[!duplicated(annot$hgnc_symbol),]
rnatypeall <- merge(rnatype,uniannot,by="hgnc_symbol",all.x=T)
rnatypeall$gene_biotype[is.na(rnatypeall$gene_biotype)] <- "undefined"
table(rnatypeall$gene_biotype)
table(rnatypeall$gene_biotype,ifelse(rnatypeall$TE.fixed >0,"UP","DOWN"))
#blautia_genes_cor_pd_all_biotype <- rnatypeall
#save(blautia_genes_cor_pd_all_biotype,file = "blautia_genes_cor_pd_all_biotype.Rdata")
load("blautia_genes_cor_pd_all_biotype.Rdata")
#rnatype$gene_biotype <- uniannot$gene_biotype[ifelse(uniannot$hgnc_symbol %in% rnatype$genesymbol,T,NA)]
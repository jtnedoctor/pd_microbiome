#有差别Genus的相关性####
ps7_pd <- ps7_pd
ps7_G = tax_glom(ps7_pd, "Genus")
ps7 <- phyloseq::subset_samples(ps7_G,source %in% samplesource)
ps7 <- prune_samples(sample_sums(ps7) > 0,ps7)
ps7  = transform_sample_counts(ps7, function(x) x / sum(x))
ps7 <- prune_samples(sample_sums(ps7) > 0,ps7)
ps <- ps7
tmp_otu<-data.frame(ps@otu_table@.Data,taxo=ps@tax_table@.Data[,6]) #c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
tmp_otu<-aggregate(.~taxo,tmp_otu,sum)
rownames(tmp_otu) <- tmp_otu$taxo
otutable <- tmp_otu[,-1]
metadata <-as.data.frame(ps@sam_data@.Data)
colnames(metadata) <- ps@sam_data@names
rownames(metadata) <- ps@sam_data@row.names
metadata$sample <- rownames(metadata)
data <- otutable
design = metadata
#top100exp = data[rownames(sigDEG)[1:100],]
tdata <- t(data)
sig_data <- tdata[,paste0("Genus","_",colnames(tdata)) %in% rownames(assign(paste("res_16s_metasigdiffgenus.fixed",samplesource,sep="_"),res_16s_metasigdiffgenus.fixed))]
sig_data <- sig_data[rowSums(sig_data)!=0,]
#install.packages("corrplot")
library(corrplot)
cor_data <-  cor(sig_data,method = "spearman")
res1 <- cor.mtest(sig_data, conf.level = .95)
#corrplot(cor_data, p.mat = res1$p, sig.level = .05, addrect = 2)
#corrplot(cor_data, method = "color", order ="AOE", addCoef.col="grey")
#corrplot(cor_data, method="number", col="black", cl.pos="n")
filename=paste0("ndd_pd_all_otu_16s_corplot_Genus_metasig.",samplesource,".pdf")
pdf(filename, width = 10, height = 10)
#corrplot(cor_data, order = "hclust", addrect = 2,tl.cex=1)
corrplot(cor_data, order = "hclust", addrect = 2,tl.cex=1,p.mat = res1$p, insig = "blank") #不显示无统计差别的
dev.off()
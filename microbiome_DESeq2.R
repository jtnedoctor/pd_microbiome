setwd("~/project/16s/ndd/r_analysis")
#DESeq2组间丰度比较--不同研究同一OTU水平差异性统计
library(DESeq2)
#TCGA数据库认为是阳性存在的菌群
#The Cancer MicrobiomeAtlas (TCMA)
#tcma <- read.table("/home/shpc_6fbd50ad63/project/16s/ndd/r_analysis/bacteria.sample.relabund.genus.txt",sep = "\t",header = T,fill=T,row.names = 1)
#save(tcma,file = "tcma.Rdata")
load("tcma.Rdata")
#save(contam_all,file = "contam_all.Rdata")
load("contam_all.Rdata")   #要去除的潜在污染菌群
############
#save(ps7_G,file = "pd_ps7_G.Rdata") #PD的脑和血的菌群数据
load("pd_ps7_G.Rdata")
############
############
#save(ps7_G_otu_99_0.97,file = "ps7_G_otu_99_0.97.Rdata") #PD的粪便和血的菌群数据
# load("ps7_G_otu_99_0.97.Rdata")
# ps7_G <- ps7_G_otu_99_0.97
############
`%!in%` = Negate(`%in%`) #高级
taxlevel="Genus"
datasets="dataset_1"
samplesource="brain"  #blood,feces,brain
ps7 <- phyloseq::subset_samples(ps7_G,organ %in% samplesource)  #rnaseq为organ，16s为source
ps7 <- phyloseq::subset_samples(ps7,dataset %!in% datasets) #datasets
ps <- ps7
ps <- prune_taxa(taxa_sums(ps) > 0, ps) #为了保证targetgenus都在后面的研究中至少出现一次
#取出所有未污染genera
# feces_posi <- targetgenus #粪便中有的，为血液和脑提供过滤
# save(feces_posi,file = "feces_posi_genus.Rdata")
load("feces_posi_genus.Rdata")
targetgenus <- unique(ps@tax_table@.Data[,6]) 
targetgenus <- targetgenus[!targetgenus %in% contam_all]
targetgenus <- targetgenus[targetgenus %in% feces_posi]
metadata <-as.data.frame(ps@sam_data@.Data)
colnames(metadata) <- ps@sam_data@names
rownames(metadata) <- ps@sam_data@row.names
metadata$sample <- rownames(metadata)
prjlist <- as.vector(as.character(unique(metadata$project)))
#重复样本多次测序的去重样本id为sam_index
#save(genes_com,otus_com,sam_index,metadata,file = "common_genus_genes.pd.brain.Rdata")
#load("common_genus_genes.pd.brain.Rdata")  #重复的prjna557205
#save(sam_index,file = "sam_index.Rdata")
load("sam_index.Rdata")  #重复的prjna557205
#prjlist <- c("prjna391524")
#prjlist <- c("prjdb8639","prjeb22977","prjeb27564","prjeb30615","prjeb4927","prjna381395","prjna510730","prjna601994")
res_16s_prjs <- list() 
k=1
for (j in prjlist) {
  ps7 <- phyloseq::subset_samples(ps7_G,organ %in% samplesource)  #rnaseq为organ，16s为source
  ps7 <- phyloseq::subset_samples(ps7,project %in% j)
  ps <- ps7
  tmp_otu<-data.frame(ps@otu_table@.Data,taxo=ps@tax_table@.Data[,6]) #c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  tmp_otu<-aggregate(.~taxo,tmp_otu,sum)
  rownames(tmp_otu) <- tmp_otu$taxo
  tmp_otu <- tmp_otu[,-1]
  tmp_otu <- tmp_otu[!rownames(tmp_otu) %in% contam_all,]
  tmp_otu <- tmp_otu[!rowSums(tmp_otu)==0,]
  if (j == "prjna557205") {
    tmp_otu <- tmp_otu[,colnames(tmp_otu) %in% sam_index]
  }else{
    tmp_otu <- tmp_otu
  }
  ps@otu_table@.Data<-as.matrix(tmp_otu)
  ps@sam_data$group <- factor(ps@sam_data$group) #先因子化才行
  ps@sam_data$group <- relevel(ps@sam_data$group, ref = "Control") #前后顺序很重要
  #ps <- prune_taxa(taxa_sums(ps) > 0, ps) 
  ps <- prune_samples(sample_sums(ps) > 0,ps)
  psDE <- phyloseq_to_deseq2(ps, ~group)
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  geoMeans = apply(counts(psDE), 1, gm_mean)
  psDE = estimateSizeFactors(psDE, geoMeans = geoMeans)
  ddpsDE <- DESeq(psDE)
  reslist <- resultsNames(ddpsDE)
  alpha = 0.05
  res = results(ddpsDE,name = reslist[2])
  res = res[order(rownames(res), na.last=NA), ]
  #sigres = res[(res$padj < alpha), ]
  res <- as.data.frame(res)
  res$u95ci <- res$log2FoldChange + 1.96*res$lfcSE
  res$l95ci <- res$log2FoldChange - 1.96*res$lfcSE
  res_16s_prjs[[j]] <- res
  k=k+1
}

#assign(paste("res_16s_otu_97_prjs",samplesource,sep="_"),res_16s_prjs) 
#assign(paste("res_16s_otu_99_0.97_prjs_nocontam",samplesource,sep="_"),res_16s_prjs) 
assign(paste("res_rnaseq_prjs_nocontam",samplesource,sep="_"),res_16s_prjs) 
#assign(paste("res_16s_dada2_99_prjs",samplesource,sep="_"),res_16s_prjs) 

#meta分析####
#install.packages("meta")
library(meta)
#准备数据
#load("res_ndd.Rdata")
res_16s_allgenus <- list()
res_16s_metaallgenus <- data.frame()
#临时救急
#genes <- Reduce(intersect, list(rownames(res_16s_prjs[[1]]),rownames(res_16s_prjs[[2]]),rownames(res_16s_prjs[[3]]),rownames(res_16s_prjs[[4]]),rownames(res_16s_prjs[[5]])))
#targetgenus <- genes[genes %in% feces_posi]
#res_prjs
prjlist <- as.vector(as.character(unique(metadata$project)))
for (i in targetgenus){
  tmp <- data.frame()
  for (j in prjlist){
    tmp[toupper(j),1:2] <- as.vector(res_16s_prjs[[j]][i,])[2:3] #只要log2FoldChange和lfcSE
  }
  tmp <- na.omit(tmp)
  stulabs <- rownames(tmp)
  res_16s_allgenus[[paste0("Genus","_", i)]] <- tmp #此处添加的Genus命名list的数据框
  pfs=metagen(log2FoldChange,lfcSE,sm="SMD",data=res_16s_allgenus[[paste0("Genus","_",i)]],
              studlab=toupper(stulabs)) #studlab=paste("Genus",toupper(prjlist),sep="-")
  #p <- forest(pfs)
  res_16s_metaallgenus[paste0("Genus","_", i),1:12] <- c(pfs$TE.fixed,pfs$lower.fixed,pfs$upper.fixed,pfs$pval.fixed,pfs$TE.random,
                                                         pfs$lower.random,pfs$upper.random,pfs$pval.random,pfs$tau2,pfs$tau,pfs$I2,pfs$pval.Q)
}

colnames(res_16s_metaallgenus) <- c("TE.fixed","lower.fixed","upper.fixed","pval.fixed","TE.random",
                                    "lower.random","upper.random","pval.random","tau2","tau","I2","pval.Q")
alpha = 0.05
res_16s_metaallgenus <- res_16s_metaallgenus[,1:8]
res_16s_metaNAgenus <- na.omit(res_16s_metaallgenus)
res_16s_metasigdiffgenus.fixed <- res_16s_metaNAgenus[(res_16s_metaNAgenus$pval.fixed < alpha), ]
res_16s_metasigdiffgenus.random <- res_16s_metaNAgenus[(res_16s_metaNAgenus$pval.random < alpha), ]
####
#仅包含在feces存在的菌群
#16s的粪便和血的meta数据
assign(paste("res_16s_otu_99_0.97_metaallgenus_nocontam",samplesource,sep="_"),res_16s_metaallgenus) 
assign(paste("res_16s_otu_99_0.97_metasigdiffgenus.fixed_nocontam",samplesource,sep="_"),res_16s_metasigdiffgenus.fixed) 
assign(paste("res_16s_otu_99_0.97_metasigdiffgenus.random_nocontam",samplesource,sep="_"),res_16s_metasigdiffgenus.random) 
save(res_16s_otu_99_0.97_prjs_nocontam_feces,res_16s_otu_99_0.97_prjs_nocontam_blood,res_16s_otu_99_0.97_metaallgenus_nocontam_feces,res_16s_otu_99_0.97_metaallgenus_nocontam_blood,
     res_16s_otu_99_0.97_metasigdiffgenus.fixed_nocontam_feces,res_16s_otu_99_0.97_metasigdiffgenus.fixed_nocontam_blood,
     res_16s_otu_99_0.97_metasigdiffgenus.random_nocontam_feces,res_16s_otu_99_0.97_metasigdiffgenus.random_nocontam_blood,file = "res_16s_otu_99_0.97_fecesblood_nocontam.Rdata")
load("res_16s_otu_99_0.97_fecesblood_nocontam.Rdata")
write.table(file='res_16s_otu_99_0.97_metasigdiffgenus.fixed_nocontam_feces.xls',res_16s_otu_99_0.97_metasigdiffgenus.fixed_nocontam_feces,sep="\t",quote=F,row.names=T)
write.table(file='res_16s_otu_99_0.97_prjs_nocontam_blood.xls',res_16s_otu_99_0.97_prjs_nocontam_blood[[1]],sep="\t",quote=F,row.names=T)

#仅包含在feces存在的菌群
#Rna-seq的粪便和血的meta数据
assign(paste("res_rnaseq_metaallgenus_nocontam",samplesource,sep="_"),res_16s_metaallgenus) 
assign(paste("res_rnaseq_metasigdiffgenus.fixed_nocontam",samplesource,sep="_"),res_16s_metasigdiffgenus.fixed) 
assign(paste("res_rnaseq_metasigdiffgenus.random_nocontam",samplesource,sep="_"),res_16s_metasigdiffgenus.random) 
save(res_rnaseq_prjs_nocontam_brain,res_rnaseq_prjs_nocontam_blood,res_rnaseq_metaallgenus_nocontam_brain,res_rnaseq_metaallgenus_nocontam_blood,
     res_rnaseq_metasigdiffgenus.fixed_nocontam_brain,res_rnaseq_metasigdiffgenus.fixed_nocontam_blood,
     res_rnaseq_metasigdiffgenus.random_nocontam_brain,res_rnaseq_metasigdiffgenus.random_nocontam_blood,file = "res_rnaseq_fecesblood_nocontam.Rdata")
load("res_rnaseq_fecesblood_nocontam.Rdata")

write.table(file='res_rnaseq_metasigdiffgenus.fixed_nocontam_brain.xls',res_rnaseq_metasigdiffgenus.fixed_nocontam_brain,sep="\t",quote=F,row.names=T)


#三个血的数据放在一起meta####
#一个16s，两个rna-seq
res_16s_prjs <- list(prjna675864=res_rnaseq_prjs_nocontam_blood[[1]],prjna349023=res_rnaseq_prjs_nocontam_blood[[2]],prjna391524=res_16s_otu_99_0.97_prjs_nocontam_blood[[1]]) 
assign(paste("res_rnaseq_res_16s_otu_99_0.97_prjs_nocontam",samplesource,sep="_"),res_16s_prjs) 
#meta分析
#install.packages("meta")
library(meta)
#准备数据
#load("res_ndd.Rdata")
res_16s_allgenus <- list()
res_16s_metaallgenus <- data.frame()
#临时救急
genes <- Reduce(intersect, list(rownames(res_16s_prjs[[1]]),rownames(res_16s_prjs[[2]]),rownames(res_16s_prjs[[3]])))
targetgenus <- genes
#res_prjs
#prjlist <- as.vector(as.character(unique(metadata$project)))
prjlist <- as.vector(as.character(names(res_16s_prjs)))
for (i in targetgenus){
  tmp <- data.frame()
  for (j in prjlist){
    tmp[toupper(j),1:2] <- as.vector(res_16s_prjs[[j]][i,])[2:3] #只要log2FoldChange和lfcSE
  }
  tmp <- na.omit(tmp)
  stulabs <- rownames(tmp)
  res_16s_allgenus[[paste0("Genus","_", i)]] <- tmp #此处添加的Genus命名list的数据框
  pfs=metagen(log2FoldChange,lfcSE,sm="SMD",data=res_16s_allgenus[[paste0("Genus","_",i)]],
              studlab=toupper(stulabs)) #studlab=paste("Genus",toupper(prjlist),sep="-")
  #p <- forest(pfs)
  res_16s_metaallgenus[paste0("Genus","_", i),1:12] <- c(pfs$TE.fixed,pfs$lower.fixed,pfs$upper.fixed,pfs$pval.fixed,pfs$TE.random,
                                                         pfs$lower.random,pfs$upper.random,pfs$pval.random,pfs$tau2,pfs$tau,pfs$I2,pfs$pval.Q)
}

colnames(res_16s_metaallgenus) <- c("TE.fixed","lower.fixed","upper.fixed","pval.fixed","TE.random",
                                    "lower.random","upper.random","pval.random","tau2","tau","I2","pval.Q")
alpha = 0.05
res_16s_metaallgenus <- res_16s_metaallgenus[,1:8]
res_16s_metaNAgenus <- na.omit(res_16s_metaallgenus)
res_16s_metasigdiffgenus.fixed <- res_16s_metaNAgenus[(res_16s_metaNAgenus$pval.fixed < alpha), ]
res_16s_metasigdiffgenus.random <- res_16s_metaNAgenus[(res_16s_metaNAgenus$pval.random < alpha), ]


assign(paste("res_rnaseq_16s_otu_99_0.97_metaallgenus_nocontam",samplesource,sep="_"),res_16s_metaallgenus) 
assign(paste("res_rnaseq_16s_otu_99_0.97_metasigdiffgenus.fixed_nocontam",samplesource,sep="_"),res_16s_metasigdiffgenus.fixed) 
assign(paste("res_rnaseq_16s_otu_99_0.97_metasigdiffgenus.random_nocontam",samplesource,sep="_"),res_16s_metasigdiffgenus.random) 
save(res_rnaseq_res_16s_otu_99_0.97_prjs_nocontam_blood,res_rnaseq_16s_otu_99_0.97_metaallgenus_nocontam_blood,res_rnaseq_16s_otu_99_0.97_metasigdiffgenus.fixed_nocontam_blood,res_rnaseq_16s_otu_99_0.97_metasigdiffgenus.random_nocontam_blood,file = "res_rnaseq_16s_otu_99_0.97_blood_nocontam.Rdata")
load("res_rnaseq_16s_otu_99_0.97_blood_nocontam.Rdata")
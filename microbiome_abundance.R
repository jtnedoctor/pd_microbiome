#根据不同source(blood,brain,feces)来源生成新的ps对象用于作图####
#先整理数据
setwd("~/project/16s/ndd/r_analysis")
#save(contam_all,file = "contam_all.Rdata")
load("contam_all.Rdata")   #要去除的潜在污染菌群
#save(feces_posi,file = "feces_posi_genus.Rdata")
load("feces_posi_genus.Rdata")  #值分析在肠道中存在的菌群
#save(sam_index,file = "sam_index.Rdata")
load("sam_index.Rdata")  #脑菌群重复的prjna557205
############
#save(ps7_G,file = "pd_ps7_G.Rdata") #PD的脑和血的菌群数据
#load("pd_ps7_G.Rdata")
############
############
#save(ps7_G_otu_99_0.97,file = "ps7_G_otu_99_0.97.Rdata") #PD的粪便和血的菌群数据
load("ps7_G_otu_99_0.97.Rdata")
ps7_G <- ps7_G_otu_99_0.97
############
`%!in%` = Negate(`%in%`) #高级
taxlevel="Genus"
datasets="dataset_1"
samplesource="blood"  #blood,feces,brain
ps7 <- phyloseq::subset_samples(ps7_G,source %in% samplesource)  #rnaseq为organ，16s为source
ps7 <- phyloseq::subset_samples(ps7,dataset %!in% datasets) #datasets
#ps7 <- phyloseq::subset_samples(ps7,sampleid %in% sam_index)  #脑,brain菌群重复的prjna557205
#自定义%in%的否定表达
`%!in%` = Negate(`%in%`) #高级
ps7 <- phyloseq::subset_taxa(ps7,Genus %!in% contam_all) #此处的feces_posi已经是去除过污染菌群的genera
#ps7 <- phyloseq::subset_taxa(ps7,Genus %in% feces_posi) #此处的feces_posi已经是去除过污染菌群的genera
ps7 <- prune_samples(sample_sums(ps7) > 0,ps7)
ps7 <- prune_taxa(taxa_sums(ps7) > 0, ps7) #为了保证targetgenus都在后面的研究中至少出现一次
ps <- ps7
#Phylum层面
library(amplicon)
metadata <-as.data.frame(ps@sam_data@.Data)
colnames(metadata) <- ps@sam_data@names
rownames(metadata) <- ps@sam_data@row.names
metadata$sample <- rownames(metadata)
ps.phylum = tax_glom(ps, taxrank="Phylum", NArm=FALSE)
otutable.phylum <- as.data.frame(ps.phylum@otu_table@.Data)
taxa.phylum <- as.data.frame(ps.phylum@tax_table@.Data)
taxa.phylum$Phylum <- as.character(taxa.phylum$Phylum)
#taxa.phylum$Phylum[2] <- "Unassigned"
taxa.phylum <- taxa.phylum %>% mutate(Phylum = ifelse(grepl("Unassigned",Kingdom),"Unassigned",Phylum),
                                      Phylum = ifelse(is.na(Phylum), "Others", Phylum))
otutable.phylum$tmp <- taxa.phylum$Phylum
otutable.phylum <- aggregate(otutable.phylum[,1:length(colnames(otutable.phylum))-1],by=list(otutable.phylum$tmp),FUN=sum)
rownames(otutable.phylum) <- otutable.phylum[,1]
otutable.phylum_forall <- otutable.phylum
metadata_forall <- metadata
otutable.phylum <- otutable.phylum[,-1]
taxa.phylum <- as.matrix(taxa.phylum)
taxa.phylum_forall <- taxa.phylum
#单独作图#
#filename=paste0("plot_relabund_Phylum_tax.group.",samplesource,".pdf")
#pdf(filename)
#tax_stackplot(otutable.phylum,metadata,topN = 10,groupID = "group",style = "group",sorted = "abundance")
#dev.off()

assign(paste("otutable.phylum_16s_otu_99_0.97_nocontam",samplesource,sep="_"),otutable.phylum_forall)
assign(paste("metadata_16s_otu_99_0.97_nocontam",samplesource,sep="_"),metadata_forall)
assign(paste("taxa.phylum_16s_otu_99_0.97_nocontam",samplesource,sep="_"),taxa.phylum_forall)
# assign(paste("otutable.phylum_res_rnaseq_nocontam",samplesource,sep="_"),otutable.phylum_forall)
# assign(paste("metadata_res_rnaseq_nocontam",samplesource,sep="_"),metadata_forall)
# assign(paste("taxa.phylum_res_rnaseq_nocontam",samplesource,sep="_"),taxa.phylum_forall)

#汇总后的数据
#save(otutable.phylum_16s_otu_99_0.97_nocontam_feces,metadata_16s_otu_99_0.97_nocontam_feces,taxa.phylum_16s_otu_99_0.97_nocontam_feces,otutable.phylum_16s_otu_99_0.97_nocontam_blood,metadata_16s_otu_99_0.97_nocontam_blood,taxa.phylum_16s_otu_99_0.97_nocontam_blood,otutable.phylum_res_rnaseq_nocontam_brain,metadata_res_rnaseq_nocontam_brain,taxa.phylum_res_rnaseq_nocontam_brain,file = "res_rnaseq_and_otu_97_16s_for_allsource_plot_nocontam.Rdata")
#load("res_rnaseq_and_otu_97_16s_for_allsource_plot_nocontam_fecesposi.Rdata") #仅包含去污染且仅在肠道出现的菌群，直接加载此数据用于后续分析
load("res_rnaseq_and_otu_97_16s_for_allsource_plot_nocontam.Rdata") #仅去除污染菌属

otutable.phylumlist <- list(otutable.phylum_res_rnaseq_nocontam_brain,otutable.phylum_16s_otu_99_0.97_nocontam_feces,otutable.phylum_16s_otu_99_0.97_nocontam_blood)
otutable.phylum_all <- Reduce(function(x, y) merge(x, y, all=TRUE),otutable.phylumlist)
otutable.phylum_all[is.na(otutable.phylum_all)] <- 0
rownames(otutable.phylum_all) <- otutable.phylum_all[,1]
otutable.phylum_all <- otutable.phylum_all[,-1]
metadata_rnaseq_brain_1 <- metadata_res_rnaseq_nocontam_brain[,-1]
metadata_all <- rbind(metadata_rnaseq_brain_1,metadata_16s_otu_99_0.97_nocontam_feces,metadata_16s_otu_99_0.97_nocontam_blood)
#metadata_rnaseq_blood_1 <- metadata_rnaseq_blood[,-1] #包含blood的RNA-Seq数据
#metadata_all <- rbind(metadata_rnaseq_brain_1,metadata_rnaseq_blood_1,metadata_16s_otu_97_feces,metadata_16s_otu_97_blood)
metadata_all$source <- as.character(metadata_all$source)
metadata_all %>%
  mutate(methods=ifelse(grepl("RNA-Seq",region),"RNA-Seq","16S"),
         organs_group=paste0(organ,"_", diseasetype),
         methods_organs_group=paste0(methods,"_", organ,"_", diseasetype)) %>%
  select(source,diseasetype,group,gender,bmi,age,country,project,dataset,platform,avgSpotLen,
         region,DiseaseState,city,sample,organ,methods,organs_group,methods_organs_group) -> Sampledata
rownames(Sampledata) <- Sampledata$sample
#pdf("plot_relabund_Phylum_tax.blood.feces.brain_nocontam_fecesposi.pdf")
pdf("plot_relabund_Phylum_tax.blood.feces.brain_nocontam.pdf")
tax_stackplot(otutable.phylum_all,Sampledata,topN = 10,groupID = "methods_organs_group",style = "group",sorted = "abundance")
dev.off()
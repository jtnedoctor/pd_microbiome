setwd("~/project/16s/ndd/r_analysis")
library(RColorBrewer)
library(pheatmap)
library(ggplot2)

#save(res_rnaseq_prjs_nocontam_brain,res_rnaseq_prjs_nocontam_blood,res_rnaseq_metaallgenus_nocontam_brain,res_rnaseq_metaallgenus_nocontam_blood,res_rnaseq_metasigdiffgenus.fixed_nocontam_brain,res_rnaseq_metasigdiffgenus.fixed_nocontam_blood,res_rnaseq_metasigdiffgenus.random_nocontam_brain,res_rnaseq_metasigdiffgenus.random_nocontam_blood,file = "res_rnaseq_fecesblood_nocontam.Rdata")
load("res_rnaseq_fecesblood_nocontam.Rdata")

#save(res_16s_otu_99_0.97_prjs_nocontam_feces,res_16s_otu_99_0.97_prjs_nocontam_blood,res_16s_otu_99_0.97_metaallgenus_nocontam_feces,res_16s_otu_99_0.97_metaallgenus_nocontam_blood,res_16s_otu_99_0.97_metasigdiffgenus.fixed_nocontam_feces,res_16s_otu_99_0.97_metasigdiffgenus.fixed_nocontam_blood,res_16s_otu_99_0.97_metasigdiffgenus.random_nocontam_feces,res_16s_otu_99_0.97_metasigdiffgenus.random_nocontam_blood,file = "res_16s_otu_99_0.97_fecesblood_nocontam.Rdata")
load("res_16s_otu_99_0.97_fecesblood_nocontam.Rdata")
#仅包含在feces存在的菌群
union_genus <- c(rownames(res_16s_otu_99_0.97_metasigdiffgenus.fixed_nocontam_feces),rownames(res_16s_otu_99_0.97_metasigdiffgenus.fixed_nocontam_blood),rownames(res_rnaseq_metasigdiffgenus.fixed_nocontam_brain))
union_genus <- unique(union_genus)
feces=data.frame(genus_name=rownames(res_16s_otu_99_0.97_metaallgenus_nocontam_feces),Feces_ES=res_16s_otu_99_0.97_metaallgenus_nocontam_feces$TE.fixed,feces_pvalue=res_16s_otu_99_0.97_metaallgenus_nocontam_feces$pval.fixed)

blood=data.frame(genus_name=rownames(res_16s_otu_99_0.97_metaallgenus_nocontam_blood),Blood_ES=res_16s_otu_99_0.97_metaallgenus_nocontam_blood$TE.fixed,blood_pvalue=res_16s_otu_99_0.97_metaallgenus_nocontam_blood$pval.fixed)

brain=data.frame(genus_name=rownames(res_rnaseq_metasigdiffgenus.fixed_nocontam_brain),Brain_ES=res_rnaseq_metasigdiffgenus.fixed_nocontam_brain$TE.fixed,brain_pvalue=res_rnaseq_metasigdiffgenus.fixed_nocontam_brain$pval.fixed)
feces <- feces[feces$genus_name %in% union_genus,]
blood <- blood[blood$genus_name %in% union_genus,]
brain <- brain[brain$genus_name %in% union_genus,]
#NA赋值为-0.000001，为了后续和不存在的区分
feces[is.na(feces)] <- -0.000001
blood[is.na(blood)] <- -0.000001
brain[is.na(brain)] <- -0.000001

feces_blood_brain_list <- list(feces=feces,blood=blood,brain=brain)
merge_feces_blood_brain <- Reduce(function(x, y) merge(x, y, all=TRUE),feces_blood_brain_list)
merge_feces_blood_brain[is.na(merge_feces_blood_brain)] <- 0
# merge_taxo <- rbind(brain_taxo,feces_taxo,blood_taxo)
#save(merge_feces_blood_brain,merge_taxo,merge_taxo_uni,file = "merge_feces_blood_brain.Rdata")
load("merge_feces_blood_brain.Rdata")
merge_taxo_uni <- merge_taxo[!duplicated(merge_taxo$Genus),] #去除某一列重复的行
merge_taxo_uni$genus <- paste0("Genus","_", merge_taxo_uni$Genus)
merge_feces_blood_brain$Phylum <- merge_taxo_uni$Phylum[merge_taxo_uni$genus %in% merge_feces_blood_brain$genus_name]
merge_feces_blood_brain_nocontam <- merge_feces_blood_brain  #已经去除了污染菌群
comsig_genus <- merge_feces_blood_brain_nocontam
rownames(comsig_genus) <- comsig_genus[,1]
#row <- c(1,1,1,1,-0.250074656,-0.250074656)
comsig_genus$Phylum <- factor(comsig_genus$Phylum)
annotation_r <- data.frame(Phylum = comsig_genus$Phylum)
rownames(annotation_r) <- rownames(comsig_genus)

#ann_colors <- rainbow(length(unique(comsig_genus$Phylum)))
ann_colors <- list(Phylum = c( 'Actinobacteria' = "#550A46",'Bacteroidetes' = "#902044",
                               'Chloroflexi' = "#CB414B", 'Deinococcus-Thermus' = "#D16F7C",
                               'Euryarchaeota' =  "#E56B46", 'Firmicutes' ="#e59346",
                               'Planctomycetes' = "#F4A862", 'Proteobacteria' = "#F6DB86",
                               'Synergistetes' = "#e9f686", 'Verrucomicrobia' = "#DFE899"))
#breaks，把0设置为白色
bk <- c(seq(-4,-0.1,by=0.01),seq(0,4,by=0.01))
#merge_feces_blood_brain的NA值给赋值0，则能聚类
g <- pheatmap(as.data.frame(comsig_genus[,c(2,4,6)]),color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)), #把0设置为白色
              cluster_cols = FALSE,fontsize = 6,  cellwidth=50, cellheight = 7, border_color = NA,
              cluster_rows = TRUE,  scale = 'none',
              clustering_method = "complete",
              annotation_row  = annotation_r, treeheight_row = 150,annotation_colors = ann_colors,
              display_numbers= matrix(ifelse(comsig_genus[c(3,5,7)] < 0.01 & comsig_genus[c(3,5,7)] > 0, "**", ifelse(comsig_genus[c(3,5,7)] < 0.05 & comsig_genus[c(3,5,7)] > 0, "*",ifelse(comsig_genus[c(3,5,7)] >= 0.05 | comsig_genus[c(3,5,7)] == -0.000001,"ns",""))), ncol = 3),number_color="black",fontsize_number=6,breaks=bk)
ggsave(g, file = "merge_feces_blood_brain_diff_genus_nocontam_fposi.pdf", width = 20, height = 20, unit = 'cm')

#merge_feces_blood_brain的NA值不赋值0，则不能聚类
g <- pheatmap(as.data.frame(comsig_genus[,c(2,4,6)]),color = colorRampPalette(c("navy", "white", "firebrick3"))(150), 
              cluster_cols = FALSE,fontsize = 6,  cellwidth=50, cellheight = 7, border_color = NA,
              scale = 'none',
              cluster_rows = F,
              #clustering_method = "complete",
              annotation_row  = annotation_r, treeheight_row = 150,annotation_colors = ann_colors,
              display_numbers= matrix(ifelse(comsig_genus[c(3,5,7)] < 0.01 , "**", ifelse(comsig_genus[c(3,5,7)] < 0.05 , "*",ifelse(comsig_genus[c(3,5,7)] >= 0.05,"ns","X"))), 
                                      ncol = 3),number_color="black",fontsize_number=6)
ggsave(g, file = "merge_feces_blood_brain_diff_genus_nocontam_fposi.pdf", width = 20, height = 40, unit = 'cm')
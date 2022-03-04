setwd("~/project/16s/ndd/r_analysis")
#1.载入包
#BiocManager::install("SIAMCAT")
#BiocManager::install("mlr3")
library("SIAMCAT")
library(mlr3)
#2.载入数据
#先整理数据
setwd("/home/shpc_6fbd50ad63/project/16s/ndd/r_analysis")
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
samplesource="feces"  #blood,feces,brain
ps7 <- phyloseq::subset_samples(ps7_G,source %in% samplesource)  #rnaseq为organ，16s为source
ps7 <- phyloseq::subset_samples(ps7,dataset %!in% datasets) #datasets
#ps7 <- phyloseq::subset_samples(ps7,sampleid %in% sam_index)  #脑,brain菌群重复的prjna557205
#自定义%in%的否定表达
`%!in%` = Negate(`%in%`) #高级
ps7 <- phyloseq::subset_taxa(ps7,Genus %!in% contam_all) #此处的feces_posi已经是去除过污染菌群的genera
ps7 <- phyloseq::subset_taxa(ps7,Genus %in% feces_posi) #此处的feces_posi已经是去除过污染菌群的genera
ps7 <- prune_samples(sample_sums(ps7) > 0,ps7)
ps7 <- prune_taxa(taxa_sums(ps7) > 0, ps7) #为了保证targetgenus都在后面的研究中至少出现一次
ps <- ps7

ps  = transform_sample_counts(ps, function(x) x / sum(x))
tmp_otu<-data.frame(ps@otu_table@.Data,taxo=ps@tax_table@.Data[,6]) #c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
tmp_otu<-aggregate(.~taxo,tmp_otu,sum)
rownames(tmp_otu) <- tmp_otu$taxo
otutable <- tmp_otu[,-1]
metadata <-as.data.frame(ps@sam_data@.Data)
colnames(metadata) <- ps@sam_data@names
rownames(metadata) <- ps@sam_data@row.names
metadata$sample <- rownames(metadata)
#metadata <- metadata[,c("diseasetype","gender","bmi","age","country","project","platform","region", "organ")]
metadata <- metadata[,c("diseasetype","project","platform","region")]
metadata$project <- as.numeric(factor(metadata$project,labels = (1:length(unique(metadata$project)))))
metadata$platform <- as.numeric(factor(metadata$platform,labels = (1:length(unique(metadata$platform)))))
metadata$region <- as.numeric(factor(metadata$region,labels = (1:length(unique(metadata$region)))))
#metadata <- metadata[,c("diseasetype")]
metadata[metadata==""] <- NA
#定义健康组和正常对照组
#label.crc.zeller <- create.label(meta=meta.crc.zeller,label='Group', case='CRC')
#生成SIAMCAT需要的数据对象
#sc.obj <- siamcat(feat=feat.crc.zeller,label=label.crc.zeller,meta=meta.crc.zeller)
#或者
#sc.obj <- siamcat(phyloseq=ps,label='group', case='NDD')
sc.obj <- siamcat(feat=otutable,meta=metadata,label='diseasetype', case='PD')
#过滤table的feature数据
sc.obj <- filter.features(sc.obj,
                          filter.method = 'abundance',
                          cutoff = 0.001
)

#3.相关分析及作图
file_name=paste0("PD","_",taxlevel,"_",samplesource,"_","Association_Testing_test_nocontam.pdf")
#pdf(file=file_name,width = 20,height = 16)
check.associations(
  sc.obj,
  sort.by = 'fc',
  alpha = 0.05,
  max.show = 50,
  mult.corr = "fdr",
  detect.lim = 10 ^-6,
  plot.type = "quantile.box",    #bean,box,quantile.rect,quantile.box
  panels = c("fc", "prevalence", "auroc"),
  fn.plot=file_name)
#Yes
#dev.off()

#4.Confounder Testing混杂因素分析及作图
file_name=paste0("PD","_",taxlevel,"_",samplesource,"_","confounder_plots_nocontam.pdf")
check.confounders(
  sc.obj,
  fn.plot = file_name,
  meta.in = NULL,
  feature.type = 'filtered'
)

#5.Model Building建模分析
#数据标准化
sc.obj <- normalize.features(
  sc.obj,
  norm.method = "log.unit",
  norm.param = list(
    log.n0 = 1e-06,
    n.p = 2,
    norm.margin = 1
  )
)
#添加需要加入模型的meta变量
sc.obj <- add.meta.pred(sc.obj,pred.names=c('platform', 'region',"project"),std.meta=FALSE)

#Prepare Cross-Validation准备交叉验证
sc.obj <-  create.data.split(
  sc.obj,
  num.folds = 5,
  num.resample = 2
)
#Model Training训练模型
sc.obj <- train.model(
  sc.obj,
  perform.fs=T,
  method = "randomForest"      #"lasso", "enet", "ridge", "lasso_ll","ridge_ll", "randomForest"
)
# get information about the model type
model_type(sc.obj)
# access the models
models <- models(sc.obj)
models[[1]]
#Make Predictions预测
sc.obj <- make.predictions(sc.obj)
pred_matrix <- pred_matrix(sc.obj)
head(pred_matrix)
#Model Evaluation and Interpretation
sc.obj <-  evaluate.predictions(sc.obj)
#Evaluation Plot
file_name=paste0("PD","_",taxlevel,"_",samplesource,"_","model.evaluation.plot_nocontam.pdf")
#file_name=paste0("PD","_",taxlevel,"_",samplesource,"_","model.evaluation.plot_nocontam_lasso.pdf")
model.evaluation.plot(sc.obj,colours=c("red","green"),show.all=T,fn.plot=file_name)
#model.evaluation.plot('FR-CRC'=sc.obj.1,'CN-CRC'=sc.obj.2,colours=c('dimgrey', 'orange'),fn.plot="model.evaluation.plot.pdf") #同时画两个组的ROC曲线

#Interpretation Plot变量重要度排序
file_name=paste0("PD","_",taxlevel,"_",samplesource,"_","interpretation_nocontam.pdf")
#file_name=paste0("PD","_",taxlevel,"_",samplesource,"_","interpretation_nocontam_lasso.pdf")
model.interpretation.plot(
  sc.obj,
  fn.plot = file_name,
  consens.thres = 0.5,
  color.scheme = "PiYG",  #BrBG, PiYG, PRGn, PuOr, RdBu, RdGy, RdYlBu, RdYlGn,Spectral
  prompt=TRUE, 
  max.show = 30,
  verbose = 1,
  #norm.models = TRUE,
  limits = c(-3, 3),
  heatmap.type = 'zscore',
)
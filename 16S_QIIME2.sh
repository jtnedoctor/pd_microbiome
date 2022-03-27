#trim_galore
cd ~/project/16s/ndd/studies/prjna510730.pd.pair/rawseq/
dir='~/project/16s/ndd/studies/prjna510730.pd.pair/clean/'
for fq1 in *_1.fastq.gz
do
fq2=${fq1%%_1.fastq.gz}"_2.fastq.gz"
nohup trim_galore -q 20 -j 2 --phred33 --length 100 -e 0.1 --stringency 3 --paired --fastqc -o $dir $fq1 $fq2 &
done

#单端数据
cd ~/project/16s/ndd/studies/prjna510730.pd.pair/rawseq/tmp/
dir='~/project/16s/ndd/studies/prjna510730.pd.pair/clean/'
for fq in *fastq.gz
do
nohup trim_galore -q 20 -j 2 --phred33 --length 100 -e 0.1 --stringency 3 --fastqc -o $dir $fq &
done

#生成manifest文件
cd ~/project/16s/ndd/
mkdir -p manifests
cd ~/project/16s/ndd/manifests/
metapath="~/project/16s/ndd/metadata/"
ls ${metapath}|grep pair|cut -d"_" -f1 | while read id; do
  pwd="~/project/16s/ndd/studies/${id}/clean/"
  ls -lh $pwd|grep val_1.fq|awk -F" " '{print $NF}'>manifest1.tsv
  ls -lh $pwd|grep val_2.fq|awk -F" " '{print $NF}'>manifest2.tsv
#对于paired end文件，首行添加所有数据绝对路径
  sed -i "s#^#$pwd#g" manifest1.tsv
  sed -i "s#^#$pwd#g" manifest2.tsv
#提取metadata中的sample-id信息
  cut -d$'\t' -f1 ${metapath}${id}_metadata.tsv>sample-id.tsv
  sed -i '1d' sample-id.tsv
#分别生成forward和reverse的manifest文件
  paste -d"," sample-id.tsv manifest1.tsv manifest2.tsv > manifest3.tsv
#为三列添加一行为变量名
  awk  -F ','  'BEGIN {print "sample-id,forward-absolute-filepath,reverse-absolute-filepath"} {print $1","$2","$3}' manifest3.tsv |awk -F ',' ' {print $1"\t"$2"\t"$3}' > ${id}_manifest.tsv
  rm manifest1.tsv manifest2.tsv manifest3.tsv sample-id.tsv
done

#生成demux文件
cd ~/project/16s/ndd/
mkdir -p qiimefiles/qzvfiles
metapath="~/project/16s/ndd/metadata/"
manipath="~/project/16s/ndd/manifests/"
qiimepath="~/project/16s/ndd/qiimefiles/"
ls ${metapath}|grep pair|cut -d"_" -f1 | while read id; do
  #生成demux的qza文件
  qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ${manipath}${id}_manifest.tsv \
  --output-path ${qiimepath}${id}_end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2
   #summarize生成qzv文件
  qiime demux summarize \
  --i-data ${qiimepath}${id}_end-demux.qza \
  --o-visualization ${qiimepath}qzvfiles/${id}_end-demux.qzv  
done

#DADA2去噪文件
cd ~/project/16s/ndd/
mkdir -p qiimefiles/ASVfiles
qiimepath="~/project/16s/ndd/qiimefiles/"
metapath="~/project/16s/ndd/metadata/"
id="prjna391524.pd.pair"
#生成qza文件
  qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ${qiimepath}${id}_end-demux.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 295 \
  --p-trunc-len-r 220 \
  --o-table ${qiimepath}ASVfiles/${id}_table.qza \
  --o-representative-sequences ${qiimepath}ASVfiles/${id}_rep-seqs.qza \
  --o-denoising-stats ${qiimepath}ASVfiles/${id}_stats.qza \
  --p-n-threads 6  #cpu numbers
#生成qzv文件
  #table
  qiime feature-table summarize \
  --i-table ${qiimepath}ASVfiles/${id}_table.qza \
  --o-visualization ${qiimepath}qzvfiles/${id}_table.qzv \
  --m-sample-metadata-file ${metapath}${id}_metadata.tsv 
  #rep-seq
  qiime feature-table tabulate-seqs \
  --i-data ${qiimepath}ASVfiles/${id}_rep-seqs.qza \
  --o-visualization ${qiimepath}qzvfiles/${id}_rep-seqs.qzv
  #stat
  qiime metadata tabulate \
  --m-input-file ${qiimepath}ASVfiles/${id}_stats.qza \
  --o-visualization ${qiimepath}qzvfiles/${id}_stats.qzv


#物种注释,16s全序列
#vsearch 方法
qiime feature-classifier classify-consensus-vsearch \
  --i-query ${otuprjpath}ndd_merge_rep-seqs-cr-green-99.qza \
  --i-reference-reads ~/project/16s/microbioReference/green_99_otus.qza \
  --i-reference-taxonomy ~/project/16s/microbioReference/green_99_otus_taxonomy.qza \
  --p-perc-identity 0.9 \
  --p-query-cov 0.9 \
  --p-threads 10 \
  --o-classification ${otuprjpath}ndd_merge_taxonomy-green-99.qza


#构建进化树
cd ~/project/16s/ndd/
mkdir -p qiimefiles/otufiles/prjs
otupath="~/project/16s/ndd/qiimefiles/otufiles/"
otuprjpath="~/project/16s/ndd/qiimefiles/otufiles/prjs/"
asvpath="~/project/16s/ndd/qiimefiles/ASVfiles/"
metapath="~/project/16s/ndd/metadata/"
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ${otuprjpath}ndd_merge_rep-seqs-cr-green-99.qza \
  --o-alignment ${otuprjpath}ndd_merge_aligned-rep-seqs-halfnum-thirdtotal-ten.qza \
  --o-masked-alignment ${otuprjpath}ndd_merge_masked-aligned-rep-seqs.qza \
  --o-tree ${otuprjpath}ndd_merge_unrooted-tree.qza \
  --o-rooted-tree ${otuprjpath}ndd_merge_rooted-tree.qza \
  --p-n-threads 10





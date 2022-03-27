cd ~/project/16s/ndd/
mkdir -p ~/project/16s/ndd/revision_pd_picrust
revisionpath="~/project/16s/ndd/revision_pd_picrust/"
otuprjpath="~/project/16s/ndd/qiimefiles/otufiles/prjs/"
#生成dna-sequences.fasta和feature-table.biom文件
mkdir -p ${revisionpath}/picrust
mkdir -p ${revisionpath}/ftable
unzip -q ${otuprjpath}/ndd_merge_table-cr-green-99.qza -d ${revisionpath}/ftable/
tablepath="${revisionpath}/ftable/*/data/"
mkdir -p ${revisionpath}/dnaseq
unzip -q ${otuprjpath}/ndd_merge_rep-seqs-cr-green-99.qza -d ${revisionpath}/dnaseq/
dnapath="${revisionpath}/dnaseq/*/data/"
#正式数据
#unstratified output
time picrust2_pipeline.py -s ${dnapath}dna-sequences.fasta -i ${tablepath}feature-table.biom -o ${revisionpath}/picrust/picrust2_out_pipeline -p 8
#pathway stratified output based on per-sequence contributions
#picrust2_pipeline.py -s ${dnapath}dna-sequences.fasta -i ${tablepath}feature-table.biom -o ${otuprjpath}/picrust/picrust2_out_pipeline -p 6 --stratified --per_sequence_contrib
#dna-sequences.fasta来自rep-seqs.qza,而feature-table.biom来自table.qza
#输出结果均在picrust2_out_pipeline文件夹
#在第一列后面添加一列注释，结果表也方便在STAMP中进行差异比较。
#添加EC的注释
add_descriptions.py -i ${revisionpath}/picrust/picrust2_out_pipeline/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
  -o ${revisionpath}/picrust/picrust2_out_pipeline/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz
#KO添加注释
add_descriptions.py -i ${revisionpath}/picrust/picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO \
  -o ${revisionpath}/picrust/picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz
#pathway添加注释
add_descriptions.py -i ${revisionpath}/picrust/picrust2_out_pipeline/pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
  -o ${revisionpath}/picrust/picrust2_out_pipeline/pathways_out/path_abun_unstrat_descrip.tsv.gz

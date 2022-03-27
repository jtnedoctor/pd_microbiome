#pairend测序
ids=projectid
ls ~/project/16s/ndd/kraken2/studies/${ids}/rawseq/*gz|cut -d"/" -f 11|cut -d"_" -f 1|sort -u |while read id; do
  kraken2 --db ~/project/16s/microbioReference/minikraken_8GB_20200312  --threads 6  --report ${dir}/${ids}/${id}.test.report --output ${dir}/${ids}/${id}.test.output  --paired --gzip-compressed ~/project/16s/ndd/kraken2/studies/${ids}/rawseq/${id}_1.fastq.gz  ~/project/16s/ndd/kraken2/studies/${ids}/rawseq/${id}_2.fastq.gz
  bracken -d ~/project/16s/microbioReference/minikraken_8GB_20200312 -i ${dir}/${ids}/${id}.test.report -o ${dir}/${ids}/${id}.test.s.bracken -w ${dir}/${ids}/${id}.test.s.bracken.report -r 50 -l S
  kreport2mpa.py -r ${dir}/${ids}/${id}.test.s.bracken.report  -o ${dir}/${ids}/${id}.test.new.report
  #生成含各个level的数据框
  cat ${dir}/${ids}/${id}.test.new.report|grep d_Bacteria|grep p_|grep c_|grep o_|grep f_|grep g_|grep s_|sed 's/|/\t/g'|sed 's/[a-z]_//g'\
  |sed 's/\t/,/g'|awk  -F ','  'BEGIN {print "Kingdom,Phylum,Class,Order,Family,Genus,Species,samid"} {print $1","$2","$3","$4","$5","$6","$7","$8}'\
  |awk -F ',' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}'> ${dir}/${ids}/${id}_kraken2.txt
  sed -i "1s#samid#$id#g" ${dir}/${ids}/${id}_kraken2.txt; done
  rm ~/project/16s/ndd/kraken2/studies/${ids}/rawseq/*fastq.gz
  rm ${dir}/${ids}/*output
  let from=${from}+2
  let to=${to}+2
  done
done
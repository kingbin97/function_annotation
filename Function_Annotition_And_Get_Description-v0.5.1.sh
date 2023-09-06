#!/bin/bash
# 使用说明
usage() {
  echo "Usage: Scripts <-i input fasta> <-t thread> "
  echo "       -i          input pep"
  echo "       -t          -thread [Default:20]"
}
if [ $# -eq 0 ] || [ "$1" = "-h" ]; then
  usage
  exit 0
elif [ $# -lt 2 ]; then
  echo "Error: Missing required arguments."
  usage
  exit 1
fi
#### 传参
while getopts ":i:t:" opt; do
  case $opt in
    i) INPUT="$OPTARG";;
    t) THREADERS="$OPTARG";;
    \?) echo "Invalid option -$OPTARG" >&2;;
  esac
done
###################################################################################################################
pep_prefix=${INPUT/.*/}
###################################################################################################################
# description文件，这里无需修改
nr_description=/pub/Databases/00.Gene_Description/NR_Description.tab
uniprto_description=/pub/Databases/00.Gene_Description/UniprotID_Description.tab
GO_idmapping_db=/pub/Databases/00.Gene_Description/GOID_UniprotID.tab
GO_description=/pub/Databases/00.Gene_Description/GOID_Description.tab
COG_description=/pub/Databases/00.Gene_Description/COG_Description.tab
KOG_description=/pub/Databases/00.Gene_Description/KOG_Description.tab
KEGG_description=/pub/Databases/00.Gene_Description/kegg_total_map.tab
# scripts脚本
EXTRACT=/home/dahome/students/linjinbin/script/exract_description_from_results.py
GO_split=/home/dahome/students/linjinbin/script/GO_split.py
MERGE_ROWS=/home/dahome/students/linjinbin/tem_annotation/merge_rows_by_gene.py
INTERPRO=/pub/Databases/07.InterProscan/interproscan-5.61-93.0/interproscan.sh
##############DB##################################################################################################
NR_DB=/pub/Databases/01.NR/nr.dmnd
UNIPROT_DB=/pub/Databases/03.Uniprot_Swiss/uniprot_sprot.dmnd
COG_DB=/pub/Databases/05.KOG_COG/2022-07-11/COG2020/tem_COG.fa.dmnd
KOG_DB=/pub/Databases/05.KOG_COG/2022-07-11/KOG/KOG.dmnd
EGGNOG_DB=/pub/Databases/06.EggNOG
########################################################################################################################################################################################################################################################
# run map
########################################################################################################################################################################################################################################################
function run_nr_map() {
  echo "NR map"
  if ! command -v diamond &> /dev/null; then
    echo "未找到 diamond  命令，请先安装 diamond ."
    exit 1
  fi
  # 继续执行
  echo "已找到 diamond 命令，继续执行"
  ###################################
  diamond blastp --query ${INPUT} --db $NR_DB --max-target-seqs 1 --max-hsps 1 --evalue 1e-5 --very-sensitive --outfmt 6 --threads ${THREADERS} --memory-limit 64G --quiet --out ${pep_prefix}_nr_diamond.xls # 比对
  awk 'BEGIN{FS="\t";OFS="\t"}{print $1"\t"$2}' ${pep_prefix}_nr_diamond.xls | sed '1i Gene_ID\tnr_ID' > ${pep_prefix}_nr_diamond_geneid.list # 获取基因列表
  python ${EXTRACT} ${pep_prefix}_nr_diamond_geneid.list ${nr_description} ${pep_prefix}_nr_extract_description.tem # database使用表格处理效率太低，因此专门用这个脚本提取比对到的 nr 的基因的 description
  sort ${pep_prefix}_nr_extract_description.tem |uniq > ${pep_prefix}_nr_extract_description.tem.uniq # 去重
  csvtk -t join -f "nr_ID;nr_ID" ${pep_prefix}_nr_diamond_geneid.list ${pep_prefix}_nr_extract_description.tem.uniq --na NA --left-join > ${pep_prefix}_nr_description.results # join
  rm ${pep_prefix}_nr_extract_description.tem ${pep_prefix}_nr_extract_description.tem.uniq
  python $MERGE_ROWS ${pep_prefix}_nr_description.results ${pep_prefix}_nr_description.merge.results
}
########################################################################################################################################################################################################################################################
function run_swiss_prot_map() {
  echo "SwissProt map"
  if ! command -v diamond &> /dev/null; then
    echo "未找到 diamond  命令，请先安装 diamond ."
    exit 1
  fi
  # 继续执行
  echo "已找到 diamond 命令，继续执行"
  ###################################
  diamond blastp --query ${INPUT} --db $UNIPROT_DB --max-target-seqs 1 --max-hsps 1 --evalue 1e-5 --very-sensitive --outfmt 6 --threads ${THREADERS} --memory-limit 64G --quiet --out ${pep_prefix}_uniprot_diamond.xls # 比对
  awk 'BEGIN{FS="\t";OFS="\t"}{print $1"\t"$2}' ${pep_prefix}_uniprot_diamond.xls | sed 's/sp|//' | sed 's/|/\t/' | awk '{print $1"\t"$2}' | sort | uniq | sed '1i Gene_ID\tUniprot_ID' > ${pep_prefix}_uniprot_diamond_geneid.list # 提取基因信息
  csvtk -t join -f "Uniprot_ID;Uniprot_ID" ${pep_prefix}_uniprot_diamond_geneid.list ${uniprto_description} --na NA --left-join | uniq > ${pep_prefix}_uniprot_description.results
  python $MERGE_ROWS ${pep_prefix}_uniprot_description.results ${pep_prefix}_uniprot_description.merge.results
  #
  python $EXTRACT ${pep_prefix}_uniprot_diamond_geneid.list ${GO_idmapping_db} ${pep_prefix}_GOIDs_of_uniprot.txt #也可以用csvtk -t join
  python $GO_split ${pep_prefix}_GOIDs_of_uniprot.txt ${pep_prefix}_GOIDs_of_uniprot.split.txt
  csvtk -t join -f "GO_ID;GO_ID" ${pep_prefix}_GOIDs_of_uniprot.split.txt ${GO_description} --na NA --left-join | uniq > ${pep_prefix}_GO_uniprot_description.split.results
  python $MERGE_ROWS ${pep_prefix}_GO_uniprot_description.split.results ${pep_prefix}_GO_uniprot_description.split.merge_GO.results
}
########################################################################################################################################################################################################################################################
function run_COG_map() {
  echo "COG map"
  if ! command -v diamond &> /dev/null; then
    echo "未找到 diamond  命令，请先安装 diamond ."
    exit 1
  fi
  # 继续执行
  echo "已找到 diamond 命令，继续执行"
  diamond blastp --query ${INPUT} --db $COG_DB --max-target-seqs 1 --max-hsps 1 --evalue 1e-5 --very-sensitive --outfmt 6 --threads ${THREADERS} --memory-limit 64G --quiet --out ${pep_prefix}_COG_diamond.xls # 比对
  awk 'BEGIN{FS="\t";OFS="\t"}{print $1"\t"$2}' ${pep_prefix}_COG_diamond.xls | sort | uniq | sed '1i Gene_ID\tCOG_Protein' > ${pep_prefix}_COG_diamond_geneid.list # 提取基因信息
  csvtk -t join -f "COG_Protein;COG_Protein" ${pep_prefix}_COG_diamond_geneid.list ${COG_description} --na NA --left-join | uniq > ${pep_prefix}_COG_description.txt # 一个基因一行对应一个COG
  python $MERGE_ROWS ${pep_prefix}_COG_description.txt ${pep_prefix}_COG_description.merge.results # 一个基因一行对应合并后的COG
}
########################################################################################################################################################################################################################################################
function run_KOG_map() {
  echo "KOG map"
  if ! command -v diamond &> /dev/null; then
    echo "未找到 diamond  命令，请先安装 diamond ."
    exit 1
  fi
  # 继续执行
  echo "已找到 diamond 命令，继续执行"
  diamond blastp --query ${INPUT} --db $KOG_DB --max-target-seqs 1 --max-hsps 1 --evalue 1e-5 --very-sensitive --outfmt 6 --threads ${THREADERS} --memory-limit 64G --quiet --out ${pep_prefix}_KOG_diamond.xls # 比对
  awk 'BEGIN{FS="\t";OFS="\t"}{print $1"\t"$2}' ${pep_prefix}_KOG_diamond.xls | sort | uniq | sed '1i Gene_ID\tKOG_Protein' > ${pep_prefix}_KOG_diamond_geneid.list
  csvtk -t join -f "KOG_Protein;KOG_Protein" ${pep_prefix}_KOG_diamond_geneid.list ${KOG_description} --na NA --left-join |uniq > ${pep_prefix}_KOG_description.results
  python $MERGE_ROWS ${pep_prefix}_KOG_description.results ${pep_prefix}_KOG_description.merge.results
}
########################################################################################################################################################################################################################################################
function run_interproscan_map() {
  echo "InterproScan map, include Pfam "
  echo "检测 Java 版本 "
  java_version=$(java -version 2>&1 | awk -F '"' '/version/ {print $2}')
  # 提取主要版本号
  major_version=$(echo "$java_version" | awk -F '.' '{print $1}')
  if [[ $major_version -lt 11 ]]; then # 检查版本是否小于11
    echo "Java版本小于1.11，退出"
    exit 1
  fi
    echo "Java版本大于等于11，继续执行"
    bash $INTERPRO -i ${INPUT} --iprlookup --highmem --cpu ${THREADERS} --goterms --pathways -f tsv -o ${pep_prefix}_interproscan.xls >interproscan.log
    # IPR
    grep IPR ${pep_prefix}_interproscan.xls | awk 'BEGIN{FS="\t";OFS="\t"}{print $1"\t"$12"\t"$13}' | sort | uniq | sed '1i Gene_ID\tIPR_ID\tIPR_description' > ${pep_prefix}_IPR_interproscan.xls
    python $MERGE_ROWS ${pep_prefix}_IPR_interproscan.xls ${pep_prefix}_IPR_interproscan.merge.results
    # GO
    grep "GO:" ${pep_prefix}_interproscan.xls | awk '!seen[$0]++' | cut -f1,4,14 | sed '1i Gene_ID\tDB\t_GO_ID' | sed 's/|/;/g' > ${pep_prefix}_GOIDs_of_interproscan.txt
    python $GO_split ${pep_prefix}_GOIDs_of_interproscan.txt ${pep_prefix}_GOIDs_of_interproscan.split.tem
    sed '1d' ${pep_prefix}_GOIDs_of_interproscan.split.tem | sort | uniq | sed '1i Gene_ID\tGO_ID' > ${pep_prefix}_GOIDs_of_interproscan.split.txt
    rm ${pep_prefix}_GOIDs_of_interproscan.split.tem
    csvtk -t join -f "GO_ID;GO_ID" ${pep_prefix}_GOIDs_of_interproscan.split.txt ${GO_description} --na NA --left-join | uniq > ${pep_prefix}_GO_Interproscan_description.results
    python $MERGE_ROWS ${pep_prefix}_GO_Interproscan_description.results ${pep_prefix}_GO_Interproscan_description.merge.results
    # 其他
    cut -f4 ${pep_prefix}_interproscan.xls | sort | uniq | while read i
      do
      awk -v var="$i" '{ if ($4 == var) { print } }' ${pep_prefix}_interproscan.xls | awk 'BEGIN{FS="\t";OFS="\t"}{print $1"\t"$5"\t"$6}' | sort | uniq | sed "1i Gene_ID\t${i}_ID\t${i}_description" > ${pep_prefix}_${i}_interproscan.xls
      python $MERGE_ROWS ${pep_prefix}_${i}_interproscan.xls ${pep_prefix}_${i}_interproscan.merge.results
    done
}
########################################################################################################################################################################################################################################################
function run_eggnog_map() {
  echo "EggNOG map"
  if ! command -v emapper.py &> /dev/null; then
    echo "未找到 emapper.py  命令，请先安装 emapper.py ."
    exit 1
  fi
  # 继续执行
  echo "已找到 emapper.py 命令，继续执行"
  emapper.py -i ${INPUT} --cpu ${THREADERS} -o ${pep_prefix}_EggNOG --data_dir $EGGNOG_DB --override -m diamond --dmnd_ignore_warnings --dmnd_algo ctg --evalue 0.00001 --score 60 --pident 40 --query_cover 30 --subject_cover 30 --tax_scope auto --target_orthologs one2one --go_evidence non-electronic --pfam_realign none --report_orthologs --decorate_gff yes --go_evidence experimental --excel >emapper.out 2> emapper.err
  mkdir eggnog && mv ${pep_prefix}_EggNOG* eggnog/. 
  #GO
  grep "GO:" eggnog/${pep_prefix}_EggNOG.emapper.annotations | cut -f1,5,10 | sed 's/\,GO/\;GO/g' | sed '1i Gene_ID\tEGGnog_Ortho\tGO_ID' > ${pep_prefix}_GOIDs_of_EggNOG.txt
  python $GO_split ${pep_prefix}_GOIDs_of_EggNOG.txt   ${pep_prefix}_GOIDs_of_EggNOG.split.txt
  csvtk -t join -f "GO_ID;GO_ID" ${pep_prefix}_GOIDs_of_EggNOG.split.txt ${GO_description} --na NA --left-join | uniq > ${pep_prefix}_EggNOG_GO_description.results
  python $MERGE_ROWS ${pep_prefix}_EggNOG_GO_description.results ${pep_prefix}_EggNOG_GO_description.merge.results
}
########################################################################################################################################################################################################################################################
function run_kofam_map() {
  echo "KEGG map"
  if ! command -v kofamscan &> /dev/null; then
    echo "未找到 kofamscan  命令，请先安装 kofamscan ."
    exit 1
  fi
  # 继续执行
  echo "已找到 kofamscan.py 命令，继续执行"
  kofamscan --cpu ${THREADERS} --format mapper -e 1e-5 ${INPUT} -o ${pep_prefix}_kofam.xls   # map
  grep -P '\tK' ${pep_prefix}_kofam.xls | sort | uniq | sed '1i Gene_ID\tKO_ID' > ${pep_prefix}_kofam.txt   # description
  csvtk -t join -f "KO_ID;KO_ID" ${pep_prefix}_kofam.txt ${KEGG_description} --left-join --na "-" > ${pep_prefix}_kofam.description.results
  python $MERGE_ROWS ${pep_prefix}_kofam.description.results ${pep_prefix}_kofam.description.merge.results
}
############################################################################################################
# run map
# 有不需要比对的，如eggnog，kegg，nr比较慢的，可以直接注释不跑
# eg: #echo "NR map"
# eg: #run_nr_map()
############################################################################################################
if [[ -e ${INPUT} ]]; then
  grep ">" ${INPUT} | sed 's/>//' | sed 's/ /\t/' | cut -f1 |sed '1i Gene_ID' > ${pep_prefix}_PEP_GeneID.list.txt
  #run_swiss_prot_map
  #run_KOG_map
  #run_COG_map
  #run_interproscan_map
  #run_eggnog_map
  #run_kofam_map
  run_nr_map
else
  echo "请检查输入文件"
fi
############################################################################################################
# merge total function annotations 
############################################################################################################
# 或者直接采用以下脚本
############################################################################################################
python $MERGE_ANNOTATION

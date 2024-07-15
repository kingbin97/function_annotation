#!/bin/bash
# 使用说明
usage() {
  echo "Usage: Scripts <-i input fasta> [-t threads] [-m sensitive mode]"
  echo "       -i          input pep"
  echo "       -t          thread  [Default:40]"
  echo "       -s          species [Default:osa]"
  echo "       -m          sensitive mode  fast | mid-sensitive | sensitive | more-sensitive | very-sensitive | ultra-sensitive [Default:very-sensitive]"
  echo -e "\e[34m                                      _  __  _                   _       _          \e[0m"
  echo -e "\e[34m                                     | |/ / (_)                 | |     (_)        \e[0m"
  echo -e "\e[34m                                     | ' /   _   _ __     __ _  | |__    _   _ __  \e[0m"
  echo -e "\e[34m                                     |  <   | | | '_ \   / _  | | '_ \  | | | '_ \ \e[0m"
  echo -e "\e[34m                                     | . \  | | | | | | | (_| | | |_) | | | | | | |\e[0m"
  echo -e "\e[34m                                     |_|\_\ |_| |_| |_|  \__, | |_.__/  |_| |_| |_|\e[0m"
  echo -e "\e[34m                                                          __/ |                    \e[0m"
  echo -e "\e[34m                                                         |___/                     \e[0m"
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
while getopts ":i:t:s:m:" opt; do
  case $opt in
    i) INPUT="$OPTARG"
     ;;
    t) THREADS="$OPTARG"
     ;;
    s) SPECIES="$OPTARG"
     ;;
    m) SENSITIVE_MODE="$OPTARG"
     ;;
    \?) echo "Invalid option -$OPTARG" >&2
     ;;
  esac
done
###################################################################################################################
  if [[ -z $SPECIES ]]; then SPECIES=osa; fi
  if [[ -z $THREADS ]]; then THREADS=40; fi
  if [[ -z $SENSITIVE_MODE ]]; then SENSITIVE_MODE=very-sensitive; fi
###################################################################################################################
#pep_prefix=${INPUT/.*/}
  pep_file=$(basename "$INPUT")
  pep_prefix=${pep_file%.*}
###################################################################################################################
#######################################
# scripts脚本
#######################################
ITAK=/home/linjinbin/script/function_annotation/iTAK/iTAK.pl
EXTRACT=/home/linjinbin/script/function_annotation/exract_description_from_results.py
GO_split=/home/linjinbin/script/function_annotation/GO_split.py
MERGE_ROWS=/home/linjinbin/script/function_annotation/merge_rows_by_gene.py
INTERPRO=/pub/Databases/07.InterProscan/interproscan-5.65-97.0/interproscan.sh
GO_anno_from_inter=/home/linjinbin/script/function_annotation/GO_KEGG_annotation/GO_anno_from_tab.py
KEGG_anno_from_kofam=/home/linjinbin/script/function_annotation/GO_KEGG_annotation/kaas_kofam2pathwayAnalysis.py
CONT_GENE=/home/linjinbin/script/function_annotation/count_genes.py
#######################################
# Description文件，这里无需修改
#######################################
nr_description=/pub/Databases/00.Gene_Description/NR_Description.tab
uniprto_description=/pub/Databases/00.Gene_Description/UniprotID_Description.tab
GO_idmapping_db=/pub/Databases/00.Gene_Description/GOID_UniprotID.tab
GO_description=/pub/Databases/00.Gene_Description/GOID_Description.tab
COG_description=/pub/Databases/00.Gene_Description/COG_Description.tab
KOG_description=/pub/Databases/00.Gene_Description/KOG_Description.tab
KEGG_description=/pub/Databases/00.Gene_Description/kegg_total_map.tab
OS_description=/pub/Databases/00.Gene_Description/Os_Description.tab
ATH_description=/pub/Databases/00.Gene_Description/Ath_Description.tab
KO_KEGG_DESCRIPTION=/pub/Databases/00.Gene_Description/KO_Gene_KEGG_Pathway_Description.tab
PLANT_TFDB_DESCRIPTION=/pub/Databases/17.TF/2023-12-19/TF_Description.tab
#######################################
# DataBase                             
#######################################
NR_DB=/pub/Databases/01.NR/nr.dmnd
UNIPROT_DB=/pub/Databases/03.Uniprot_Swiss/uniprot_sprot.dmnd
COG_DB=/pub/Databases/05.KOG_COG/2022-07-11/COG2020/tem_COG.fa.dmnd
KOG_DB=/pub/Databases/05.KOG_COG/2022-07-11/KOG/KOG.dmnd
EGGNOG_DB=/pub/Databases/06.EggNOG
KO_PROFILE=/pub/Databases/08.kofam/2023-10-03/profiles
KO_LIST=/pub/Databases/08.kofam/2023-10-03/ko_list
OS_DB=/pub/Databases/15.Os/Os_IRGSP.dmnd
ATH_DB=/pub/Databases/16.Ath/Ath_TAIR11.dmnd
PLANT_TF_DB=/pub/Databases/17.TF/PlantTFDB.dmnd
###################################################################################################################
# soft check function                                                                                             #
###################################################################################################################
function check_soft() {
  ##########################
  # diamond
  ##########################
  if ! command -v diamond &> /dev/null; then
    echo "未找到软件: diamond , 请先安装 diamond ."
    exit 1
  else
    echo "已找到软件: diamond "
  fi
  ##########################
  ## Java
  ##########################
  java_version=$(java -version 2>&1 | awk -F '"' '/version/ {print $2}')
  major_version=$(echo "$java_version" | awk -F '.' '{print $1}') # 提取主要版本号
  if [[ $major_version -lt 11 ]]; then # 检查版本是否小于11
    echo "Java版本小于1.11，退出"
    exit 1
  else
    echo "Java版本大于等于11"
  fi
  ##########################
  # eggnog
  ##########################
  if ! command -v emapper.py &> /dev/null; then
    echo "未找到软件: emapper.py , 请先安装 emapper.py ."
    exit 1
  else
    echo "已找到软件: emapper.py "
  fi
  ##########################
  # kofam
  ##########################
  if ! command -v kofamscan &> /dev/null; then
    echo "未找到软件: kofamscan , 请先安装 kofamscan ."
    exit 1
  else
     echo "已找到软件: kofamscan "
  fi
  ##########################
  # csvtk
  ##########################
  if ! command -v csvtk &> /dev/null; then
    echo "未找到软件: csvtk , 请先安装 csvtk ."
    exit 1
  else
     echo "已找到软件: csvtk "
  fi
  ##########################
  # seqkit
  ##########################
  if ! command -v seqkit &> /dev/null; then
    echo "未找到软件: seqkit , 请先安装 seqkit ."
    exit 1
  else
     echo "已找到软件: seqkit "
  fi
}
##############################################################
# run map function                                           #
##############################################################
###############################
# NR                          #
###############################
function run_nr_map() {
  echo -e "\e[33m ################################\e[0m"
  echo -e "\e[33m ## run NR map                 ##\e[0m"
  echo -e "\e[33m ################################\e[0m"
  # 比对
  diamond blastp \
    --query $INPUT \
    --db $NR_DB \
    --max-target-seqs 1 \
    --max-hsps 1 \
    --evalue 1e-5 \
    --${SENSITIVE_MODE} \
    --outfmt 6 \
    --threads $THREADS \
    --quiet \
    --out ${pep_prefix}_nr_diamond.xls # 比对
  # 统计一下score的数量分布
    count100=$(awk '$12 > 100 { count++ } END { print count }' ${pep_prefix}_nr_diamond.xls)
    count60=$(awk '$12 > 60 { count++ } END { print count }' ${pep_prefix}_nr_diamond.xls)
    echo -e "NR60\t$count60"
    echo -e "NR100\t$count100"
  # 获取基因列表
  grep -v ^# ${pep_prefix}_nr_diamond.xls \
    | awk 'BEGIN{FS="\t";OFS="\t"}{print $1"\t"$2}' \
    | sort | uniq \
    | sed '1i Gene_ID\tnr_ID' > ${pep_prefix}_nr_diamond_geneid.list
  # 基因description
  python $EXTRACT ${pep_prefix}_nr_diamond_geneid.list $nr_description ${pep_prefix}_nr_extract_description.tem # database使用表格处理效率太低，因此专门用这个脚本提取比对到的 nr 的基因的 description
  sed '1d ' ${pep_prefix}_nr_extract_description.tem | sort  |uniq | sed '1i Gene_ID\tnr_ID\tnr_Description' > ${pep_prefix}_nr_extract_description.tem.uniq # 去重
  csvtk -t join -f "Gene_ID;Gene_ID" ${pep_prefix}_nr_diamond_geneid.list ${pep_prefix}_nr_extract_description.tem.uniq --na - | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$4}' > ${pep_prefix}_nr_description.results # join
  rm ${pep_prefix}_nr_extract_description.tem ${pep_prefix}_nr_extract_description.tem.uniq
  python $MERGE_ROWS ${pep_prefix}_nr_description.results ${pep_prefix}_nr_description.merge.results
  mkdir NR
  mv ${pep_prefix}_nr* NR/.
}
##############################################################
# SwissProt                   ################################
##############################################################
function run_swiss_prot_map() {
  echo -e "\e[33m #################################\e[0m"
  echo -e "\e[33m ## run SwissProt map (Uniprot) ##\e[0m"
  echo -e "\e[33m #################################\e[0m"
  # 比对
  diamond blastp \
    --query $INPUT \
    --db $UNIPROT_DB \
    --max-target-seqs 1 \
    --max-hsps 1 \
    --evalue 1e-5 \
    --${SENSITIVE_MODE} \
    --outfmt 6 \
    --threads $THREADS \
    --quiet \
    --out ${pep_prefix}_uniprot_diamond.xls # 比对      --min-score 100 \
    # 统计一下score的数量分布
  count60=$(awk '$12 > 60 { count++ } END { print count }' ${pep_prefix}_uniprot_diamond.xls)
  count100=$(awk '$12 > 100 { count++ } END { print count }' ${pep_prefix}_uniprot_diamond.xls)
  echo -e "uniprot60 : $count60"
  echo -e "uniprot100 : $count100"
  # 提取基因信息
  grep -v ^# ${pep_prefix}_uniprot_diamond.xls \
    | awk 'BEGIN{FS="\t";OFS="\t"}{print $1"\t"$2}' \
    | sed 's/sp|//' \
    | sed 's/|/\t/' \
    | awk 'BEGIN{FS="\t";OFS="\t"}{print $1"\t"$2}' \
    | sort | uniq \
    | sed '1i Gene_ID\tUniprot_ID' > ${pep_prefix}_uniprot_diamond_geneid.list
  # 提取基因description,合并
  csvtk -t join -f "Uniprot_ID;Uniprot_ID" ${pep_prefix}_uniprot_diamond_geneid.list ${uniprto_description} | uniq > ${pep_prefix}_uniprot_description.results
  python $MERGE_ROWS ${pep_prefix}_uniprot_description.results ${pep_prefix}_uniprot_description.merge.results
  # GO
  python $EXTRACT ${pep_prefix}_uniprot_diamond_geneid.list $GO_idmapping_db ${pep_prefix}_uniprot_GOIDs.tmp #也可以用csvtk -t join
  python $GO_split ${pep_prefix}_uniprot_GOIDs.tmp ${pep_prefix}_uniprot_GOIDs.split.results
  python $MERGE_ROWS ${pep_prefix}_uniprot_GOIDs.split.results ${pep_prefix}_uniprot_GOIDs.merge.results
  #
  mkdir Swissprot
  rm ${pep_prefix}_uniprot_GOIDs.tmp
  mv ${pep_prefix}_uniprot_diamond.xls Swissprot/.
  mv ${pep_prefix}_uniprot_diamond_geneid.list Swissprot/.
  mv ${pep_prefix}_uniprot_description.results Swissprot/.
  mv ${pep_prefix}_uniprot_description.merge.results Swissprot/.
}
##############################################################
# Os                          ################################
##############################################################
function run_os_map() {
  echo -e "\e[33m #################################\e[0m"
  echo -e "\e[33m ## run Os ( IRGSP )            ##\e[0m"
  echo -e "\e[33m #################################\e[0m"
  # 比对
  diamond blastp \
    --query $INPUT \
    --db $OS_DB \
    --max-target-seqs 1 \
    --max-hsps 1 \
    --evalue 1e-5 \
    --${SENSITIVE_MODE} \
    --outfmt 6 \
    --threads $THREADS \
    --quiet \
    --out ${pep_prefix}_Os_diamond.xls
    # 提取基因信息
  grep -v ^# ${pep_prefix}_Os_diamond.xls \
    | awk 'BEGIN{FS="\t";OFS="\t"}{print $1"\t"$2}' \
    | sort | uniq \
    | sed '1i Gene_ID\tOs_ID' > ${pep_prefix}_Os_diamond_geneid.list
  # 提取基因description,合并
  csvtk -t join -f "Os_ID;Os_ID" ${pep_prefix}_Os_diamond_geneid.list $OS_description | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$3,$4,$5,$6}' | uniq > ${pep_prefix}_Os_description.results
  python $MERGE_ROWS ${pep_prefix}_Os_description.results ${pep_prefix}_Os_description.merge.results
  # 
  mkdir Os
  mv ${pep_prefix}_Os*  Os/.
}
##############################################################
# Ath                          ################################
##############################################################
function run_ath_map() {
  echo -e "\e[33m #################################\e[0m"
  echo -e "\e[33m ## run Ath ( IRGSP )           ##\e[0m"
  echo -e "\e[33m #################################\e[0m"
  # 比对
  diamond blastp \
    --query $INPUT \
    --db $ATH_DB \
    --max-target-seqs 1 \
    --max-hsps 1 \
    --evalue 1e-5 \
    --${SENSITIVE_MODE} \
    --outfmt 6 \
    --threads $THREADS \
    --quiet \
    --out ${pep_prefix}_Ath_diamond.xls
  # 提取基因信息
  grep -v ^# ${pep_prefix}_Ath_diamond.xls \
    | awk 'BEGIN{FS="\t";OFS="\t"}{print $1"\t"$2}' \
    | sort | uniq \
    | sed '1i Gene_ID\tAth_ID' > ${pep_prefix}_Ath_diamond_geneid.list
  # 提取基因description,合并
  csvtk -t join -f "Ath_ID;Ath_ID" ${pep_prefix}_Ath_diamond_geneid.list $ATH_description | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$3,$4,$8}' |uniq > ${pep_prefix}_Ath_description.results
  python $MERGE_ROWS ${pep_prefix}_Ath_description.results ${pep_prefix}_Ath_description.merge.results
  # 
  mkdir Ath
  mv ${pep_prefix}_Ath*  Ath/.
}
##############################################################
# TFs and PKs                   ##############################
##############################################################

function run_os_map() {
  echo -e "\e[33m #################################\e[0m"
  echo -e "\e[33m ## run TFs and PKs ( iTAK )        ##\e[0m"
  echo -e "\e[33m #################################\e[0m"
  # 比对
  perl $ITAK -p $THREADS -o ITAK  $INPUT
  cd ITAK
  ####################
  # TFs       ########
  ####################
  echo -e "\e[33m #################################\e[0m"
  echo -e "\e[33m ## get TFs                     ##\e[0m"
  echo -e "\e[33m #################################\e[0m"
  awk 'BEGIN{FS="\t";OFS="\t"} {split($4, arr, "->"); if(length(arr) > 1) print $1, arr[1], arr[2]; else print $1, $4, $4}' tf_classification.txt \
    | sort | uniq \
    | (echo -e "Gene_ID\tTF_Category\tTF_Name" && cat) >${pep_prefix}_TF.split.txt
  awk -F'\t' 'BEGIN {OFS="\t"} {print $1, "[ " $2 ": " $3 " ]"}' ${pep_prefix}_TF.split.txt \
    | sed '1d' \
    | sed '1i Gene_ID\tTF_Category_Name' >${pep_prefix}_TF.split.results
  python $MERGE_ROWS ${pep_prefix}_TF.split.results ${pep_prefix}_TF.merge.results
  ####################
  # PKs       ########
  ####################
  echo -e "\e[33m #################################\e[0m"
  echo -e "\e[33m ## get PKs                     ##\e[0m"
  echo -e "\e[33m #################################\e[0m"
  awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2}' shiu_classification.txt \
    | uniq |sort \
    | sed '1i Gene_ID\tPKs_Name' >${pep_prefix}_PK.split.txt
  awk -F'\t' 'BEGIN {OFS="\t"} {print $1, "[ " $2 " ]"}' ${pep_prefix}_PK.split.txt \
    | sed '1d' \
    | sed '1i Gene_ID\tPKs_Name' >${pep_prefix}_PK.split.results
  python $MERGE_ROWS ${pep_prefix}_PK.split.results ${pep_prefix}_PK.merge.results
  cd ..
}

##############################################################
# COG                         ################################
##############################################################
function run_COG_map() {
  echo "COG map"
  echo -e "\e[33m ################################\e[0m"
  echo -e "\e[33m ## run COG map                ##\e[0m"
  echo -e "\e[33m ################################\e[0m"
  # 比对
  diamond blastp \
    --query $INPUT \
    --db $COG_DB \
    --max-target-seqs 1 \
    --max-hsps 1 \
    --evalue 1e-5 \
    --${SENSITIVE_MODE} \
    --outfmt 6 \
    --threads $THREADS \
    --quiet \
    --out ${pep_prefix}_COG_diamond.xls
    # 统计一下score的数量分布
    count100=$(awk '$12 > 100 { count++ } END { print count }' ${pep_prefix}_COG_diamond.xls)
    count60=$(awk '$12 > 60 { count++ } END { print count }' ${pep_prefix}_COG_diamond.xls)
    echo -e "COG60 : $count60"
    echo -e "COG100 : $count100"
  # 提取基因信息
  grep -v ^# ${pep_prefix}_COG_diamond.xls \
    | awk 'BEGIN{FS="\t";OFS="\t"}{print $1"\t"$2}' \
    | sort | uniq \
    | sed '1i Gene_ID\tCOG_Protein' > ${pep_prefix}_COG_diamond_geneid.list
  # 合并
  csvtk -t join -f "COG_Protein;COG_Protein" ${pep_prefix}_COG_diamond_geneid.list ${COG_description} \
    | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$3,$4,$5,$6}' \
    | uniq > ${pep_prefix}_COG_description.txt
  awk -F'\t' 'BEGIN {OFS="\t"} {print $1, "[ " $2 ": " $3 "; " $4 " ]", $5}' ${pep_prefix}_COG_description.txt \
    | sed '1d' \
    | sed "1i Gene_ID\tCOG_ID_description\tCOG_Functional_category_ID" >${pep_prefix}_COG_description.split.results
  python $MERGE_ROWS ${pep_prefix}_COG_description.split.results ${pep_prefix}_COG_description.merge.results
  # 
  mkdir COG
  mv ${pep_prefix}_COG* COG/.
}
##############################################################
# KOG                         ################################
##############################################################
function run_KOG_map() {
  echo -e "\e[33m ################################\e[0m"
  echo -e "\e[33m ## run KOG map                ##\e[0m"
  echo -e "\e[33m ################################\e[0m"
  # 比对
  diamond blastp \
    --query $INPUT \
    --db $KOG_DB \
    --max-target-seqs 1 \
    --max-hsps 1 \
    --evalue 1e-5 \
    --${SENSITIVE_MODE} \
    --outfmt 6 \
    --threads $THREADS \
    --quiet \
    --out ${pep_prefix}_KOG_diamond.xls
    # 统计一下score的数量分布
    count100=$(awk '$12 > 100 { count++ } END { print count }' ${pep_prefix}_KOG_diamond.xls)
    count60=$(awk '$12 > 60 { count++ } END { print count }' ${pep_prefix}_KOG_diamond.xls)
    echo -e "KOG60 : $count60"
    echo -e "KOG100 : $count100"
  ## 获取基因列表
  grep -v ^# ${pep_prefix}_KOG_diamond.xls \
    | awk 'BEGIN{FS="\t";OFS="\t"}{print $1"\t"$2}' \
    | sort | uniq \
    | sed '1i Gene_ID\tKOG_Protein' > ${pep_prefix}_KOG_diamond_geneid.list
  # 提取基因description,合并
  csvtk -t join -f "KOG_Protein;KOG_Protein" ${pep_prefix}_KOG_diamond_geneid.list ${KOG_description} | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$4,$5,$3}' |uniq > ${pep_prefix}_KOG_description.txt
  awk -F'\t' 'BEGIN {OFS="\t"} {print $1, "[ " $2 ": " $3 " ]", $4}' ${pep_prefix}_KOG_description.txt \
    | sed '1d' \
    | sed "1i Gene_ID\tKOG_ID_description\tKOG_Functional_category_ID" >${pep_prefix}_KOG_description.split.results
  python $MERGE_ROWS ${pep_prefix}_KOG_description.split.results ${pep_prefix}_KOG_description.merge.results
  # 
  mkdir KOG
  mv ${pep_prefix}_KOG* KOG/.
}
##############################################################
# InterProscan                ################################
##############################################################
function run_interproscan_map() {
  echo -e "\e[33m #########################################\e[0m"
  echo -e "\e[33m ## run InterproScan map, include Pfam  ##\e[0m"
  echo -e "\e[33m #########################################\e[0m"
  # 比对
  bash $INTERPRO \
    -i $INPUT \
    --disable-precalc \
    --iprlookup \
    --cpu $THREADS \
    --goterms \
    --pathways \
    -f tsv \
    -o ${pep_prefix}_interproscan.xls >interproscan.log
  ########################################################################################
  # IPR                                                                                  #
  ########################################################################################
  mkdir interproscan
  echo -e "\e[33m## InterproScan : IPR     ##\e[0m"
  grep -v ^# ${pep_prefix}_interproscan.xls \
    | grep IPR \
    | awk 'BEGIN{FS="\t";OFS="\t"}{print $1"\t"$12"\t"$13}' \
    | sort | uniq \
    | sed 's/\"//g' \
    | sed '1i Gene_ID\tIPR_ID\tIPR_description' > ${pep_prefix}_IPR_interproscan.xls
  awk -F'\t' 'BEGIN {OFS="\t"} {print $1, "[ " $2 ": " $3 " ]"}' ${pep_prefix}_IPR_interproscan.xls \
    | sed '1d' \
    | sed '1i Gene_ID\tIPR_ID_description' >${pep_prefix}_IPR_interproscan.results
  python $MERGE_ROWS ${pep_prefix}_IPR_interproscan.results ${pep_prefix}_IPR_interproscan.merge.results
  grep -v ^# ${pep_prefix}_IPR_interproscan.xls \
    | sed '1d' \
    | cut -f1 \
    | sort | uniq \
    | sed '1i Gene_ID' >${pep_prefix}_IPR_interproscan_geneid.list
  mv interproscan.log interproscan/.
  mv ${pep_prefix}_IPR* interproscan/.
  ########################################################################################
  # 其他库                                                                               #
  ########################################################################################
  awk 'BEGIN{FS="\t";OFS="\t"} {print $1, $4, $5, $6}' ${pep_prefix}_interproscan.xls | sort | uniq > processed_interproscan.tsv
  # 获取唯一的 InterproScan 值
  cut -f2 processed_interproscan.tsv | tail -n +2 | sort | uniq | while read i
  do
    echo -e "\e[33m## InterproScan : $i      ##\e[0m"
    # 对于每个唯一值，处理保存的文件
    awk -v var="$i" -F'\t' 'BEGIN {OFS="\t"} { if ($2 == var) { print $1, $3, $4 } }' processed_interproscan.tsv \
      | sed 's/\"//g' \
      | sed "1i Gene_ID\t${i}_ID\t${i}_description" > ${pep_prefix}_${i}_interproscan.xls
    awk -F'\t' 'BEGIN {OFS="\t"} {print $1, "[ " $2 ": " $3 " ]"}' ${pep_prefix}_${i}_interproscan.xls \
      | sed '1d' \
      | sed "1i Gene_ID\t${i}_ID_description" > ${pep_prefix}_${i}_interproscan.results
    python $MERGE_ROWS ${pep_prefix}_${i}_interproscan.results ${pep_prefix}_${i}_interproscan.merge.results
    grep -v ^# ${pep_prefix}_${i}_interproscan.merge.results \
      | sed '1d' \
      | cut -f1 \
      | sort | uniq \
      | sed '1i Gene_ID' >${pep_prefix}_${i}_interproscan_geneid.list
    mv ${pep_prefix}_${i}_interproscan* interproscan/.
  done
  # 清理中间文件
  rm processed_interproscan.tsv
  ########################################################################################
  # InterProscan GO                                                                      #
  ########################################################################################
  grep "GO:" ${pep_prefix}_interproscan.xls | awk '!seen[$0]++' | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$4,$14}' | sed '1i Gene_ID\tDB\tGO_ID' | sed 's/|/;/g' > ${pep_prefix}_interproscan_GOIDs.tsv
  python $GO_split ${pep_prefix}_interproscan_GOIDs.tsv ${pep_prefix}_interproscan_GOIDs.split.tem
  sed '1d' ${pep_prefix}_interproscan_GOIDs.split.tem | sed 's/(/\t/' | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2}' | sort | uniq | sed '1i Gene_ID\tGO_ID' > ${pep_prefix}_interproscan_GOIDs.split.results
  python $MERGE_ROWS ${pep_prefix}_interproscan_GOIDs.split.results ${pep_prefix}_interproscan_GOIDs.merge.results
  rm ${pep_prefix}_interproscan_GOIDs.split.tem ${pep_prefix}_interproscan_GOIDs.tsv
  # 
  mv ${pep_prefix}_interproscan.xls interproscan/.
}
##############################################################
# EggNOG                      ################################
##############################################################
function run_eggnog_map() {
  echo -e "\e[33m ################################\e[0m"
  echo -e "\e[33m ## run EggNOG map             ##\e[0m"
  echo -e "\e[33m ################################\e[0m"
  # 比对
  emapper.py \
    -i $INPUT \
    --cpu $THREADS \
    -o ${pep_prefix}_EggNOG \
    --data_dir $EGGNOG_DB \
    --override \
    -m diamond \
    --dmnd_ignore_warnings \
    --dmnd_algo auto \
    --sensmode ${SENSITIVE_MODE} \
    --evalue 0.001 \
    --score 60 \
    --pident 40 \
    --query_cover 20 \
    --subject_cover 20 \
    --tax_scope auto \
    --target_orthologs all \
    --go_evidence non-electronic \
    --pfam_realign none \
    --report_orthologs \
    --decorate_gff yes \
    --excel >emapper.out 2> emapper.err
    
  # 基因列表
  grep -v ^# ${pep_prefix}_EggNOG.emapper.annotations \
    | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$8,$5}' \
    | sed '1i Gene_ID\tEggNOG_OGs\tEggNOG_Description' >${pep_prefix}_EggNOG.split.txt
  awk -F'\t' 'BEGIN {OFS="\t"} {print $1, "[ " $2 ": " $3 " ]" }' ${pep_prefix}_EggNOG.split.txt \
    | sed '1d' \
    | sed "1i Gene_ID\tEggNOG_Descripton_OGs" >${pep_prefix}_EggNOG.split.results
  python $MERGE_ROWS ${pep_prefix}_EggNOG.split.results ${pep_prefix}_EggNOG.merge.results
  grep -v ^# ${pep_prefix}_EggNOG.emapper.annotations \
    | awk 'BEGIN{FS="\t";OFS="\t"}{print $1}' \
    | sort |uniq \
    | sed '1i Gene_ID' >${pep_prefix}_EggNOG_geneid.list
  #GO
  grep "GO:" ${pep_prefix}_EggNOG.emapper.annotations | cut -f1,5,10 | sed 's/\,GO/\;GO/g' | sed '1i Gene_ID\tEGGnog_Ortho\tGO_ID' > ${pep_prefix}_EggNOG_GOIDs.tsv
  python $GO_split ${pep_prefix}_EggNOG_GOIDs.tsv ${pep_prefix}_EggNOG_GOIDs.split.results.tmp
  sed '1d' ${pep_prefix}_EggNOG_GOIDs.split.results.tmp | sort | uniq | sed '1i Gene_ID\tGO_ID' >${pep_prefix}_EggNOG_GOIDs.split.results
  python $MERGE_ROWS ${pep_prefix}_EggNOG_GOIDs.split.results ${pep_prefix}_EggNOG_GOIDs.merge.results
  #
  #rm ${pep_prefix}_EggNOG_GOIDs.tsv ${pep_prefix}_EggNOG_GOIDs.split.results.tmp
  mkdir eggnog
  mv ${pep_prefix}_EggNOG* eggnog/.
  mv emapper* eggnog/.
  mv eggnog/${pep_prefix}_EggNOG.emapper.annotations eggnog/${pep_prefix}_EggNOG.emapper.annotations.xls
  mv eggnog/${pep_prefix}_EggNOG_GOIDs.split.results .
  mv eggnog/${pep_prefix}_EggNOG_GOIDs.merge.results .
}
##############################################################
# GO                          ################################
##############################################################
function run_GO_online() {
  echo -e "\e[33m #################################################\e[0m"
  echo -e "\e[33m ## 合并 SwissProt、InterProscan、EggNOG GO注释 ##\e[0m"
  echo -e "\e[33m #################################################\e[0m"
  mkdir GO
  mv *results GO/.
  cd GO
  ls *split.results | while read if
  do
    sed -i '1d' $if
  done
  cat *split.results | sort | uniq | sed '1i Gene_ID\tGO_ID' > ${pep_prefix}_GOIDs.split.txt
  mv ${pep_prefix}_EggNOG_GOIDs.merge.results ${pep_prefix}_EggNOG_GOIDs.merge.results.txt
  mv ${pep_prefix}_uniprot_GOIDs.merge.results ${pep_prefix}_uniprot_GOIDs.merge.results.txt
  mv ${pep_prefix}_interproscan_GOIDs.merge.results ${pep_prefix}_interproscan_GOIDs.merge.results.txt
  mv ${pep_prefix}_EggNOG_GOIDs.split.results ${pep_prefix}_EggNOG_GOIDs.split.results.txt
  mv ${pep_prefix}_interproscan_GOIDs.split.results ${pep_prefix}_interproscan_GOIDs.split.results.txt
  mv ${pep_prefix}_uniprot_GOIDs.split.results ${pep_prefix}_uniprot_GOIDs.split.results.txt
  
  # description
  # all，注意，这一步有的GO是被舍弃的，所以在这一步的GO数量是少于原来输入的GO，这是因为GO数据库在更新，会废弃一些GO term
  csvtk -t join -f 'GO_ID;GO_ID' ${pep_prefix}_GOIDs.split.txt $GO_description > ${pep_prefix}_GOIDs.description.split.txt
  awk -F'\t' 'BEGIN {OFS="\t"} {print $1, "[ " $2 ": " $3 " ]", $4}' ${pep_prefix}_GOIDs.description.split.txt \
    | sed '1d' \
    |sed '1i Gene_ID\tGO_Term\tGO_Category' > ${pep_prefix}_GOIDs.description.split.results

  # 处理BP
  # awk -v prefix="$pep_prefix" 'NR==1{print > prefix "_GOIDs_BP.split.description.results"} NR > 1 && $4=="BP"{print > prefix "_GOIDs_BP.split.description.results"}' ${pep_prefix}_GOIDs.split.description_for_analisys.results1
  # 处理CC
  # awk -v prefix="$pep_prefix" 'NR==1{print > prefix "_GOIDs_CC.split.description.results"} NR > 1 && $4=="CC"{print > prefix "_GOIDs_CC.split.description.results"}' ${pep_prefix}_GOIDs.split.description_for_analisys.results1
  # 处理MF
  # awk -v prefix="$pep_prefix" 'NR==1{print > prefix "_GOIDs_MF.split.description.results"} NR > 1 && $4=="MF"{print > prefix "_GOIDs_MF.split.description.results"}' ${pep_prefix}_GOIDs.split.description_for_analisys.results1
  # 合并
  python $MERGE_ROWS ${pep_prefix}_GOIDs.description.split.results ${pep_prefix}_GOIDs.split.description.merge.results
  # 处理BP
  # python $MERGE_ROWS ${pep_prefix}_GOIDs_BP.split.description.results ${pep_prefix}_GOIDs_BP.split.description.merge.results
  # 处理CC
  # python $MERGE_ROWS ${pep_prefix}_GOIDs_CC.split.description.results ${pep_prefix}_GOIDs_CC.split.description.merge.results
  # 处理MF
  # python $MERGE_ROWS ${pep_prefix}_GOIDs_MF.split.description.results ${pep_prefix}_GOIDs_MF.split.description.merge.results
  # GO 基因list
  sed '1d' ${pep_prefix}_GOIDs.split.description.merge.results    | awk 'BEGIN{FS="\t";OFS="\t"}{print $1}' | sort | uniq | sed '1i Gene_ID' > ${pep_prefix}_GO_geneid.list
  # sed '1d' ${pep_prefix}_GOIDs_BP.split.description.merge.results | cut -f1 | sort | uniq | sed '1i Gene_ID' > ${pep_prefix}_GO_BP_geneid.list
  # sed '1d' ${pep_prefix}_GOIDs_CC.split.description.merge.results | cut -f1 | sort | uniq | sed '1i Gene_ID' > ${pep_prefix}_GO_CC_geneid.list
  # sed '1d' ${pep_prefix}_GOIDs_MF.split.description.merge.results | cut -f1 | sort | uniq | sed '1i Gene_ID' > ${pep_prefix}_GO_MF_geneid.list
  cd ..
}
##############################################################
# kofam                       ################################
##############################################################
function run_kofam_map() {
  echo "KEGG map"
  echo -e "\e[33m #########################\e[0m"
  echo -e "\e[33m ## run kofam           ##\e[0m"
  echo -e "\e[33m #########################\e[0m"
  kofamscan \
    -p $KO_PROFILE \
    -k $KO_LIST \
    --cpu $THREADS \
    --format mapper \
    -e 1e-10 \
    $INPUT \
    -o ${pep_prefix}_kofam.xls \
    --no-report-unannotated
    #####
  grep -v ^# ${pep_prefix}_kofam.xls \
    | sort |uniq \
    | sed '1i Gene_ID\tKO_ID' >${pep_prefix}_EggNOG_geneid.list
  csvtk -t -f "KO_ID;KO_ID" join ${pep_prefix}_EggNOG_geneid.list $KO_KEGG_DESCRIPTION >${pep_prefix}_KEGG.split.txt
  awk -F'\t' 'BEGIN {OFS="\t"} {
    print $1, 
          ($2 == "-" ? "[ - ]" : "[ " $2 ": " $3 " ]"),
          ($4 == "-" ? "[ - ]" : "[ " $4 ": " $5 " ]"),
          ($6 == "-" ? "[ - ]" : "[ " $6 ": " $7 " ]")
  }' ${pep_prefix}_KEGG.split.txt \
  | sed '1d' \
  | sed "1i Gene_ID\tKO_Gene_Description\tKEGG_PathwayID_Description\tKEGG_Pathway_Class" > ${pep_prefix}_EggNOG.split.results

  python $MERGE_ROWS ${pep_prefix}_EggNOG.split.results ${pep_prefix}_KEGG.merge.results
  # cut -f1,2 ${pep_prefix}_kofam.xls > ${pep_prefix}_kofam.only1.xls
  # python $KEGG_anno_from_kofam ${pep_prefix}_kofam.only1.xls $SPECIES
  # python $MERGE_ROWS ${pep_prefix}_kofam.only1.kegg_anno.xls ${pep_prefix}_kofam.merge.results
  mkdir kofam
  mv ${pep_prefix}_KEGG* kofam/.
  mv ${pep_prefix}_kofam* kofam/.
  mv tmp kofam/kofam_tmp
}

############################################################################################################
# run map
# 有不需要比对的，如eggnog，kegg，nr比较慢的，可以直接注释不跑
# eg: #echo "NR map"
# eg: #run_nr_map
###########################################################################################################################
if [[ -e ${INPUT} ]]; then
  echo -e "\e[34m                                      _  __  _                   _       _             \e[0m"
  echo -e "\e[34m                                     | |/ / (_)                 | |     (_)            \e[0m"
  echo -e "\e[34m                                     | ' /   _   _ __     __ _  | |__    _   _ __      \e[0m"
  echo -e "\e[34m                                     |  <   | | | '_ \   / _  | | '_ \  | | | '_ \     \e[0m"
  echo -e "\e[34m                                     | . \  | | | | | | | (_| | | |_) | | | | | | |    \e[0m"
  echo -e "\e[34m                                     |_|\_\ |_| |_| |_|  \__, | |_.__/  |_| |_| |_|    \e[0m"
  echo -e "\e[34m                                                          __/ |                        \e[0m"
  echo -e "\e[34m                                                         |___/                         \e[0m"
###########################################################################################################################
  # 获取geneid
  grep '^>' $INPUT | cut -d ' ' -f 1 | tr -d '>' | sort | uniq | sed '1i Gene_ID' > ${pep_prefix}_gene_id.list
  # 1.检查软件
  check_soft
  # 2.Swiss-Prot库比对
  run_swiss_prot_map
  # 3.KOG数据库比对
  run_KOG_map
  # 4. COG数据库比对
  run_COG_map
  # 5. InterProscan数据库比对
  run_interproscan_map
  # 6. EggNOG数据库比对
  run_eggnog_map
  # 7.根据结果获取GO注释
  run_GO_online
  # 8.kofam比对
  run_kofam_map
  # 9. NR数据库比对
  run_nr_map
  ############################################################################################################
  #合并
  csvtk -t join -f1 \
    ${pep_prefix}_gene_id.list \
    Swissprot/${pep_prefix}_uniprot_description.merge.results \
    Os/${pep_prefix}_Os_description.merge.results \
    Ath/${pep_prefix}_Ath_description.merge.results \
    NR/${pep_prefix}_nr_description.merge.results \
    ITAK/${pep_prefix}_TF.merge.results \
    ITAK/${pep_prefix}_PK.merge.results \
    kofam/${pep_prefix}_KEGG.merge.results \
    GO/${pep_prefix}_GOIDs.split.description.merge.results \
    COG/${pep_prefix}_COG_description.merge.results \
    KOG/${pep_prefix}_KOG_description.merge.results \
    eggnog/${pep_prefix}_EggNOG.merge.results \
    interproscan/${pep_prefix}_Pfam_interproscan.merge.results \
    interproscan/${pep_prefix}_IPR_interproscan.merge.results \
    interproscan/${pep_prefix}_CDD_interproscan.merge.results \
    interproscan/${pep_prefix}_PANTHER_interproscan.merge.results \
    interproscan/${pep_prefix}_Coils_interproscan.merge.results \
    interproscan/${pep_prefix}_FunFam_interproscan.merge.results \
    interproscan/${pep_prefix}_PRINTS_interproscan.merge.results \
    interproscan/${pep_prefix}_Gene3D_interproscan.merge.results \
    interproscan/${pep_prefix}_ProSitePatterns_interproscan.merge.results \
    interproscan/${pep_prefix}_Hamap_interproscan.merge.results \
    interproscan/${pep_prefix}_ProSiteProfiles_interproscan.merge.results \
    interproscan/${pep_prefix}_SMART_interproscan.merge.results \
    interproscan/${pep_prefix}_MobiDBLite_interproscan.merge.results \
    interproscan/${pep_prefix}_SUPERFAMILY_interproscan.merge.results \
    interproscan/${pep_prefix}_NCBIfam_interproscan.merge.results \
    --left-join --na - >${pep_prefix}_total.annatation.xls
###################################################################################################################
###################################################################################################################
  #################
  ##  统计       ##
  #################
  #python $CONT_GENE -i ${pep_prefix}_gene_id.list -o ${pep_prefix}_counts.xls
  ####################
  ## 统计各个数据库 ##
  ####################
  echo -e "\e[33m 统计数据库注释数量 \e[0m"
  find . -type f -name "*_merge.results" | while read -r file
  do
    lines=$(($(wc -l < "$file") - 1))                                                                        # 计算行数并减去1
    filename=$(basename "$file")                                                                             # 获取文件名
    clean_name=${filename#${pep_prefix}_}                                                                    # 去除文件名中的 ${pep_prefix}_
    echo -e "${clean_name}\t${lines}"
  done > ${pep_prefix}_Gene_Count.xls
###################################################################################################################
###################################################################################################################
  mkdir temp_first_col_files                                                                                 # 创建一个临时目录用于存储处理后的第一列文件
  echo -e "Filename\tCount" > ${pep_prefix}_Gene_Count.xls                                                   # 用于保存统计量的文件
  # 查找所有以 *merge.results 结尾的文件并处理
  find . -type f -name "*merge.results" | while read -r file; do
    filename=$(basename "$file") 
    clean_name=${filename#${pep_prefix}_}
    clean_name=${clean_name%.merge.results}                                                                  # 获取简化后的文件名，去除前缀和 .merge.results 后缀
    col_count=$(tail -n +2 "$file" | cut -f1 | sort | uniq | wc -l)                                          # 获取第一列的唯一元素个数（排除标题行）
    echo -e "${clean_name}\t${col_count}" >> ${pep_prefix}_Gene_Count.xls                                    # 输出文件名和列数到统计文件
    echo "$clean_name" > "temp_first_col_files/${clean_name}.txt"                                            # 提取第一列的唯一元素（排除标题行）并保存到临时文件，添加处理后的文件名作为列名
    tail -n +2 "$file" | cut -f1 | sort | uniq >> "temp_first_col_files/${clean_name}.txt"
  done
  paste temp_first_col_files/*.txt > ${pep_prefix}_All_annotation_combine.xls                                # 使用 paste 命令合并所有提取的第一列到一个文件中
  rm -r temp_first_col_files                                                                                 # 清理临时文件和目录
  #####
  cat ${pep_prefix}_All_annotation_combine.xls | sed '1d' | tr '\t' '\n' | grep -v '^$' | sort | uniq | sed '1i Gene_ID' > ${pep_prefix}_All_merge_Gene_list.xls
  line_count=$(($(wc -l < ${pep_prefix}_All_merge_Gene_list.xls) - 1))                                       # 计算文件行数并减去一
  echo -e "Total_annotation\t${line_count}" >> ${pep_prefix}_Gene_Count.xls                                  # 将有注释到的基因数量结果追加到 Gene_Count 文件
  count=$(grep -c '^>' "$INPUT")                                                                             # 统计基因数量
  echo -e "Total_Genes\t${count}" >> ${pep_prefix}_Gene_Count.xls                                            # 将所有的基因数量追加到 Gene_Count 文件
  #########
  ####统计
  #########
  awk -v count="$count" 'BEGIN {FS="\t"; OFS="\t"; print "Database", "Counts", "Percent"} 
    NR>1 {print $1, $2, ($2/count)*1}' ${pep_prefix}_Gene_Count.xls > ${pep_prefix}_Gene_Count_Updated.xls   # 读取 ${pep_prefix}_Gene_Count.xls 文件并添加新的计算列
  rm ${pep_prefix}_Gene_Count.xls
###################################################################################################################
###################################################################################################################
else
  echo "请检查输入文件"
  exit 1
fi

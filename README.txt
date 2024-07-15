数据库包括：nr、uniprot、GO、COG、KOG、KEGG、OS（水稻）、ATH（拟南芥）、PLANT_TFDB、STRINGDB、iTAK
具体等后续更新，文件较大如有急需分享，可私聊QQ：719837929

更新历史：
# Function_Annotition_And_Get_Description-v0.5.1.sh  注释功能函数化
# Function_Annotition_And_Get_Description-v0.5.2.sh  将KEGG与EGGNOG注释的结果直接实现在线爬虫，并且舍去uniprot的结果对应关系。kofam比对若出现一个基因对应多个KO，只保留第一个KO。
# Function_Annotition_And_Get_Description-v0.6.1     使用累计一致性确定确定同源基因；使用swissprot+interproscan+eggnog并集，构建GO注释。如果KGEE有多个Kgene，则保留。
# Function_Annotition_And_Get_Description-v0.6.2     添加转录因子（TF）数据库以及激酶数（PK）据库

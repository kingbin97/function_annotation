#!/usr/bin/env python
import sys

def before_GOstat(f1, f2):
    #f2.write("Gene_ID\tGO_ID\n")  # 写入输出文件头部
    
    for i in f1.readlines():
        j = i.rstrip().split('\t')
        gene = j[0]  # 基因
        kos  = j[2].split(';')  # KO编号列表
        
        for ko in kos:
            m = f"{gene}\t{ko}\n"  # 按制表符分割基因和KO编号，并添加换行符
            f2.write(m)  # 写入输出文件

# 从命令行参数获取输入和输出文件路径
input_file = sys.argv[1]
output_file = sys.argv[2]

# 打开输入和输出文件
with open(input_file, 'r') as f1, open(output_file, 'w') as f2:
    before_GOstat(f1, f2)
    
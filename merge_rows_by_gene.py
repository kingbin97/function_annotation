#!/usr/bin/env python

import sys
import csv

csv.field_size_limit(sys.maxsize)  # 增加字段大小限制

def merge_info(input_file, output_file, delimiter=';'):
    data = {}
    with open(input_file, 'r', encoding='utf-8', errors='ignore') as file:
        reader = csv.reader(file, delimiter='\t')
        header = next(reader)  # 读取标题行
        id_column_name = header[0]  # 第一列的标题
        other_columns = range(1, len(header))  # 除第一列外的其它所有列
        
        for row in reader:
            key_value = row[0]  # 第一列的值作为键
            
            if key_value not in data:
                data[key_value] = {column: set() for column in other_columns}
            
            for i in other_columns:
                data[key_value][i].add(row[i])  # 使用add代替append添加到集合中

    with open(output_file, 'w') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerow(header)
        
        for key_value, info in data.items():
            row = [key_value]
            for i in other_columns:
                # 使用join将集合转换回字符串，集合自动去除重复项
                row.append(delimiter.join(sorted(info[i])))  # 如果需要按顺序排序可以使用sorted()
            
            writer.writerow(row)

if __name__ == '__main__':
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        print("Usage:")
        print("  python script.py <input_file.tsv> <output_file.tsv> [delimiter]")
        print("\nDescription:")
        print("  Merges rows in a TSV file based on the first column.")
        print("  Concatenates values of other columns, separating them with the provided delimiter (default is ';').")
        print("\nArguments:")
        print("  <input_file.tsv>  - Path to the input TSV file.")
        print("  <output_file.tsv> - Path where the output TSV file will be saved.")
        print("  [delimiter]       - Optional delimiter for joining values (default ';').")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    delimiter = sys.argv[3] if len(sys.argv) == 4 else ';'
    merge_info(input_file, output_file, delimiter)

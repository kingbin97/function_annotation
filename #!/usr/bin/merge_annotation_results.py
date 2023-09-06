#!/usr/bin/env python
import os
import pandas as pd

files = os.listdir()
results_list = []
merged_df = pd.DataFrame()
#print(merged_df)
for results in files:
    # print(results)
    if results.strip().endswith("merge.results"):
        #print(results)
        results_list.append(results.strip())
    elif results.strip().endswith("GeneID.list.txt"):
        df_id = pd.read_csv(results.strip(), delimiter='\t')
        merged_df = df_id
        #print(df_id)
# print(results_list)
for results_id in results_list:
    results_df = pd.read_csv(results_id.strip(), delimiter='\t')
    #print(results_df)
    merged_df = pd.merge(merged_df, results_df, on='Gene_ID', how='left')
    merged_df.fillna('NA', inplace=True)
# print(merged_df)
merged_df.to_csv('total_annotation_results.tsv', sep='\t', index=False)
merged_df.to_excel('total_annotation_results.xlsx', index=False)

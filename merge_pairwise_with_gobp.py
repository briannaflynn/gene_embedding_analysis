#!/usr/bin/env python3
#Date: Dec 08, 2024
#Author: Brianna Flynn
#Description:
#This script is for merging the correlation or cosine sim pairwise data with GOBP pairwise data, this format is compatible with the scikit learn metric calculation class (see calculate_benchmark_metrics.py)

import pandas as pd
import sys
import time

def merge_pairwise_dataframes(df1, df2, col_a='gene1', col_b='gene2', suffixes=('_left', '_right')):
    """
    Merges two dataframes on pairwise columns, treating (A, B) as equivalent to (B, A).
    
    Parameters:
    - df1: First dataframe
    - df2: Second dataframe
    - col_a: Name of the first column in the pair, defaults to gene1
    - col_b: Name of the second column in the pair, defaults to gene2
    - suffixes: Suffixes to append to overlapping column names in the resulting dataframe
    
    Returns:
    - Merged dataframe
    """
    """
    # Normalize the pairwise columns in both dataframes
    df1['normalized_a'] = df1[[col_a, col_b]].min(axis=1)
    print(df1)
    df1['normalized_b'] = df1[[col_a, col_b]].max(axis=1)
    print(df1)
    
    df2['normalized_a'] = df2[[col_a, col_b]].min(axis=1)
    print(df2)
    df2['normalized_b'] = df2[[col_a, col_b]].max(axis=1)
    print(df2)

    print('\n\n') 
    print(df1)
    print(df2)
    
    # Merge on the normalized columns
    merged_df = pd.merge(
        df1, 
        df2, 
        on=['normalized_a', 'normalized_b'], 
        suffixes=suffixes
    )
    
    # Drop the helper columns used for normalization
    merged_df = merged_df.drop(columns=['normalized_a', 'normalized_b'])
    """

    merged_df = pd.merge(df1, df2, on=[col_a, col_b], suffixes=suffixes)
    return merged_df

def df_rename_by_idx(df):
    df.columns.values[0] = 'gene1'
    df.columns.values[1] = 'gene2'

def filter_by_obs(df, filter_df):
    
    print(df)
    filter_df_gene1 = filter_df.rename(columns={'Domain':'gene1'})
    data = filter_df_gene1.merge(df, on = 'gene1')
    filter_df_gene2 = filter_df.rename(columns={'Domain':'gene2'})
    
    result = data.merge(filter_df_gene2, on = 'gene2')
    print(result)
    print(df.shape[0] - result.shape[0])
    return result

# start = time.time()
# print('reading data, could take a while')
# df_pairwise = pd.read_pickle(sys.argv[1])
# df_rename_by_idx(df_pairwise)
# filter_obs_df = pd.read_csv(sys.argv[2])
# df_pair = filter_by_obs(df_pairwise, filter_obs_df)

# df_gobp = pd.read_pickle(sys.argv[3])
# df_rename_by_idx(df_gobp)

# print('startng merge')

# merged = merge_pairwise_dataframes(df_pair, df_gobp) 
# print(merged)
# end_main_merge = time.time()
# # merged.to_pickle('pairwise_150_gobp_100_lvl2_debug.pkl')
# end_export = time.time()

# print('from start to merge', end_main_merge - start)
# print('from start to export', end_export - start)
#merged.to_csv('pairwise_150_gobp_100_lvl2.csv')

#!/usr/bin/env python
# coding: utf-8
"""
dual_guide_aggregate.py

This script combines sgRNA and clonal barcode information for each mouse sample
from a collection of output files generated in previous steps. It searches for 
barcode and cluster files within an input directory structure, merges these with 
reference sgRNA data, deduplicates the merged data by grouping, and writes the 
results to CSV files. 
- Searches recursively for files in the `Clonal_barcode` directories,
- Merges barcode, cluster, and bartender output files,
- Combines this merged data with reference sgRNA information,
- Deduplicates the data by grouping on key columns, and
- Outputs per-sample CSV files as well as a final combined CSV file for downstream analysis.

Usage:
    python dual_guide_aggregate.py --a <input_folder> --o <output_prefix>
"""

import pandas as pd
import argparse
import glob
import os

def merge_barcode_and_sgRNA_output(barcode_file, cluster_file, bartender_input_file):
    """
    Merge data from barcode, cluster, and bartender files.

    Parameters:
      barcode_file (str): Path to the barcode CSV file.
      cluster_file (str): Path to the cluster CSV file.
      bartender_input_file (str): Path to the bartender file.

    Returns:
      pd.DataFrame: Merged DataFrame.
    """
    # Load barcode data (dropping the 'Frequency' column)
    barcode_df = pd.read_csv(barcode_file).drop(columns=['Frequency'])
    # Load cluster data (dropping unnecessary columns)
    cluster_df = pd.read_csv(cluster_file).drop(columns=['Cluster.Score', 'time_point_1'])
    # Load bartender data; no header and two columns: Clonal_barcode and Read_ID
    bartender_df = pd.read_csv(bartender_input_file, sep=',', names=['Clonal_barcode', 'Read_ID'], header=None)
    
    # Merge barcode and cluster data on 'Cluster.ID'
    merged_df = pd.merge(barcode_df, cluster_df, how='inner', on='Cluster.ID')
    merged_df.rename(columns={'Unique.reads': 'Clonal_barcode', 'Center': 'Clonal_barcode_center'}, inplace=True)
    # Drop the 'Cluster.ID' and merge with bartender data on 'Clonal_barcode'
    merged_df = merged_df.drop(columns=['Cluster.ID']).merge(bartender_df, on='Clonal_barcode', how='right')
    
    return merged_df

def combine_sgRNA_barcode_from_same_mouse(input_folder):
    """
    Combine sgRNA and clonal barcode information from files in a given folder.

    Parameters:
      input_folder (str): The path to the folder containing sample data.

    Returns:
      tuple: (deduplicated_raw, deduplicated_complete) DataFrames.
    """
    # Find all cluster files within Clonal_barcode subfolders (recursive)
    pattern = '**/Clonal_barcode/*_cluster.csv'
    barcode_files = glob.glob(os.path.join(input_folder, pattern), recursive=True)
    # The reference sgRNA data is assumed to be in 'Intermediate_df.csv' at the input folder level.
    ref_file = os.path.join(input_folder, 'Intermediate_df.csv')
    barcode_dataframes = []
    
    for cluster_file in barcode_files:
        # Construct paths for the barcode and bartender files
        barcode_file = cluster_file.replace('cluster', 'barcode')
        bartender_input_file = cluster_file.replace('_cluster.csv', '.bartender')
        
        # Merge data from the three files
        merged_df = merge_barcode_and_sgRNA_output(barcode_file, cluster_file, bartender_input_file)
        barcode_dataframes.append(merged_df)
    
    # Combine data from all files into one DataFrame
    combined_df = pd.concat(barcode_dataframes).reset_index(drop=True)
    sgRNA_ref_df = pd.read_csv(ref_file)
    final_df = combined_df.merge(sgRNA_ref_df, on=['Read_ID', 'Clonal_barcode'])
    
    # Deduplicate data with two grouping strategies
    group_cols_raw = ['gRNA_combination', 'Clonal_barcode_center', 'gRNA1', 'gRNA2', 'Clonal_barcode', 'Sample_ID']
    deduplicated_raw = final_df.groupby(group_cols_raw, as_index=False)['Read_ID'].count()
    deduplicated_raw.rename(columns={'Read_ID': 'Frequency'}, inplace=True)
    
    group_cols_complete = ['gRNA_combination', 'Clonal_barcode_center', 'gRNA1', 'gRNA2', 'Sample_ID']
    deduplicated_complete = final_df.groupby(group_cols_complete, as_index=False)['Read_ID'].count()
    deduplicated_complete.rename(columns={'Read_ID': 'Frequency', 'Clonal_barcode_center': 'Clonal_barcode'}, inplace=True)
    
    return deduplicated_raw, deduplicated_complete

def main():
    parser = argparse.ArgumentParser(description='Combine sgRNA and clonal barcode information')
    parser.add_argument("--a", required=True, help="Input bartender folder containing sample subdirectories")
    parser.add_argument("--o", required=True, help="Prefix for output files")
    args = parser.parse_args()

    output_prefix = args.o
    input_folder = args.a
    # Get a list of subfolders (each corresponding to a sample/mouse)
    subfolders = [os.path.join(input_folder, name) for name in os.listdir(input_folder) 
                  if os.path.isdir(os.path.join(input_folder, name))]
    final_dataframes = []

    for subfolder in subfolders:
        sample_id = os.path.basename(subfolder)
        output_file1 = f"{output_prefix}{sample_id}/Combined_ND_df.csv"
        output_file2 = f"{output_prefix}{sample_id}/Combined_deduplexed_df.csv"
        df1, df2 = combine_sgRNA_barcode_from_same_mouse(subfolder)
        # Ensure output directories exist
        os.makedirs(os.path.dirname(output_file1), exist_ok=True)
        os.makedirs(os.path.dirname(output_file2), exist_ok=True)
        df1.to_csv(output_file1, index=False)
        df2.to_csv(output_file2, index=False)
        
        final_dataframes.append(df2)
    
    # Combine deduplicated data from all samples into one file
    final_combined_df = pd.concat(final_dataframes)
    final_combined_df.to_csv(f"{output_prefix}gRNA_clonalbarcode_combined.csv", index=False)

if __name__ == "__main__":
    main()

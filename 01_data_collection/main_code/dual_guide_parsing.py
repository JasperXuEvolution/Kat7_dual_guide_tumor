#!/usr/bin/env python
# coding: utf-8
"""
dual_guide_parsing.py

This script extracts gRNA and clonal barcode information from paired-end
FASTQ files (gzipped) using regular expressions and reference sgRNA data.
It outputs summary statistics and writes separate files for "Expected" and
"Unexpected" reads, as well as bartender input file for downstream processing.

Usage:
    python dual_guide_parsing.py --a1 <R1_fastq.gz> --a2 <R2_fastq.gz> --b <guide_reference.csv> --o <output_directory>

Input:
- The input FASTQ files must be gzipped.
- A reference CSV file containing sgRNA data with at least these two columns:
  - `Position` (with values such as `G1` and `G2`)
  - `gRNA_complete` (the sgRNA sequence that designed, it may have first 'G' excluded)    
    
"""

import gzip
import regex
import argparse
import pandas as pd

def find_reverse_complementary(input_string):
    """
    Returns the reverse complementary sequence of the given DNA string.
    Supports uppercase and lowercase letters.
    """
    temp_dic = {
        'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G', 'N': 'N',
        'a': 't', 'g': 'c', 't': 'a', 'c': 'g', 'n': 'n'
    }
    # Reverse the string and convert each base using the dictionary
    return ''.join(temp_dic.get(x, x) for x in input_string[::-1])

def main():
    parser = argparse.ArgumentParser(
        description='Extract gRNA and clonal barcode information from paired-end FASTQ gz files.'
    )
    parser.add_argument("--a1", required=True, help="Input FASTQ gz file R1")
    parser.add_argument("--a2", required=True, help="Input FASTQ gz file R2")
    parser.add_argument("--b", required=True, help="Input CSV file of reference sgRNA")
    parser.add_argument("--o", required=True, help="Output directory")
    args = parser.parse_args()
    
    fastqgz_input_address1 = args.a1
    fastqgz_input_address2 = args.a2
    ref_address = args.b
    output_dir = args.o
    
    # Load reference sgRNA data
    ref_sgRNA_df = pd.read_csv(ref_address)
    gRNA1_list = ref_sgRNA_df[ref_sgRNA_df.Position == 'G1']['gRNA_complete'].to_list()
    gRNA2_list = ref_sgRNA_df[ref_sgRNA_df.Position == 'G2']['gRNA_complete'].to_list()
    
    # Compile regex patterns for sequence extraction
    # Pattern for read1: Extract a 16 bp barcode and a sgRNA (16-21 bp) between fixed sequence markers.
    temp_pattern1 = regex.compile('TAGTT(.{16})TATGG(.{16,21})GTTTA')
    # Pattern for read2: Extract a sgRNA (16-21 bp) from the reverse complemented read.
    temp_pattern2 = regex.compile('TGTTG(.{16,21})GTTTG')
    
    # Initialize counters and lists for collecting output data.
    total_reads = 0
    extracted_reads = 0
    matched_reads = 0
    sample_id = output_dir.split('/')[-1]  # Use the last part of output directory as sample ID
    
    gRNA1_list_out, gRNA2_list_out = [], []
    read_id_list, clonal_barcode_list, label_list = [], [], []
    
    with gzip.open(fastqgz_input_address1, 'rt') as handler1, gzip.open(fastqgz_input_address2, 'rt') as handler2:
        # Read the first record from both files.
        read_id = handler1.readline().rstrip()  # R1 ID
        handler2.readline().rstrip()            # R2 ID (not used)
        seq1 = handler1.readline().rstrip()     # R1 sequence
        seq2 = find_reverse_complementary(handler2.readline().rstrip())  # R2 sequence (reverse complemented)
        handler1.readline()  # skip '+'
        handler1.readline()  # skip quality
        handler2.readline()  # skip '+'
        handler2.readline()  # skip quality
        
        while read_id:
            total_reads += 1
            
            # Apply regex patterns to extract data.
            match1 = temp_pattern1.search(seq1)
            match2 = temp_pattern2.search(seq2)
            if match1 and match2:
                extracted_reads += 1
                clonal_barcode = match1.group(1)
                gRNA1 = match1.group(2)
                gRNA2 = match2.group(1)
                gRNA1_list_out.append(gRNA1)
                gRNA2_list_out.append(gRNA2)
                read_id_list.append(read_id)
                clonal_barcode_list.append(clonal_barcode)
                
                # Classify the read as 'Expected' if both sgRNAs match the reference.
                if (gRNA1 in gRNA1_list) and (gRNA2 in gRNA2_list):
                    matched_reads += 1
                    label_list.append('Expected')
                else:
                    label_list.append('Unexpected')
            
            # Read next record from both FASTQ files.
            read_id = handler1.readline().rstrip()  # R1 ID
            handler2.readline().rstrip()            # R2 ID (not used)
            seq1 = handler1.readline().rstrip()       # R1 sequence
            seq2 = find_reverse_complementary(handler2.readline().rstrip())  # R2 sequence (reverse complemented)
            handler1.readline()  # skip '+'
            handler1.readline()  # skip quality
            handler2.readline()  # skip '+'
            handler2.readline()  # skip quality
            
    # Print a summary of the processing.
    print(f"Sample {sample_id} has a total of {total_reads} reads. "
          f"{extracted_reads} reads ({extracted_reads/total_reads:.3f}) have barcode and sgRNA. "
          f"{matched_reads} reads ({matched_reads/total_reads:.3f}) have expected sgRNA.")
    
    # Create a DataFrame with the extraction results.
    Final_df = pd.DataFrame({
        'gRNA1': gRNA1_list_out,
        'gRNA2': gRNA2_list_out,
        'Clonal_barcode': clonal_barcode_list,
        'Read_ID': read_id_list,
        'Sample_ID': sample_id,
        'Class': label_list
    })
    Final_df['gRNA_combination'] = Final_df['gRNA1'] + '_' + Final_df['gRNA2']
    
    # Write out Unexpected reads.
    unexpected_file = f"{output_dir}/Unexpected_reads.csv"
    Final_df[Final_df.Class == 'Unexpected'].to_csv(unexpected_file, index=False)
    
    # Write out Expected reads to an intermediate file.
    intermediate_file = f"{output_dir}/Intermediate_df.csv"
    Final_df[Final_df.Class != 'Unexpected'].to_csv(intermediate_file, index=False)
    
    # Group expected reads by gRNA combination and create files for downstream processing.
    expected_df = Final_df[Final_df.Class != 'Unexpected']
    sgRNA_groups = expected_df.groupby("gRNA_combination")
    bartender_input_file = f"{output_dir}/Bartender_input_address"
    with open(bartender_input_file, 'w') as file_a:
        for group_name, group_df in sgRNA_groups:
            out_filename = f"{output_dir}/Clonal_barcode/{group_name}.bartender"
            file_a.write(out_filename + '\n')
            group_df[['Clonal_barcode', 'Read_ID']].to_csv(out_filename, sep=',', header=False, index=False)

if __name__ == "__main__":
    main()

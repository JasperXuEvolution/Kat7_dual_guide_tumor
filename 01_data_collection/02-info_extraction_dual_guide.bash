#!/bin/bash
# 02-info_extraction_dual_guide.bash

#SBATCH --job-name=UltraSeq_dual_guide_pipeline
#SBATCH --mail-user=xhq@stanford.edu
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=50g 
#SBATCH --time=2-00:00:00
#SBATCH --account=mwinslow
#SBATCH --partition=batch
#SBATCH --output=./log/%x_%j.out              # Standard output log file
#SBATCH --error=./log/%x_%j.err               # Standard error log file

# Source the configuration file 
source ../config.sh

# Load required modules and activate the Conda environment
module load adapterremoval/2.3.1
source ~/miniconda3/etc/profile.d/conda.sh 
conda activate UltraSeq

# --- Directories and Input Files ---
working_dir="$PROJECT_DIR/01_data_collection"
input_data_info_address="$working_dir/data/NGS_address"
guide_ref="$working_dir/data/guide_reference-Kat7_dual_guide.csv"

# Explanation for input files:
#   input_data_info_address:
#     This file should be a text file where each line is comma-separated with three fields:
#       1. The path to the FASTQ file for read 1.
#       2. The path to the FASTQ file for read 2.
#       3. The sample ID (used to name subdirectories for output).
#
#   guide_ref:
#     This CSV file is used as a reference for guide sequences. It is required for the dual guide parsing step.

# Define base directories for outputs:
step1_address="$working_dir/data/Trimmed_reads"     # For trimmed FASTQ files (Step 1)
step2_address="$working_dir/data/Bartender"           # Output from dual_guide_parsing.py (Step 2)
step4_address="$working_dir/data/Processed_data"      # Final processed/aggregated data (Step 4)

# Create base output directories if they don't exist
mkdir -p "$step1_address"
mkdir -p "$step2_address"
mkdir -p "$step4_address"

# --- Step 1: Sequence Trimming using AdapterRemoval ---
# For each sample, trim adapter sequences and store trimmed FASTQ files in step1_address/$sample_id
while read -r line; do
    # Extract original FASTQ paths and sample ID (comma-separated values)
    r1=$(echo "$line" | cut -d',' -f1)
    r2=$(echo "$line" | cut -d',' -f2)
    sample_id=$(echo "$line" | cut -d',' -f3)
    
    # Create a subdirectory for the sample in the trimmed reads folder
    sample_trim_dir="$step1_address/$sample_id"
    mkdir -p "$sample_trim_dir"
    
    # Define final output file names for trimmed reads (gzipped)
    trimmed_r1="$sample_trim_dir/${sample_id}_R1.trimmed.fastq.gz"
    trimmed_r2="$sample_trim_dir/${sample_id}_R2.trimmed.fastq.gz"

    echo "Trimming adapters for sample $sample_id..."
    
    # Run AdapterRemoval with --basename so that *all* output goes into $sample_trim_dir
    # The final paired reads will be named something like:
    #   ${sample_id}.pair1.truncated.gz
    #   ${sample_id}.pair2.truncated.gz
    # We then rename them to match your desired naming convention.
    AdapterRemoval \
        --file1 "$r1" \
        --file2 "$r2" \
        --basename "$sample_trim_dir/${sample_id}" \
        --adapter1 "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG" \
        --adapter2 "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT" \
        --trimns \
        --trimqualities \
        --gzip

    # Rename the final paired read outputs to your preferred names
    mv "$sample_trim_dir/${sample_id}.pair1.truncated.gz" "$trimmed_r1"
    mv "$sample_trim_dir/${sample_id}.pair2.truncated.gz" "$trimmed_r2"

   # --- Step 2: Generate Bartender Input (Dual Guide Parsing) ---
   # Create a sample-specific folder in the Bartender directory
   sample_bartender_dir="$step2_address/$sample_id"
   mkdir -p "$sample_bartender_dir"
   mkdir -p "$sample_bartender_dir/Clonal_barcode"

   # Run dual_guide_parsing.py using the trimmed files
   python3 "$working_dir/main_code/dual_guide_parsing.py" --a1 "$trimmed_r1" --a2 "$trimmed_r2" \
       --b "$guide_ref" --o "$sample_bartender_dir"

   # --- Step 3: Run Bartender Clustering ---
   # For each file listed in the Bartender_input_address file, run the clustering command
   while read -r line2; do 
      new_name=${line2/.bartender/}
      bartender_single_com -z -1 -d 1 -l 5 -f "$line2" -o "$new_name"
   done < "$sample_bartender_dir/Bartender_input_address"

   # --- Step 4: Create Processed Data Folder for the Sample ---
   sample_processed_dir="$step4_address/$sample_id"
   mkdir -p "$sample_processed_dir"

done < "$input_data_info_address"

# --- Final Step: Aggregate Data Across Samples ---
python3 "$working_dir/main_code/dual_guide_aggregate.py" --a "$step2_address" --o "$step4_address/"

# Report SLURM job statistics
sacct --format=JobID,JobName,Submit,Start,End,State,Partition,ReqTRES%30,CPUTime,MaxRSS,NodeList%30 --units=M -j $SLURM_JOBID

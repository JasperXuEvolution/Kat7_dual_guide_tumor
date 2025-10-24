#!/bin/bash
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=01-Bootstrapping-Kat7_DualGuide-KT
#SBATCH --mail-user=xhq@stanford.edu
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=100g 
#SBATCH --time=10-24:00:00
#SBATCH --account=mwinslow
#SBATCH --partition=batch
#SBATCH --output=./log/%x_%j.out              # Standard output log file
#SBATCH --error=./log/%x_%j.err               # Standard error log file

# general input and output address
source ../config.sh
source ~/miniconda3/etc/profile.d/conda.sh 
conda activate UltraSeq
input_tumor_data_address="${PROJECT_DIR}/02_data_cleaning_and_QC/data/Kat7_DualGuide_final_df.parquet"
working_dir="${PROJECT_DIR}/03_bootstrapping"

# command
python3 "$working_dir/main_code/UltraSeq_Bootstrapping_Multiplex.py" \
--a0 "$input_tumor_data_address" \
--a2 300 --a3 100 --a4 1000 --a5 KT --a6 KT \
--o1 "$working_dir/data/Kat7_DualGuide_BT" \
--o2 "$working_dir/data/Kat7_DualGuide_BT" \
--l1 50 60 70 80 90 95 96 97 98 99 \
--m 'N' --c 'No' \
--gplist gRNA_combination gRNA_combination_unordered gene_combination gene_combination_unordered

sacct --format=JobID,JobName,Submit,Start,End,State,Partition,ReqTRES%30,CPUTime,MaxRSS,NodeList%30 --units=M -j $SLURM_JOBID
#!/bin/bash
# 01-data_download.bash

#“#SBATCH” directives that convey submission options:
#SBATCH --job-name=Data_Download_HX_250217
#SBATCH --mail-user=xhq@stanford.edu
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=1g 
#SBATCH --time=2-00:00:00
#SBATCH --account=mwinslow
#SBATCH --partition=batch

# Source the configuration file 
source ../config.sh

# paramters for download data
# PORT:3022
# Account:X202SC25022737-Z01-F001
# Password:zt3n67ya

# mkdir directory for NGS data
mkdir -p $NGS_DIR

# download data
lftp -c "open sftp://X202SC25022737-Z01-F001:zt3n67ya@usftp23.novogene.com:3022; mirror --verbose --use-pget-n=8 -c . $NGS_DIR"

sacct --format=JobID,JobName,Submit,Start,End,State,Partition,ReqTRES%30,CPUTime,MaxRSS,NodeList%30 --units=M -j $SLURM_JOBID
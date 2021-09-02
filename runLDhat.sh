#!/usr/bin/bash

#-----------------------------Output files-----------------------------
#SBATCH --output=LD_output_%A.%a.txt
#SBATCH --error=LD_error_%A.%a.txt
#-----------------------------Required resources-----------------------
#SBATCH --time=0-200:00:00
#SBATCH --cpus-per-task=10 #10
#SBATCH --mem=50000 #50000
#----------------------------------------------------------------------

while getopts i:l:x:y:p:b: option
do
  case "${option}"
  in
    i) inputpref=${OPTARG};;
    l) lookup=${OPTARG};;
    x) iterations=${OPTARG};;
    y) sampling=${OPTARG};;
    p) prefix=${OPTARG};;
    b) blockpen=${OPTARG};;
  esac
done

idx=$SLURM_ARRAY_TASK_ID
echo -e "${inputpref}_${idx}\n\n"
LDhat/interval -seq ${inputpref}_${idx}.sites -loc ${inputpref}_${idx}.locs -lk $lookup -its $iterations -samp $sampling -bpen $blockpen -prefix ${prefix}_${idx}.

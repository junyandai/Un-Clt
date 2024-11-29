#!/bin/bash
#SBATCH --job-name=pairwise-distance
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --qos=highmem
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --constraint=EPYC-7313
#SBATCH --output=pairwise-distance-output_%j.log  
#SBATCH --error=pairwise-distance-error_%j.log  

module load Python3/3.9.15
character_file=$1
tree1_file=$2
tree2_file=$3
output_file=$4


python3 8_pairwise_distance.py -i $character_file -t1 $tree1_file -t2 $tree2_file -o $output_file
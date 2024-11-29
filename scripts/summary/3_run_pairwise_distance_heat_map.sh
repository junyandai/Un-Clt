#!/bin/bash
#SBATCH --job-name=pairwise-distance
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --qos=highmem
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --constraint=EPYC-7313
#SBATCH --output=pairwise-distance-output_%j.log  
#SBATCH --error=pairwise-distance-error_%j.log  

module load Python3/3.9.15

# python3 3_pairwise_distance_heat_map.py -i ../../result/3513_NT_T1_Fam/pairwise_distance.npy -o 3513_NT_T1_Fam_distance_heat_map.pdf
# python3 3_pairwise_distance_heat_map.py -i ../../result/3435_NT_T1/pairwise_distance.npy -o 3435_NT_T1_distance_heat_map.pdf
# python3 3_pairwise_distance_heat_map.py -i ../../result/3515_Lkb1_T1_Fam/pairwise_distance.npy -o 3515_Lkb1_T1_Fam_distance_heat_map.pdf
# python3 3_pairwise_distance_heat_map.py -i ../../result/3724_NT_All/pairwise_distance.npy -o 3724_NT_All_distance_heat_map.pdf
python3 3_pairwise_distance_heat_map.py -i ../../result/3726_NT_T1/pairwise_distance.npy -o 3726_NT_T1_distance_heat_map.pdf

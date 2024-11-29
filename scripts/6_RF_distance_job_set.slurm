#!/bin/bash
#SBATCH --job-name=RF-distance
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --qos=highmem
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --constraint=EPYC-7313
#SBATCH --output=RF-distance-output_%j.log  
#SBATCH --error=RF-distance-error_%j.log  

module load Python3/3.9.15

EXE="6_RF_distance.py"
T1=$1
T2=$2
OUTFILE=$3

python3 $EXE -t1 $T1 -t2 $T2 -o $OUTFILE
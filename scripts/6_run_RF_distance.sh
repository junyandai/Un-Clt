#!/bin/bash

# DATAS=("3513_NT_T1_Fam" \
#         "3724_NT_All"  \
#         "3515_Lkb1_T1_Fam" \
#         "3726_NT_T1" \
#         "3435_NT_T1")
DATAS=("3515_Lkb1_T1_Fam" \
        "3726_NT_T1" \
        "3435_NT_T1")

EXE="6_RF_distance.py"

for DATA in ${DATAS[@]}; do
    T1="../result/$DATA/star_cdp_one_sol.tre"
    T2="../result/$DATA/star_cdp_random_sol_trees.tre"
    OUTFILE="..//result/$DATA/RF_distances.csv"
    sbatch 6_RF_distance_job_set.slurm $T1 $T2 $OUTFILE
    
    done
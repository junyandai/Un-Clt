#!/bin/bash

export PROJDIR="/fs/cbcb-lab/ekmolloy/jdai123/un-clt-study"
DATAS=("3513_NT_T1_Fam" \
        "3724_NT_All"   \
        "3515_Lkb1_T1_Fam"  \
        "3435_NT_T1"    \
        "3726_NT_T1")

OUT="$PROJDIR/result"

for DATA in ${DATAS[@]}; do
    INPUT="$PROJDIR/data/$DATA/${DATA}_pruned_character_matrix.csv"
    TREE1="$PROJDIR/result/$DATA/star_cdp_one_sol.tre"
    TREE2="$PROJDIR/result/$DATA/star_cdp_random_sol_trees.tre"
    OUTPUT="$OUT/$DATA/pairwise_distance.npy"
    sbatch 8_pairwise_distance_job_set.slurm $INPUT $TREE1 $TREE2 $OUTPUT
    done
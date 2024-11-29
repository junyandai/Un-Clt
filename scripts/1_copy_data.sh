#!/bin/bash
exit
export PROJDIR="/fs/cbcb-lab/ekmolloy/jdai123/un-clt-study"

DATAS=("3513_NT_T1_Fam" \
        "3724_NT_All"  \
        "3515_Lkb1_T1_Fam" \
        "3726_NT_T1" \
        "3435_NT_T1")

OUT="$PROJDIR/data"
SOURCE="/fs/cbcb-lab/ekmolloy/jdai123/star-study/A_data/KPTracer-Data-Full"

for DATA in ${DATAS[@]}; do
    TARGET="$SOURCE/$DATA"
    cp -r $TARGET $OUT
    done

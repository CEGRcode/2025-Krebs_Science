#!/bin/bash

# Organize data from X_Bulk_Processing into Z_Figures for F6

### CHANGE ME
WRK=/storage/group/bfp2/default/hxc585_HainingChen/2025_Chen_TF-Nuc/Z_Figures
#WRK=/storage/home/owl5022/scratch/2024-Krebs_Science/Z_Figures
###

set -exo
module load anaconda3
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

LIBRARY=$WRK/../X_Bulk_Processing/Library
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar
VIOLIN=$WRK/../bin/make_violin_plot.py


[ -d S6 ] || mkdir S6

# ===============================================================================================================================
[ -d S6/a ] || mkdir S6/a 
cp $LIBRARY/NFIA_SORT-Occupancy_GROUP-Q4/Composites/BNase-seq_50U-10min_merge_hg38_NFIA_SORT-Occupancy_GROUP-Q4_5read1_TotalTag.out S6/a/

## need to add BNase-seq_50U-10min_merge_hg38_NFIA_SORT-Occupancy_GROUP-Q4_read1.out
# ===============================================================================================================================
[ -d S6/b ] || mkdir S6/b 
##heatmap
cp $LIBRARY/10phase/NFIA_SORT-Occupancy_GROUP-Q4/BNase-seq_50U-10min_merge_hg38_NFIA_SORT-Occupancy_GROUP-Q4_nearNuc_read1_shuffle_250bp_phase_sort_merge.svg S6/b/
cp $LIBRARY/10phase/NFIA_SORT-Occupancy_GROUP-Q4/BNase-seq_50U-10min_merge_hg38_NFIA_SORT-Occupancy_GROUP-Q4_nearNuc_read1_original_250bp_phase_sort_merge.svg S6/b/
# Composites
cp $LIBRARY/10phase/NFIA_SORT-Occupancy_GROUP-Q4/10xplot/BNase-seq_50U-10min_merge_hg38_NFIA_SORT-Occupancy_GROUP-Q4_nearNuc_read1_original_phase_4.out S6/b/
cp $LIBRARY/10phase/NFIA_SORT-Occupancy_GROUP-Q4/10xplot/BNase-seq_50U-10min_merge_hg38_NFIA_SORT-Occupancy_GROUP-Q4_nearNuc_read1_shuffle_phase_4.out S6/b/
cp $LIBRARY/10phase/NFIA_SORT-Occupancy_GROUP-Q4/10xplot/BNase-seq_50U-10min_merge_hg38_NFIA_SORT-Occupancy_GROUP-Q4_nearNuc_read1_original_phase_9.out S6/b/
cp $LIBRARY/10phase/NFIA_SORT-Occupancy_GROUP-Q4/10xplot/BNase-seq_50U-10min_merge_hg38_NFIA_SORT-Occupancy_GROUP-Q4_nearNuc_read1_shuffle_phase_9.out S6/b/
cp $LIBRARY/10phase/NFIA_SORT-Occupancy_GROUP-Q4/10xplot/WW_NFIA_nearNuc_original_phase_9.out S6/b/
cp $LIBRARY/10phase/NFIA_SORT-Occupancy_GROUP-Q4/10xplot/RR_NFIA_nearNuc_original_phase_9.out S6/b/
cp $LIBRARY/10phase/NFIA_SORT-Occupancy_GROUP-Q4/10xplot/YY_NFIA_nearNuc_original_phase_9.out S6/b/
cp $LIBRARY/10phase/NFIA_SORT-Occupancy_GROUP-Q4/10xplot/SS_NFIA_nearNuc_original_phase_9.out S6/b/


for file in S6/b/WW_*.out S6/b/RR_*.out S6/b/YY_*.out S6/b/SS_*.out; do
    filename=$(basename "$file" ".out")
    cut -f 1-481 "$file" > S6/b/${filename}_sense_left.out
    echo "$(seq -19 19 | tr '\n' ' ')" > S6/b/${filename}_sense_mid.out
    for i in $(seq 1 $(($(wc -l < "$file") - 1))); do
        echo -e "$(printf '\t%.0s-' {1..39})" >> S6/b/${filename}_sense_mid_placeholder.out
    done
    cat S6/b/${filename}_sense_mid.out S6/b/${filename}_sense_mid_placeholder.out > S6/b/${filename}_sense_mid_combined.out
    rm S6/b/${filename}_sense_mid.out S6/b/${filename}_sense_mid_placeholder.out
    cut -f 521-1001 "$file" > S6/b/${filename}_sense_right.out
    paste S6/b/${filename}_sense_left.out S6/b/${filename}_sense_mid_combined.out S6/b/${filename}_sense_right.out > S6/b/${filename}_sense_modified.out
    rm S6/b/${filename}_sense_left.out S6/b/${filename}_sense_mid_combined.out S6/b/${filename}_sense_right.out
done



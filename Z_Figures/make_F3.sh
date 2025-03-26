#!/bin/bash

# Organize data from X_Bulk_Processing into Z_Figures for F4

### CHANGE ME
WRK=/storage/group/bfp2/default/hxc585_HainingChen/2025_Chen_TF-Nuc/Z_Figures
#WRK=/storage/home/owl5022/scratch/2024-Krebs_Science/Z_Figures
###

set -exo
module load anaconda3
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

LIBRARY=$WRK/../X_Bulk_Processing/Library
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar

[ -d F3 ] || mkdir F3

# ===============================================================================================================================

[ -d F3/a ] || mkdir F3/a

# heatmap
cp $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile4/BNase-seq_50U-10min_merge_hg38_CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile4_+1Nuc_read1_shuffle_250bp_phase_sort_merge.svg F3/a/
cp $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile4/BNase-seQ4_50U-10min_merge_hg38_CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile4_+1Nuc_read1_original_250bp_phase_sort_merge.svg F3/a/
cp $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile4/BNase-seq_50U-10min_merge_hg38_CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile4_-1Nuc_read1_shuffle_250bp_phase_sort_merge.svg F3/a/
cp $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile4/BNase-seQ4_50U-10min_merge_hg38_CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile4_-1Nuc_read1_original_250bp_phase_sort_merge.svg F3/a/
# Composites
cp $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile4/10xplot/BNase-seq_50U-10min_merge_hg38_CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile4_-1Nuc_read1_shuffle_phase_2.out F3/a/
cp $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile4/10xplot/BNase-seq_50U-10min_merge_hg38_CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile4_-1Nuc_read1_original_phase_2.out F3/a/
cp $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile4/10xplot/BNase-seq_50U-10min_merge_hg38_CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile4_+1Nuc_read1_shuffle_phase_2.out F3/a/
cp $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile4/10xplot/BNase-seq_50U-10min_merge_hg38_CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile4_+1Nuc_read1_original_phase_2.out F3/a/

# ===============================================================================================================================

[ -d F3/b ] || mkdir F3/b

# heatmap
cp $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile1/BNase-seq_50U-10min_merge_hg38_CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile1_nearNuc_read1_shuffle_250bp_phase_sort_merge.svg F3/b/
cp $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile1/BNase-seq_50U-10min_merge_hg38_CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile1_nearNuc_read1_original_250bp_phase_sort_merge.svg F3/b/

# Composites
cp $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile1/10xplot/BNase-seq_50U-10min_merge_hg38_CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile1_nearNuc_read1_shuffle_phase_2.out F3/b/
cp $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile1/10xplot/BNase-seq_50U-10min_merge_hg38_CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile1_nearNuc_read1_shuffle_phase_7.out F3/b/
cp $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile1/10xplot/BNase-seq_50U-10min_merge_hg38_CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile1_nearNuc_read1_original_phase_2.out F3/b/
cp $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile1/10xplot/BNase-seq_50U-10min_merge_hg38_CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile1_nearNuc_read1_original_phase_7.out F3/b/


# ===============================================================================================================================
[ -d F3/c ] || mkdir F3/c

cp $LIBRARY/10phase/CTCF_phase_aligned/RR*.out F3/c/
cp $LIBRARY/10phase/CTCF_phase_aligned/YY*.out F3/c/
cp $LIBRARY/10phase/CTCF_phase_aligned/SS*.out F3/c/
cp $LIBRARY/10phase/CTCF_phase_aligned/WW*.out F3/c/

for file in F3/c/*.out; do
    filename=$(basename "$file" ".out")
    cut -f 1-231 "$file" > F3/c/${filename}_sense_left.out
    echo "$(seq -19 19 | tr '\n' ' ')" > F3/c/${filename}_sense_mid.out
    for i in $(seq 1 $(($(wc -l < "$file") - 1))); do
        echo -e "$(printf '\t%.0s-' {1..39})" >> F3/c/${filename}_sense_mid_placeholder.out
    done
    cat F3/c/${filename}_sense_mid.out F3/c/${filename}_sense_mid_placeholder.out > F3/c/${filename}_sense_mid_combined.out
    rm F3/c/${filename}_sense_mid.out F3/c/${filename}_sense_mid_placeholder.out
    cut -f 271-501 "$file" > F3/c/${filename}_sense_right.out
    paste F3/c/${filename}_sense_left.out F3/c/${filename}_sense_mid_combined.out F3/c/${filename}_sense_right.out > F3/c/${filename}_sense_modified.out
    rm F3/c/${filename}_sense_left.out F3/c/${filename}_sense_mid_combined.out F3/c/${filename}_sense_right.out
done



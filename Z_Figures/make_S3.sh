#!/bin/bash

# Organize data from X_Bulk_Processing into Z_Figures for S3c

### CHANGE ME
WRK=/storage/group/bfp2/default/hxc585_HainingChen/2025_Chen_TF-Nuc/Z_Figures
#WRK=/storage/home/owl5022/scratch/2024-Krebs_Science/Z_Figures
###

set -exo
module load anaconda3
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

LIBRARY=$WRK/../X_Bulk_Processing/Library
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar

[ -d S3 ] || mkdir S3

# ===============================================================================================================================

[ -d S3/a ] || mkdir S3/a

# Composites
cp $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile4/10xplot/BNase-seq_50U-10min_merge_hg38_CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile4_-1Nuc_read1_original_phase_*.out S3/a/
cp $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile4/10xplot/BNase-seq_50U-10min_merge_hg38_CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile4_+1Nuc_read1_original_phase_*.out S3/a/
cp $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile4/10xplot/RR*.out S3/a/
cp $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile4/10xplot/YY*.out S3/a/
cp $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile4/10xplot/WW*.out S3/a/
cp $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile4/10xplot/SS*.out S3/a/

for file in S3/a/RR*.out S3/a/YY*.out S3/a/WW*.out S3/a/SS*.out; do
    filename=$(basename "$file" ".out")
    cut -f 1-481 "$file" > S3/a/${filename}_sense_left.out
    echo "$(seq -19 19 | tr '\n' ' ')" > S3/a/${filename}_sense_mid.out
    for i in $(seq 1 $(($(wc -l < "$file") - 1))); do
        echo -e "$(printf '\t%.0s-' {1..39})" >> S3/a/${filename}_sense_mid_placeholder.out
    done
    cat S3/a/${filename}_sense_mid.out S3/a/${filename}_sense_mid_placeholder.out > S3/a/${filename}_sense_mid_combined.out
    rm S3/a/${filename}_sense_mid.out S3/a/${filename}_sense_mid_placeholder.out
    cut -f 521-1001 "$file" > S3/a/${filename}_sense_right.out
    paste S3/a/${filename}_sense_left.out S3/a/${filename}_sense_mid_combined.out S3/a/${filename}_sense_right.out > S3/a/${filename}_sense_modified.out
    rm S3/a/${filename}_sense_left.out S3/a/${filename}_sense_mid_combined.out S3/a/${filename}_sense_right.out
    rm $file
done

[ -d S3/b ] || mkdir S3/b

cp $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile1/10xplot/BNase-seq_50U-10min_merge_hg38_CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile1_nearNuc_read1_original_phase_*.out S3/b/
cp $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile1/10xplot/BNase-seq_50U-10min_merge_hg38_CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile1_nearNuc_read1_original_phase_*.out S3/b/
cp $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile1/10xplot/RR*.out S3/b/
cp $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile1/10xplot/YY*.out S3/b/
cp $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile1/10xplot/WW*.out S3/b/
cp $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile1/10xplot/SS*.out S3/b/


for file in S3/b/RR*.out S3/b/YY*.out S3/b/WW*.out S3/b/SS*.out; do
    filename=$(basename "$file" ".out")
    cut -f 1-481 "$file" > S3/b/${filename}_sense_left.out
    echo "$(seq -19 19 | tr '\n' ' ')" > S3/b/${filename}_sense_mid.out
    for i in $(seq 1 $(($(wc -l < "$file") - 1))); do
        echo -e "$(printf '\t%.0s-' {1..39})" >> S3/b/${filename}_sense_mid_placeholder.out
    done
    cat S3/b/${filename}_sense_mid.out S3/b/${filename}_sense_mid_placeholder.out > S3/b/${filename}_sense_mid_combined.out
    rm S3/b/${filename}_sense_mid.out S3/b/${filename}_sense_mid_placeholder.out
    cut -f 521-1001 "$file" > S3/b/${filename}_sense_right.out
    paste S3/b/${filename}_sense_left.out S3/b/${filename}_sense_mid_combined.out S3/b/${filename}_sense_right.out > S3/b/${filename}_sense_modified.out
    rm S3/b/${filename}_sense_left.out S3/b/${filename}_sense_mid_combined.out S3/b/${filename}_sense_right.out
    rm $file
done

[ -d S3/c ] || mkdir S3/c

cp $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile1_-1Nuc/10plot/BNase-seq_50U-10min_merge_hg38_CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile1_-1Nuc_1000bp_read1_original_phase_*.out S3/c/
cp $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile1_+1Nuc/10plot/BNase-seq_50U-10min_merge_hg38_CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile1_+1Nuc_1000bp_read1_original_phase_*.out S3/c/
cp $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile1_-1Nuc/10plot/RR*.out S3/c/
cp $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile1_-1Nuc/10plot/YY*.out S3/c/
cp $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile1_-1Nuc/10plot/WW*.out S3/c/
cp $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile1_-1Nuc/10plot/SS*.out S3/c/
cp $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile1_+1Nuc/10plot/RR*.out S3/c/
cp $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile1_+1Nuc/10plot/YY*.out S3/c/
cp $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile1_+1Nuc/10plot/WW*.out S3/c/
cp $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile1_+1Nuc/10plot/SS*.out S3/c/

for file in S3/c/RR*.out S3/c/YY*.out S3/c/WW*.out S3/c/SS*.out; do
    filename=$(basename "$file" ".out")
    cut -f 1-481 "$file" > S3/c/${filename}_sense_left.out
    echo "$(seq -19 19 | tr '\n' ' ')" > S3/c/${filename}_sense_mid.out
    for i in $(seq 1 $(($(wc -l < "$file") - 1))); do
        echo -e "$(printf '\t%.0s-' {1..39})" >> S3/c/${filename}_sense_mid_placeholder.out
    done
    cat S3/c/${filename}_sense_mid.out S3/c/${filename}_sense_mid_placeholder.out > S3/c/${filename}_sense_mid_combined.out
    rm S3/c/${filename}_sense_mid.out S3/c/${filename}_sense_mid_placeholder.out
    cut -f 521-1001 "$file" > S3/c/${filename}_sense_right.out
    paste S3/c/${filename}_sense_left.out S3/c/${filename}_sense_mid_combined.out S3/c/${filename}_sense_right.out > S3/c/${filename}_sense_modified.out
    rm S3/c/${filename}_sense_left.out S3/c/${filename}_sense_mid_combined.out S3/c/${filename}_sense_right.out
    rm $file
done


[ -d S3/d ] || mkdir S3/d
cp $LIBRARY/DNAshape/ATF7_MA0834.1/Composites/*.out S3/d
cp $LIBRARY/DNAshape/CTCF_MA1929.1/Composites/*.out S3/d
cp $LIBRARY/DNAshape/ZKSCAN1_MA1585.1/Composites/*.out S3/d
cp $LIBRARY/DNAshape/ESRRA_MA0592.1/Composites/*.out S3/d
cp $LIBRARY/DNAshape/JUND_MA0492.2/Composites/*.out S3/d
cp $LIBRARY/DNAshape/NRF1_MA0506.3/Composites/*.out S3/d
cp $LIBRARY/DNAshape/REST_MA0138.3/Composites/*.out S3/d
cp $LIBRARY/DNAshape/SP1_MA0079.5/Composites/*.out S3/d
cp $LIBRARY/DNAshape/SPI1_MA0080.7/Composites/*.out S3/d
cp $LIBRARY/DNAshape/ZNF263_MA0528.1/Composites/*.out S3/d
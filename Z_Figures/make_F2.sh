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

[ -d F2 ] || mkdir F2

# ===============================================================================================================================

cp $LIBRARY/WebLogos/CTCF_MA1929.1_logo.eps F2/
cp $LIBRARY/WebLogos/ZKSCAN1_MA1585.1_logo.eps F2/
cp $LIBRARY/WebLogos/REST_MA0138.3_logo.eps F2/
cp $LIBRARY/WebLogos/ESRRA_MA0592.1_logo.eps F2/
cp $LIBRARY/WebLogos/ATF7_MA0834.1_logo.eps F2/
cp $LIBRARY/WebLogos/JUND_MA0492.2_logo.eps F2/
cp $LIBRARY/WebLogos/NRF1_MA0506.3_logo.eps F2/
cp $LIBRARY/WebLogos/SP1_MA0079.5_logoRC.eps F2/
cp $LIBRARY/WebLogos/SPI1_MA0080.7_logo.eps F2/
cp $LIBRARY/WebLogos/ZNF263_MA0528.1_logo.eps F2/

# Composites
cp $LIBRARY/BI_Pileups/MA1929/Composites/*.out F2/
cp $LIBRARY/BI_Pileups/MA1585/Composites/*.out F2/
cp $LIBRARY/BI_Pileups/MA0138/Composites/*.out F2/
cp $LIBRARY/BI_Pileups/MA0592/Composites/*.out F2/
cp $LIBRARY/BI_Pileups/MA0834/Composites/*.out F2/
cp $LIBRARY/BI_Pileups/MA0492/Composites/*.out F2/
cp $LIBRARY/BI_Pileups/MA0509/Composites/*.out F2/
cp $LIBRARY//BI_Pileups/MA0079/Composites/*.out F2/
cp $LIBRARY/BI_Pileups/MA0080/Composites/*.out F2/
cp $LIBRARY/BI_Pileups/MA0528/Composites/*.out F2/

# heatmap

cat $LIBRARY/10phase/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile1/CTCF_nearNuc_Q1_original_all.bed | sort -k5,5n > F2/CTCF_nearNuc_Q1_original_Nucsort.bed

java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 500 F2/CTCF_nearNuc_Q1_original_Nucsort.bed -o F2/CTCF_nearNuc_Q1_original_Nucsort_500bp.bed
java -jar $SCRIPTMANAGER read-analysis tag-pileup F2/CTCF_nearNuc_Q1_original_Nucsort_500bp.bed $WRK/../data/BAM/BNase-seq_50U-10min_merge_hg38.bam --cpu 4 -m -M F2/BNase-seq_50U-10min_merge_hg38_CTCF_nearNuc_Q1_original_Nucsort_500bp_midpoint
java -jar $SCRIPTMANAGER figure-generation heatmap -p .95 -c B8C400 F2/BNase-seq_50U-10min_merge_hg38_CTCF_nearNuc_Q1_original_Nucsort_500bp_midpoint_combined.cdt -o  F2/BNase-seq_50U-10min_merge_hg38_CTCF_nearNuc_Q1_original_Nucsort_500bp_midpoint.png
java -jar $SCRIPTMANAGER figure-generation label-heatmap F2/BNase-seq_50U-10min_merge_hg38_CTCF_nearNuc_Q1_original_Nucsort_500bp_midpoint.png -f 20 -l -250 -m 0 -r +250  -o F2/BNase-seq_50U-10min_merge_hg38_CTCF_nearNuc_Q1_original_Nucsort_500bp_midpoint.svg
rm F2/CTCF_nearNuc_Q1_original_Nucsort_500bp.bed
rm F2/BNase-seq_50U-10min_merge_hg38_CTCF_nearNuc_Q1_original_Nucsort_500bp_midpoint.png
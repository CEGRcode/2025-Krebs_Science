#!/bin/bash

# Organize data from X_Bulk_Processing into Z_Figures for F4

### CHANGE ME
WRK=/path/to/2024-Chen_Nature/Z_Figures
WRK=/storage/home/owl5022/scratch/2024-Chen_Nature/Z_Figures
###

set -exo
module load anaconda3
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

LIBRARY=$WRK/../X_Bulk_Processing/Library
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar

BED=CTCF_SORT-Occupancy_1000bp

[ -d F4 ] || mkdir F4

# ===============================================================================================================================

[ -d F4/a ] || mkdir F4/a

cp $LIBRARY/WebLogos/CTCF_M1_logo.eps F4/a/

# Composites
cp $LIBRARY/$BED/Composites/K562_CTCF_BX_rep1_hg38_${BED}_5read1_NCIS.out F4/a/
cp $LIBRARY/$BED/Composites/K562_CTCF_BX_rep1_hg38_${BED}_5read2_NCIS.out F4/a/
cp $LIBRARY/$BED/Composites/BNase-seq_50U-10min_merge_hg38_${BED}_midpoint_combined.out F4/a/

# ===============================================================================================================================

[ -d F4/b ] || mkdir F4/b

cp $LIBRARY/$BED/FourColor/${BED}_31bp.svg F4/b/

# Heatmaps
cp $LIBRARY/$BED/SVG/K562_IgG_BX_merge_hg38_${BED}_5read1_Raw_merge_label.svg F4/b/
cp $LIBRARY/$BED/SVG/K562_CTCF_BX_rep1_hg38_${BED}_5read1_NCIS_merge_label.svg F4/b/
cp $LIBRARY/$BED/SVG/K562_CTCF_BX_rep1_hg38_${BED}_5read2_NCIS_merge_label.svg F4/b/
cp $LIBRARY/$BED/SVG/K562_CTCF_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS_merge_label.svg F4/b/
cp $LIBRARY/$BED/SVG/K562_CTCF_BX_rep1_hg38_${BED}_5read2-MIN100_NCIS_merge_label.svg F4/b/
cp $LIBRARY/$BED/SVG/BNase-seq_50U-10min_merge_hg38_${BED}_midpoint-MIN100_TotalTag_combined.svg F4/b/

# ===============================================================================================================================

[ -d F4/c ] || mkdir F4/c

# Composites
cp $LIBRARY/$BED/Composites/K562_CTCF_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS.out F4/c/
cp $LIBRARY/$BED/Composites/K562_CTCF_BX_rep1_hg38_${BED}_5read2-MIN100_NCIS.out F4/c/
cp $LIBRARY/$BED/Composites/BNase-seq_50U-10min_merge_hg38_${BED}_midpoint-MIN100_TotalTag_combined.out F4/c/

# ===============================================================================================================================

[ -d F4/d ] || mkdir F4/d

# Composites
cp $LIBRARY/$BED/Composites/K562_CTCF_BX_rep1_hg38_${BED}_5read1_NCIS.out F4/d/CTCF_K562_CTCF_BX_rep1_hg38_${BED}_5read1_NCIS.out
cp $LIBRARY/$BED/Composites/K562_RAD21_BX_rep1_hg38_${BED}_5read1_NCIS.out F4/d/RAD21_K562_RAD21_BX_rep1_hg38_${BED}_5read1_NCIS.out
cp $LIBRARY/$BED/Composites/K562_SMC3_BX_rep1_hg38_${BED}_5read1_NCIS.out F4/d/SMC3_K562_SMC3_BX_rep1_hg38_${BED}_5read1_NCIS.out

# ===============================================================================================================================

[ -d F4/e ] || mkdir F4/e

# Composites
cp $LIBRARY/$BED/Composites/K562_RAD21_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS.out F4/e/RAD21_K562_RAD21_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS.out
cp $LIBRARY/$BED/Composites/K562_RAD21_BX_rep1_hg38_${BED}_5read2-MIN100_NCIS.out F4/e/RAD21_K562_RAD21_BX_rep1_hg38_${BED}_5read2-MIN100_NCIS.out
cp $LIBRARY/$BED/Composites/K562_SMC3_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS.out F4/e/SMC3_K562_SMC3_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS.out
cp $LIBRARY/$BED/Composites/K562_SMC3_BX_rep1_hg38_${BED}_5read2-MIN100_NCIS.out F4/e/SMC3_K562_SMC3_BX_rep1_hg38_${BED}_5read2-MIN100_NCIS.out
cp $LIBRARY/$BED/Composites/BNase-seq_50U-10min_merge_hg38_${BED}_5read2-MIN100_TotalTag.out F4/e/

# ===============================================================================================================================

[ -d F4/f ] || mkdir F4/f

# Composites
cp $LIBRARY/$BED/Composites/K562_CTCF_BX_rep1_hg38_${BED}_5read2-MIN100_NCIS.out F4/f/CTCF_K562_CTCF_BX_rep1_hg38_${BED}_5read2-MIN100_NCIS.out
cp $LIBRARY/$BED/Composites/K562_RAD21_BX_rep1_hg38_${BED}_5read2-MIN100_NCIS.out F4/f/RAD21_K562_RAD21_BX_rep1_hg38_${BED}_5read2-MIN100_NCIS.out
cp $LIBRARY/$BED/Composites/K562_SMC3_BX_rep1_hg38_${BED}_5read2-MIN100_NCIS.out F4/f/SMC3_K562_SMC3_BX_rep1_hg38_${BED}_5read2-MIN100_NCIS.out
cp $LIBRARY/$BED/Composites/BNase-seq_50U-10min_merge_hg38_${BED}_5read2-MIN100_TotalTag.out F4/f/
# HAINING QUESTION: Are these filtered by min=100?
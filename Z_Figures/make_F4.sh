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

BED=CTCF_SORT-Occupancy_1000bp

[ -d F4 ] || mkdir F4

# ===============================================================================================================================

[ -d F4/B ] || mkdir F4/B

cp $LIBRARY/WebLogos/CTCF_M1_logo.eps F4/B/

# Composites
cp $LIBRARY/$BED/Composites/K562_CTCF_BX_rep1_hg38_${BED}_5read1_NCIS.out F4/B/
cp $LIBRARY/$BED/Composites/K562_CTCF_BX_rep1_hg38_${BED}_5read2_NCIS.out F4/B/
cp $LIBRARY/$BED/Composites/BNase-seq_50U-10min_merge_hg38_${BED}_midpoint_TotalTag_combined.out F4/B/

# ===============================================================================================================================

[ -d F4/C ] || mkdir F4/C

cp $LIBRARY/$BED/FourColor/${BED}_31bp.svg F4/C/

# Heatmaps
cp $LIBRARY/$BED/SVG/K562_IgG_BX_merge_hg38_${BED}_5read1_Raw_merge_label.svg F4/C/
cp $LIBRARY/$BED/SVG/K562_CTCF_BX_rep1_hg38_${BED}_5read1_NCIS_merge_label.svg F4/C/
cp $LIBRARY/$BED/SVG/K562_CTCF_BX_rep1_hg38_${BED}_5read2_NCIS_merge_label.svg F4/C/
cp $LIBRARY/$BED/SVG/K562_CTCF_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS_merge_label.svg F4/C/
cp $LIBRARY/$BED/SVG/K562_CTCF_BX_rep1_hg38_${BED}_5read2-MIN100_NCIS_merge_label.svg F4/C/
cp $LIBRARY/$BED/SVG/BNase-seq_50U-10min_merge_hg38_${BED}_midpoint_TotalTag_combined.svg F4/C/

# ===============================================================================================================================

[ -d F4/D ] || mkdir F4/D

# Composites
cp $LIBRARY/$BED/Composites/K562_CTCF_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS.out F4/D/
cp $LIBRARY/$BED/Composites/K562_CTCF_BX_rep1_hg38_${BED}_5read2-MIN100_NCIS.out F4/D/
cp $LIBRARY/$BED/Composites/BNase-seq_50U-10min_merge_hg38_${BED}_midpoint_TotalTag_combined.out F4/D/

# ===============================================================================================================================

[ -d F4/E ] || mkdir F4/E

# Composites
cp $LIBRARY/$BED/Composites/K562_CTCF_BX_rep1_hg38_${BED}_5read1_NCIS.out F4/E/CTCF_K562_CTCF_BX_rep1_hg38_${BED}_5read1_NCIS.out
cp $LIBRARY/$BED/Composites/K562_RAD21_BX_rep1_hg38_${BED}_5read1_NCIS.out F4/E/RAD21_K562_RAD21_BX_rep1_hg38_${BED}_5read1_NCIS.out
cp $LIBRARY/$BED/Composites/K562_SMC3_BX_rep1_hg38_${BED}_5read1_NCIS.out F4/E/SMC3_K562_SMC3_BX_rep1_hg38_${BED}_5read1_NCIS.out

# ===============================================================================================================================

[ -d F4/F ] || mkdir F4/F

# Composites
cp $LIBRARY/$BED/Composites/K562_RAD21_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS.out F4/F/RAD21_K562_RAD21_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS.out
cp $LIBRARY/$BED/Composites/K562_RAD21_BX_rep1_hg38_${BED}_5read2-MIN100_NCIS.out F4/F/RAD21_K562_RAD21_BX_rep1_hg38_${BED}_5read2-MIN100_NCIS.out
cp $LIBRARY/$BED/Composites/K562_SMC3_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS.out F4/F/SMC3_K562_SMC3_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS.out
cp $LIBRARY/$BED/Composites/K562_SMC3_BX_rep1_hg38_${BED}_5read2-MIN100_NCIS.out F4/F/SMC3_K562_SMC3_BX_rep1_hg38_${BED}_5read2-MIN100_NCIS.out
cp $LIBRARY/$BED/Composites/BNase-seq_50U-10min_merge_hg38_${BED}_midpoint_TotalTag_combined.out F4/F/

# ===============================================================================================================================

[ -d F4/G ] || mkdir F4/G

# Composites
cp $LIBRARY/$BED/Composites/BNase-seq_50U-10min_merge_hg38_${BED}_5read1_TotalTag.out F4/G/
cp $LIBRARY/$BED/Composites/K562_CTCF_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS.out F4/G/
cp $LIBRARY/$BED/Composites/K562_RAD21_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS.out F4/G/RAD21_K562_RAD21_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS.out
cp $LIBRARY/$BED/Composites/K562_SMC3_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS.out F4/G/SMC3_K562_SMC3_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS.out



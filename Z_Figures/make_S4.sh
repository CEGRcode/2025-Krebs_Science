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
BED1=CTCF_SORT-SMC3Engagement_1000
BED2=CTCF_SORT-SMC3Engagement_GROUP-High_1000bp
BED3=CTCF_SORT-SMC3Engagement_GROUP-Low_1000bp

[ -d S4/a ] || mkdir S4/a

# ===============================================================================================================================

# Composites
cp $LIBRARY/$BED/Composites/K562_CTCF_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS.out S4/a/
## tag pile up by nonxl CTCF
java -jar $SCRIPTMANAGER read-analysis tag-pileup ../data/RefPT-Motif/1000bp/CTCF_SORT-Occupancy_1000bp.bed ../data/BAM/K562_CTCF_NonXLBX_rep1_hg38.bam -n 100 --cpu 4 -o S4/a/K562_CTCF_NonXLBX_rep1_hg38_CTCF_SORT-Occupancy_1000bp_5read1-MIN100.out

# ===============================================================================================================================

[ -d F4/b ] || mkdir F4/b

# Heatmaps

cp $LIBRARY/$BED1/SVG/K562_CTCF_BX_rep1_hg38_${BED1}_5read1-MIN100_NCIS_merge_label.svg F4/b/
cp $LIBRARY/$BED1/SVG/K562_RAD21_BX_rep1_hg38_${BED1}_5read1-MIN100_NCIS_merge_label.svg F4/b/
cp $LIBRARY/$BED1/SVG/K562_SMC3_BX_rep1_hg38_${BED1}_5read1-MIN100_NCIS_merge_label.svg F4/b/

# Composites
cp $LIBRARY/$BED2/Composites/K562_CTCF_BX_rep1_hg38_${BED2}_5read1-MIN100_NCIS.out S4/b/
cp $LIBRARY/$BED2/Composites/K562_RAD21_BX_rep1_hg38_${BED2}_5read1-MIN100_NCIS.out S4/b/
cp $LIBRARY/$BED2/Composites/K562_SMC3_BX_rep1_hg38_${BED2}_5read1-MIN100_NCIS.out S4/b/
cp $LIBRARY/$BED3/Composites/K562_CTCF_BX_rep1_hg38_${BED3}_5read1-MIN100_NCIS.out S4/b/
cp $LIBRARY/$BED3/Composites/K562_RAD21_BX_rep1_hg38_${BED3}_5read1-MIN100_NCIS.out S4/b/
cp $LIBRARY/$BED3/Composites/K562_SMC3_BX_rep1_hg38_${BED3}_5read1-MIN100_NCIS.out S4/b/

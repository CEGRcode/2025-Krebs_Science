#!/bin/bash

# Organize data from X_Bulk_Processing into Z_Figures for S3c

### CHANGE ME
WRK=/path/to/2024-Krebs_Science/Z_Figures
WRK=/storage/home/owl5022/scratch/2024-Krebs_Science/Z_Figures
###

set -exo
module load anaconda3
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

LIBRARY=$WRK/../X_Bulk_Processing/Library
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar

[ -d S3 ] || mkdir S3

# ===============================================================================================================================

[ -d S3/c ] || mkdir S3/c

cp $LIBRARY/WebLogos/ZKSCAN1_M1_logoRC.eps S3/c

# Composites
BED=ZKSCAN1_SORT-Occupancy_1000bp
cp $LIBRARY/$BED/Composites/K562_ZKSCAN1_BX_rep1_hg38_$BED\_5read1-MIN100_NCIS.out S3/c
cp $LIBRARY/$BED/Composites/K562_ZKSCAN1_BX_rep1_hg38_$BED\_5read2-MIN100_NCIS.out S3/c
cp $LIBRARY/$BED/Composites/BNase-seq_50U-10min_merge_hg38_$BED\_5read2-MIN100_TotalTag.out S3/c
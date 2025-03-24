#!/bin/bash

# Format CTCF reference point file to proper BED files.

# data/RefPT-Motif
#   |--CTCF_SORT-Occupancy.bed
#   |--CTCF_SORT-SMC3Engagement.bed
#   |--CTCF_SORT-SMC3Engagement_GROUP-Low.bed
#   |--CTCF_SORT-SMC3Engagement_GROUP-High.bed
#   |--1000bp
#     |--CTCF_SORT-Occupancy_1000bp.bed
#     |--CTCF_SORT-SMC3Engagement_GROUP-Low_1000bp.bed
#     |--CTCF_SORT-SMC3Engagement_GROUP-High_1000bp.bed
#     |--CTCF_SORT-SMC3Engagement_1000bp.bed
#   |--1bp
#     |--CTCF_SORT-Occupancy_1bp.bed
#     |--CTCF_SORT-SMC3Engagement_GROUP-Low_1bp.bed
#     |--CTCF_SORT-SMC3Engagement_GROUP-High_1bp.bed

### CHANGE ME
WRK=/path/to/2024-Krebs_Science/04_Call_Motifs
WRK=/scratch/owl5022/2024-Krebs_Science/04_Call_Motifs
GENOME=$WRK/../data/hg38_files/hg38.fa
GINFO=$WRK/../data/hg38_files/hg38.chrom.sizes
SMC3BAMFILE=$WRK/../data/BAM/K562_SMC3_BX_rep1_hg38.bam
###

# Dependencies
# - java

set -exo
module load anaconda
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Inputs and outputs
MOTIF=$WRK/../data/RefPT-Motif
CTCF_BOUND=$WRK/temp-3_Filter_and_Sort_by_occupancy/CTCF_CTCF-K562_M1_100bp_7-Occupancy_BOUND_1bp.bed

[ -d $MOTIF ] || mkdir $MOTIF
[ -d $MOTIF/1000bp ] || mkdir $MOTIF/1000bp
[ -d $MOTIF/1bp ] || mkdir $MOTIF/1bp

# Script shortcuts
SCRIPTMANAGER=../bin/ScriptManager-v0.15.jar

TEMP=temp-5_CTCF
[ -d $TEMP ] || mkdir $TEMP

# Reformat BED with SCORE=tag pileup tag count (move 7th column to score column and add K562_CTCF in the seventh)
cut -f 1-6 $CTCF_BOUND > $MOTIF/CTCF_SORT-Occupancy.bed

# QC: Stat lines
wc -l $MOTIF/CTCF_SORT-Occupancy.bed
#18838 ../data/RefPT-Motif/CTCF_SORT-Occupancy.bed
## Resort CTCF ref by SMC3-downstream nucleosome engagement
# make the region 300bp downstream from CTCF motif
bedtools shift -i $MOTIF/CTCF_SORT-Occupancy.bed -g $GINFO -p 150 -m -150 > $TEMP/CTCF-d150_SORT-Occupancy.bed

java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 300 $TEMP/CTCF-d150_SORT-Occupancy.bed -o $TEMP/CTCF-d150_SORT-Occupancy_300bp.bed
# sum the SMC3 largefragments read1 
java -jar "$SCRIPTMANAGER" read-analysis tag-pileup  $TEMP/CTCF-d150_SORT-Occupancy_300bp.bed $SMC3BAMFILE -1 -n 100 --combined --cpu 4 -M  $TEMP/SMC3_CTCF-d150_SORT-Occupancy_300bp_read1

java -jar "$SCRIPTMANAGER" read-analysis aggregate-data --sum $TEMP/SMC3_CTCF-d150_SORT-Occupancy_300bp_read1_combined.cdt -o $TEMP/SMC3_CTCF-d150_SORT-Occupancy_300bp_read1_combined.out
# attatch SMC3 score to original CTCF ref and resort it
cut -f 2 $TEMP/SMC3_CTCF-d150_SORT-Occupancy_300bp_read1_combined.out | tail -n +2  | cut -f 2 | \
paste $MOTIF/CTCF_SORT-Occupancy.bed - | bedtools sort -i | sort -k7,7nr | cut -f 1-6 > $MOTIF/CTCF_SORT-SMC3Engagement.bed
# take cohesin high engagment and low engagment
tail -n +14128 $MOTIF/CTCF_SORT-SMC3Engagement.bed > $MOTIF/CTCF_SORT-SMC3Engagement_GROUP-Low.bed
head -n 4709 $MOTIF/CTCF_SORT-SMC3Engagement.bed > $MOTIF/CTCF_SORT-SMC3Engagement_GROUP-High.bed

# Expand 1000bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/CTCF_SORT-Occupancy.bed -o $MOTIF/1000bp/CTCF_SORT-Occupancy_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/CTCF_SORT-SMC3Engagement.bed -o $MOTIF/1000bp/CTCF_SORT-SMC3Engagement_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/CTCF_SORT-SMC3Engagement_GROUP-Low.bed -o $MOTIF/1000bp/CTCF_SORT-SMC3Engagement_GROUP-Low_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/CTCF_SORT-SMC3Engagement_GROUP-High.bed -o $MOTIF/1000bp/CTCF_SORT-SMC3Engagement_GROUP-High_1000bp.bed

# Expand 1bp for rotational phasing quantification
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1 $MOTIF/CTCF_SORT-Occupancy.bed -o $MOTIF/1bp/CTCF_SORT-Occupancy_1bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1 $MOTIF/CTCF_SORT-SMC3Engagement_GROUP-Low.bed -o $MOTIF/1bp/CTCF_SORT-SMC3Engagement_GROUP-Low_1bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1 $MOTIF/CTCF_SORT-SMC3Engagement_GROUP-High.bed -o $MOTIF/1bp/CTCF_SORT-SMC3Engagement_GROUP-High_1bp.bed

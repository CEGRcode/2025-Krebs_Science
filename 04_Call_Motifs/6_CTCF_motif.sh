#!/bin/bash

# Format CTCF reference point file to proper BED files.

# data/RefPT-Motif
#   |--CTCF_SORT-Occupancy.bed
#   |--1000bp
#     |--CTCF_SORT-Occupancy_1000bp.bed

### CHANGE ME
WRK=/path/to/2024-Krebs_Science/04_Call_Motifs
WRK=/scratch/owl5022/2024-Krebs_Science/04_Call_Motifs
GENOME=../data/hg38_files/hg38.fa
GINFO=../data/hg38_files/hg38.chrom.sizes
SMC3BAMFILE=../data/BAM/K562_SMC3_BX_rep1_hg38.bam
###

# Dependencies
# - java

set -exo
module load anaconda
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Inputs and outputs
MOTIF=$WRK/../data/RefPT-Motif
CTCF_BOUND=temp-3_Filter_and_Sort_by_occupancy/CTCF_CTCF-K562_M1_100bp_7-Occupancy_BOUND_1bp.bed

[ -d $MOTIF ] || mkdir $MOTIF
[ -d $MOTIF/1000bp ] || mkdir $MOTIF/1000bp
[ -d $MOTIF/1bp ] || mkdir $MOTIF/1bp

# Script shortcuts
SCRIPTMANAGER=../bin/ScriptManager-v0.15.jar

TEMP=temp-5_CTCF
[ -d $TEMP ] || mkdir $TEMP

# Reformat BED with SCORE=tag pileup tag count (move 7th column to score column and add K562_CTCF in the seventh)
cut -f1-6 $CTCF_BOUND > $MOTIF/CTCF_SORT-Occupancy.bed

# QC: Stat lines
wc -l $MOTIF/CTCF_SORT-Occupancy.bed

# Expand 1000bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/CTCF_SORT-Occupancy.bed -o $MOTIF/1000bp/CTCF_SORT-Occupancy_1000bp.bed

## Resort CTCF ref by SMC3-downstream nucleosome engagement
# make the region 300bp downstream from CTCF motif
bedtools shift -i $MOTIF/CTCF_SORT-Occupancy.bed -g $GINFO -p 150 -m -150 > $TEMP/CTCF-d150_SORT-Occupancy.bed

java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 300 $TEMP/CTCF-d150_SORT-Occupancy.bed -o $TEMP/CTCF-d150_SORT-Occupancy_300bp.bed
# sum the SMC3 largefragments read1 
java -jar "$SCRIPTMANAGER" read-analysis tag-pileup  $TEMP/CTCF-d150_SORT-Occupancy_300bp.bed $SMC3BAMFILE -1 --combined --cpu 4 -M  $TEMP/SMC3_CTCF-d150_SORT-Occupancy_300bp_read1

java -jar "$SCRIPTMANAGER" read-analysis aggregate-data --sum $TEMP/SMC3_CTCF-d150_SORT-Occupancy_300bp_read1_combined.cdt -o $TEMP/SMC3_CTCF-d150_SORT-Occupancy_300bp_read1_combined.out
# attatch SMC3 score to original CTCF ref and resort it
cut -f 2 $TEMP/SMC3_CTCF-d150_SORT-Occupancy_300bp_read1_combined.out | tail -n +2  | cut -f 2 | \
paste CTCF_Occupancy_1bp.bed - | bedtools sort -i | sort -k8,8nr | cut -f 1-7 > CTCF_Downengagesort.bed

tail -n +13486 CTCF_Downengagesort.bed > $MOTIF/1bp/CTCF_Downengagesort_bottom_1bp.bed

head -n 4495 CTCF_Downengagesort.bed > $MOTIF/1bp/CTCF_Downengagesort_top_1bp.bed

java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 CTCF_Downengagesort.bed -o $MOTIF/1000bp/CTCF_Downengagesort_1000bp.bed

tail -n +13486 $MOTIF/1000bp/CTCF_Downengagesort_1000bp.bed > $MOTIF/1000bp/CTCF_Downengagesort_bottom_1000bp.bed

head -n 4495 $MOTIF/1000bp/CTCF_Downengagesort_1000bp.bed > $MOTIF/1000bp/CTCF_Downengagesort_top_1000bp.bed


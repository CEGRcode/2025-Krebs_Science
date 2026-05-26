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

# Dependencies
# - java

set -exo
source activate bx

# Inputs and outputs
MOTIF=../data/RefPT-Motif
CTCF_BOUND=temp-3_Filter_and_Sort_by_occupancy/CTCF_CTCF-K562_M1_100bp_7-Occupancy_BOUND_1bp.bed
GENOME=../data/hg38_files/hg38.fa
GINFO=../data/hg38_files/hg38.chrom.sizes
SMC3BAMFILE=../data/BAM/K562_SMC3_BX_rep1_hg38.bam

[ -d $MOTIF ] || mkdir $MOTIF
[ -d $MOTIF/1000bp ] || mkdir $MOTIF/1000bp
[ -d $MOTIF/1bp ] || mkdir $MOTIF/1bp

# Script shortcuts
SCRIPTMANAGER=../bin/ScriptManager-v0.15.jar

TEMP=temp-6_CTCF_motif
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

## copy CTCF ENCODE dataset for nucleosome phasing analysis
bedtools sort -i  ../data/RefPT-JASPAR/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile1.bed | bedtools shift -i - -g $GINFO -p 9 -m -9  >  $TEMP/CTCF_Q1.bed
bedtools sort -i  ../data/RefPT-JASPAR/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile4.bed | bedtools shift -i - -g $GINFO -p 9 -m -9  > $TEMP/CTCF_Q4.bed
bedtools sort -i  ../data/RefPT-JASPAR/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile3.bed | bedtools shift -i - -g $GINFO -p 9 -m -9  >  $TEMP/CTCF_Q3.bed
bedtools sort -i  ../data/RefPT-JASPAR/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile2.bed | bedtools shift -i - -g $GINFO -p 9 -m -9  > $TEMP/CTCF_Q2.bed

# find closest nucleosome to CTCF sites:

java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1 $TEMP/CTCF_Q1.bed  -o $TEMP/CTCF_Q1_1bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1 $TEMP/CTCF_Q4.bed  -o $TEMP/CTCF_Q4_1bp.bed

bedtools sort  -i $TEMP/CTCF_Q1_1bp.bed  |  bedtools closest -d -D a -t first -a - -b $TEMP/../temp-4_NFIA/Nucleosomes_SORT-Genomic.bed | sort -k13,13n > $TEMP/CTCF_Q1-ClosestDyad.tsv
bedtools sort  -i $TEMP/CTCF_Q4_1bp.bed  | bedtools closest -d -D a -t first -a - -b $TEMP/../temp-4_NFIA/Nucleosomes_SORT-Genomic.bed -iu  |  sort -k13,13n > $TEMP/CTCF_Q4-ClosestDyad-right.tsv
bedtools sort  -i $TEMP/CTCF_Q4_1bp.bed  | bedtools closest -d -D a -t first -a - -b $TEMP/../temp-4_NFIA/Nucleosomes_SORT-Genomic.bed -id  |  sort -k13,13n > $TEMP/CTCF_Q4-ClosestDyad-left.tsv


cut -f 1-6 $TEMP/CTCF_Q1-ClosestDyad.tsv > $MOTIF/CTCF_Q1-ClosestDyad.bed
sort -k5,5nr $TEMP/CTCF_Q4_1bp.bed > $MOTIF/CTCF_Q4.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 500 $MOTIF/CTCF_Q1-ClosestDyad.bed -o $MOTIF/500bp/CTCF_Q1-ClosestDyad_500bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 500 $MOTIF/CTCF_Q4.bed -o $MOTIF/500bp/CTCF_Q4_500bp.bed

java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $TEMP/CTCF_Q1.bed -o $MOTIF/1000bp/CTCF_Q1_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $TEMP/CTCF_Q2.bed -o $MOTIF/1000bp/CTCF_Q2_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $TEMP/CTCF_Q3.bed -o $MOTIF/1000bp/CTCF_Q3_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $TEMP/CTCF_Q4.bed -o $MOTIF/1000bp/CTCF_Q4_1000bp.bed

java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1 $TEMP/CTCF_Q1.bed -o $MOTIF/CTCF_Q1.bed

## take ramdom site 
bedtools random -l 1 -n 3000 -g $GINFO > $TEMP/random_3000.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1 $TEMP/random_3000.bed -o $MOTIF/Ramdom_Q1.bed
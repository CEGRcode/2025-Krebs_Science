#!/bin/bash

# Format ZKSCAN1 reference point file to shift strands and flip orientation so that it better matches JASPAR motifs.

# data/RefPT-Motif
#   |--ZKSCAN1_SORT-Occupancy.bed
#   |--1000bp
#     |--ZKSCAN1_SORT-Occupancy_1000bp.bed              (ZKSCAN1_Occupancy_flip_shift2_1000bp.bed)

### CHANGE ME
WRK=/path/to/2024-Chen_Nature/02_Call_RefPT
WRK=/storage/home/owl5022/scratch/2024-Chen_Nature/02_Call_RefPT
###

# Dependencies
# - bedtools
# - java

set -exo
module load bedtools
module load anaconda
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Inputs and outputs
MOTIF=../data/RefPT-Motif/
GENOME=../data/hg38_files/hg38.fa
GINFO=../data/hg38_files/hg38.chrom.sizes.txt
ZKSCAN1_BOUND=temp-3_Filter_and_Sort_by_occupancy/ZKSCAN1_ZKSCAN1-K562_M1_100bp_7-Occupancy_BOUND_1bp.bed

[ -d $MOTIF ] || mkdir $MOTIF
[ -d $MOTIF/1000bp ] || mkdir $MOTIF/1000bp

# Script shortcuts
SCRIPTMANAGER=../bin/ScriptManager-v0.15.jar

TEMP=temp-7_ZKSCAN1
[ -d $TEMP ] || mkdir $TEMP

# Flip strands (reorient and shift to match JASPAR's)
awk 'BEGIN{OFS="\t";FS="\t"}{
    if ($6 == "-") {
        $6="+";
    } else {
        $6="-";
    }
    print $1,$2,$3,$4,$5,$6,$7;
}' $ZKSCAN1_BOUND > $TEMP/ZKSCAN1_BOUND_flipped.bed

# Shift 2bp downstream
bedtools shift -p 2 -m -2 -g $GINFO -i $TEMP/ZKSCAN1_BOUND_flipped.bed > $TEMP/ZKSCAN1_BOUND_shifted.bed

# Deduplicate (unique)
sort -uk1,3 $TEMP/ZKSCAN1_BOUND_shifted.bed | sort -rnk7,7 | cut -f1-6 > $MOTIF/ZKSCAN1_SORT-Occupancy.bed

# Expand 1000bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/ZKSCAN1_SORT-Occupancy.bed  -o $MOTIF/1000bp/ZKSCAN1_SORT-Occupancy_1000bp.bed

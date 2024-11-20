#!/bin/bash

# Format CTCF reference point file to proper BED files.

# data/RefPT-Motif
#   |--CTCF_SORT-Occupancy.bed
#   |--1000bp
#     |--CTCF_SORT-Occupancy_1000bp.bed

### CHANGE ME
WRK=/path/to/2024-Krebs_Science/04_Call_Motifs
WRK=/scratch/owl5022/2024-Krebs_Science/04_Call_Motifs
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

# Script shortcuts
SCRIPTMANAGER=../bin/ScriptManager-v0.15.jar

# Reformat BED with SCORE=tag pileup tag count (move 7th column to score column and add K562_CTCF in the seventh)
cut -f1-6 $CTCF_BOUND > $MOTIF/CTCF_SORT-Occupancy.bed

# QC: Stat lines
wc -l $MOTIF/CTCF_SORT-Occupancy.bed

# Expand 1000bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/CTCF_SORT-Occupancy.bed -o $MOTIF/1000bp/CTCF_SORT-Occupancy_1000bp.bed

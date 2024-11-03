#!/bin/bash

### CHANGE ME
WRK=/path/to/2024-Chen_Nature/02_Call_RefPT
WRK=/storage/home/owl5022/scratch/2024-Chen_Nature/02_Call_RefPT
###

# Dependencies
# - bedtools
# - java

set -exo
module load anaconda
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Inputs and outputs
MOTIF=$WRK/../data/RefPT-Motif

# Script shortcuts
SCRIPTMANAGER=../bin/ScriptManager-v0.15.jar

TEMP=temp-6_CTCF_motif
[ -d $TEMP ] || mkdir $TEMP
[ -d Motif_analysis/CTCF ] || mkdir Motif_analysis/CTCF

cd $WRK/Motif_analysis/CTCF

awk '{OFS="\t"} {print $1,$2,$3,$4,$7,$6,"K562_CTCF"}' FIMO/CTCF/CTCF_Occupancy_1bp.bed  | awk '
{
    if ($1 !~/alt/ && $1 !~/random/ && $1 !~/Un/ ) {
        print $0 > "CTCF_Occuoancy_1bp.bed";
    } 
}'

wc -l CTCF_Occupancy_1bp.bed


java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 CTCF_Occuoancy.bed -o $MOTIF/1000bp/CTCF_Occuoancy_1000bp.bed

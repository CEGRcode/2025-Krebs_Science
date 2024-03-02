#!/bin/bash

# Sort each FIMO'd scrambled motif by FIMO score and take the top 10K for each

### CHANGE ME
WRK=/path/to/2023-Chen_Benzonase-ChIP-exo/02_Call_RefPT
###

# Dependencies
# - java

set -exo

# Script shortcuts
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.14.jar

# Inputs and outputs
BLACKLIST=$WRK/../data/hg19_files/hg19_exclude.bed
MOTIF=$WRK/../data/RefPT-Motif

[ -d $MOTIF ] || mkdir $MOTIF

# Loop through each scrambled* motif
for PWM in PWM/scrambled*.meme.txt;
do
    # Parse TF from PWM filename
    TF=$(basename "$PWM" .meme.txt)
    # Sort by FIMO score and get top 10K
    awk '{print $1"\t"$2"\t"$3"\t"$1"_"$2"_"$3"\t"$6"\t"$7}' FIMO/$TF/$TF\_motif1_unsorted.bed | sort -k5,5nr | head -n 10000 >  $MOTIF/$TF\_top10k.bed
    # Expand by 1bp and 1000bp
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1 $MOTIF/$TF\_top10k.bed -o $MOTIF/1bp/$TF\_top10k_1bp.bed
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/$TF\_top10k.bed -o $MOTIF/1000bp/$TF\_top10k_1000bp.bed
done
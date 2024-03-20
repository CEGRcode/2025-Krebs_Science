#!/bin/bash

# Plot figure of WDR5 motif relative to TSS, orientation seperately


### CHANGE ME
WRK=/path/to/2024-Chen_Nature/Z_Figures
WRK=/storage/home/owl5022/scratch/2024-Chen_Nature/Z_Figures
###

# Dependencies
# - java
# - samtools

set -exo
module load samtools

# Fill in placeholder constants with your directories
BAMDIR=$WRK/../data/BAM
MOTIF=$WRK/../data/RefPT-Motif
OTHER=$WRK/../data/RefPT-Other
TSS=$WRK/../data/RefPT-Other/TSS_DistWDR5.bed

# Setup ScriptManager for job array
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.14.jar

# Set up output directories
[ -d logs ] || mkdir logs
[ -d F4/b ] || mkdir -p F4/b

# find WDR5 orientation
awk '{OFS="\t"}{if ($6==$12) print }' $TSS | cut -f7-13 > F4/b/WDR5_TSS-strand-same.bed
awk '{OFS="\t"}{if ($6!=$12) print }' $TSS | cut -f7-13 > F4/b/WDR5_TSS-strand-opposite.bed

# Re-expand WDR5 motif to make a stronger mark in peak-align heatmap
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 32 F4/b/WDR5_TSS-strand-same.bed     -o F4/b/WDR5_TSS-strand-same_32bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 32 F4/b/WDR5_TSS-strand-opposite.bed -o F4/b/WDR5_TSS-strand-opposite_32bp.bed

# Peak-align each strand match and mismatch
java -jar $SCRIPTMANAGER peak-analysis peak-align-ref F4/b/WDR5_TSS-strand-same_32bp.bed     $OTHER/TSS_DistWDR5_1000bp.bed -o F4/b/WDR5-same_TSS_DistWDR5_1000bp.cdt
java -jar $SCRIPTMANAGER peak-analysis peak-align-ref F4/b/WDR5_TSS-strand-opposite_32bp.bed $OTHER/TSS_DistWDR5_1000bp.bed -o F4/b/WDR5-opposite_TSS_DistWDR5_1000bp.cdt

# Heatmap match and mismatch (different colors)
java -jar $SCRIPTMANAGER figure-generation heatmap -a 1 --blue F4/b/WDR5-same_TSS_DistWDR5_1000bp.cdt     -o SAME.png
java -jar $SCRIPTMANAGER figure-generation heatmap -a 1 --red  F4/b/WDR5-opposite_TSS_DistWDR5_1000bp.cdt -o OPPOSITE.png

NSITES=`wc -l $OTHER/1000bp/TSS_DistWDR5_1000bp.bed |awk '{print $1}'`

# Merge and label
java -jar $SCRIPTMANAGER figure-generation merge-heatmap SAME.png OPPOSITE.png -o F4/b/WDR5-StrandedPeakAlign_TSS_DistWDR5_1000bp.png
java -jar $SCRIPTMANAGER figure-generation label-heatmap F4/b/WDR5-StrandedPeakAlign_TSS_DistWDR5_1000bp.png 
	-l "-500" -m "0" -r "+500" -x $BED -w 1 -f 20 \
	-x "Distance from TSS" -y "TSS sorted by distance to nearby WDR5 motif (${NSITES} sites)" \
	-o F4/b/WDR5-StrandedPeakAlign_TSS_DistWDR5_1000bp.svg

# Clean-up
rm SAME.png OPPOSITE.png

## Add barplot scripts
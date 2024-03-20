#!/bin/bash

# Brief shell script to make F2e (heatmap showing sense-strand exonuclease cut patterns)

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

# Setup ScriptManager for job array
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.14.jar

# Set up output directories
[ -d logs ] || mkdir logs
[ -d F2/e ] || mkdir -p F2/e

# Count sites
NSITES=`wc -l $MOTIF/1000bp/NFIA_Occupancy_1000bp.bed | awk '{print $1-1}'`

# Pileup for fat sense strand figure
BASE=NFIA-NFIA_Occupancy_500bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 500 $MOTIF/1000bp/NFIA_Occupancy_1000bp.bed -o NFIA_Occupancy_500bp.bed
java -jar $SCRIPTMANAGER read-analysis tag-pileup NFIA_Occupancy_500bp.bed $BAMDIR/K562_NFIA_BX_rep1_hg19.bam -M $BASE
java -jar $SCRIPTMANAGER figure-generation heatmap $BASE\_sense.cdt \
    -p 0.95 -x 400 -y 600 --blue -o $BASE\_sense_treeview.png
java -jar $SCRIPTMANAGER figure-generation label-heatmap $BASE\_sense_treeview.png \
    -l "-500" -m "0" -r "+500" -w 1 -f 20 \
    -x "NFIA_Occupancy_1000bp" -y "NFIA_Occupancy_1000bp (${NSITES} sites)" \
    -o F2/e/K562_NFIA_BX_rep1_hg19_NFIA_Occupancy_1000bp_read1_sense_treeview_label.svg

# Clean-up
rm NFIA_Occupancy_500bp.bed $BASE\_*.cdt $BASE\_sense_treeview.png
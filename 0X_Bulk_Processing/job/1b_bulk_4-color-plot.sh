#!/bin/bash

# Re-expand each BED file in LIST to 32bp and generate a 4-color plot
# Note: execute from /path/to/2023-Chen_Benzonase-ChIP-exo/0X_Bulk_Processing

# Dependencies
# - java

# Define the source files as LIST
LIST=(
    "CTCF_Occupancy_1000bp.bed"
    "NFIA_Occupancy_1000bp.bed"
    "GABPA_Occupancy_1000bp.bed"
    "SP1_Occupancy_1000bp.bed"
    "WDR5_Occupancy_1000bp.bed"
    "SRF_Occupancy_1000bp.bed"
    "NRF1_Occupancy_1000bp.bed"
    "MEIS2_Occupancy_1000bp.bed"
    "USF1_Occupancy_1000bp.bed"
    "USF2_Occupancy_1000bp.bed"
    "YY1_Occupancy_1000bp.bed"
    "E4F1_Occupancy_1000bp.bed"
    "scrambledCTCF_NucSort_1000bp.bed"
    "scrambledNFIA_NucSort_1000bp.bed"
)

# Script shortcuts
SCRIPTMANAGER=../bin/ScriptManager-v0.14.jar

# Inputs and outputs
GENOME=$WRK/../data/hg19_files/hg19.fa
MOTIFS=../data/RefPT-Motif
ODIR=Library-FourColor

[ -d $ODIR ] || mkdir $ODIR

# Loop through the source files and copy them to the destination directory
for FILENAME in "${LIST[@]}"; do
    BEDFILE=$MOTIFS/1000bp/$FILENAME
    [ -f $BEDFILE ] || echo "Missing $BEDFILE"

    # Parse base from BED filename
    BED=`basename $BEDFILE "_1000bp.bed"`

    # Expand, extract, and generate four-color plot
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 32 $BEDFILE -o $ODIR/$BED\_32bp.bed
    java -jar $SCRIPTMANAGER sequence-analysis fasta-extract $ODIR/$BED\_32bp.bed -o $ODIR/$BED\_32bp.fa
    java -jar $SCRIPTMANAGER figure-generation four-color $ODIR/$BED\_32bp.fa -o $ODIR/$BED\_32bp.png
    # Add SVG label
    java -jar $SCRIPTMANAGER figure-generation label-heatmap $ODIR/$BED\_32bp.png -o $ODIR/$BED\_32bp.svg \
        -x $BED -l "-15" -m "0" -r "+15" \
        -y "$BED Motif occurences (____ sites)"
done

#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 00:40:00
#SBATCH -A open
#SBATCH -o logs/1b_bulk_4-color-plot.log.out
#SBATCH -e logs/1b_bulk_4-color-plot.log.err

# Re-expand each BED file in LIST to 32bp and generate a 4-color plot

### CHANGE ME
WRK=/path/to/2024-Chen_Nature/0X_Bulk_Processing
WRK=/storage/home/owl5022/scratch/2024-Chen_Nature/0X_Bulk_Processing
###

# Dependencies
# - java
# - opencv
# - python

set -exo
module load anaconda
conda activate bx

# Define the source files as LIST
LIST=(
    "$WRK/../data/RefPT-Motif/1000bp/CTCF_Occupancy_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/NFIA_Occupancy_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/GABPA_Occupancy_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/SP1_Occupancy_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/WDR5_Occupancy_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/CTCF_NucSort_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/NFIA_NucSort_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/GABPA_NucSort_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/SP1_NucSort_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/WDR5_NucSort_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/scrambledCTCF_top10k_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/scrambledNFIA_top10k_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/scrambledWDR5_top10k_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/scrambledCTCF_NucSort_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/scrambledNFIA_NucSort_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/scrambledWDR5_NucSort_1000bp.bed"
)

# Script shortcuts
SCRIPTMANAGER=../bin/ScriptManager-v0.14.jar
RESIZE=../bin/resize_png.py

# Inputs and outputs
GENOME=$WRK/../data/hg19_files/hg19.fa
MOTIFS=../data/RefPT-Motif
OUTDIR=Library

[ -d logs ] || mkdir logs
[ -d $OUTDIR ] || mkdir $OUTDIR

# Loop through the source files and copy them to the destination directory
for BEDFILE in "${LIST[@]}";
do
    # Parse base from BED filename
    BED=`basename $BEDFILE ".bed"`

    DIR=$OUTDIR/$BED
    [ -d $DIR ] || mkdir $DIR
    [ -d $DIR/FourColor ] || mkdir $DIR/FourColor

    # Expand, extract, and generate four-color plot
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 32 $BEDFILE -o $DIR/FourColor/$BED\_32bp.bed
    java -jar $SCRIPTMANAGER sequence-analysis fasta-extract $GENOME $DIR/FourColor/$BED\_32bp.bed -o $DIR/FourColor/$BED\_32bp.fa
    java -jar $SCRIPTMANAGER figure-generation four-color $DIR/FourColor/$BED\_32bp.fa -o $DIR/FourColor/$BED\_32bp.png
    python $RESIZE -r 600 -c 200 -i $DIR/FourColor/$BED\_32bp.png -o $DIR/FourColor/$BED\_32bp_RESIZE.png

    # Count sites
    NSITES=`wc -l $BEDFILE | awk '{print $1}'`

    # Add SVG label
    java -jar $SCRIPTMANAGER figure-generation label-heatmap $DIR/FourColor/$BED\_32bp_RESIZE.png -o $DIR/FourColor/$BED\_32bp.svg \
        -x $BED -l "-15" -m "0" -r "+15" \
        -y "$BED Motif occurences (${NSITES} sites)"
done

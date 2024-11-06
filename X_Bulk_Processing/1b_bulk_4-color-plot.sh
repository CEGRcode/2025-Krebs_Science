#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 00:40:00
#SBATCH -A open
#SBATCH -o logs/1b_bulk_4-color-plot.log.out
#SBATCH -e logs/1b_bulk_4-color-plot.log.err

# Re-expand each BED file in LIST to 32bp and generate a 4-color plot

### CHANGE ME
WRK=/path/to/2024-Chen_Nature/X_Bulk_Processing
WRK=/storage/home/owl5022/scratch/2024-Chen_Nature/X_Bulk_Processing
###

# Dependencies
# - java
# - opencv
# - python

set -exo
module load anaconda
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Define the source files as LIST
LIST=(
    "$WRK/../data/RefPT-Motif/1000bp/CTCF_SORT-Occupancy_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/FOXA_SORT-ClosestDyad_STACK-K562-uHepG2_1000bp.bed"
    "$WRK/../data/RefPT-Motif/500bp/NFIA_SORT-Occupancy_500bp.bed"
)

# Script shortcuts
SCRIPTMANAGER=../bin/ScriptManager-v0.15.jar
RESIZE=../bin/resize_png.py

# Inputs and outputs
GENOME=../data/hg38_files/hg38.fa
MOTIFS=../data/RefPT-Motif

[ -d logs ] || mkdir logs
[ -d Library ] || mkdir Library

# Loop through the source files and copy them to the destination directory
for BEDFILE in "${LIST[@]}";
do
    # Parse base from BED filename
    BED=`basename $BEDFILE ".bed"`

    DIR=Library/$BED
    [ -d $DIR ] || mkdir $DIR
    [ -d $DIR/FourColor ] || mkdir $DIR/FourColor

    # Expand 32bp
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 32 $BEDFILE -o $DIR/FourColor/${BED}_32bp.bed
    # Take off 1bp (since all motifs are odd-lengthed, we can take 1bp off to center 0bp label on motif center)
    awk '{OFS="\t"}{if($6=="-") print $1,$2+1,$3,$4,$5,$6; else print $1,$2,$3-1,$4,$5,$6}' $DIR/FourColor/${BED}_32bp.bed > $DIR/FourColor/${BED}_31bp.bed
    # Extract sequence
    java -jar $SCRIPTMANAGER sequence-analysis fasta-extract $GENOME $DIR/FourColor/${BED}_31bp.bed -o $DIR/FourColor/${BED}_31bp.fa
    # Generate four-color plot from sequence
    java -jar $SCRIPTMANAGER figure-generation four-color $DIR/FourColor/${BED}_31bp.fa -o $DIR/FourColor/${BED}_31bp.png
    # Resize PNG to standard px dimensions
    python $RESIZE -r 600 -c 200 -i $DIR/FourColor/${BED}_31bp.png -o $DIR/FourColor/${BED}_31bp_RESIZE.png
    # Count sites
   # NSITES=`wc -l $BEDFILE | awk '{print $1}'`

    # Add SVG label
    java -jar $SCRIPTMANAGER figure-generation label-heatmap $DIR/FourColor/${BED}_31bp_RESIZE.png \
        -l "-15" -m "0" -r "+15" -w 1 -f 20         -o $DIR/FourColor/${BED}_31bp.svg
done




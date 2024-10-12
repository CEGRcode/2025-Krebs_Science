#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 00:40:00
#SBATCH -A open
#SBATCH -o logs/1b_bulk_4-color-plot.log.out
#SBATCH -e logs/1b_bulk_4-color-plot.log.err

# Re-expand each BED file in LIST to 32bp and generate a 4-color plot

### CHANGE ME
WRK=/Path/to/Title
###

# Dependencies
# - java
# - opencv
# - python

set -exo
module load anaconda

# Define the source files as LIST
LIST=(
    "$WRK/data/RefPT-Motif/1000bp/CTCF_Occupancy_1000bp.bed"
    "$WRK/data/RefPT-Motif/1000bp/FOXA_all_1000bp.bed"
    "$WRK/data/RefPT-Motif/1000bp/NFIA_downNuc_500bp.bed"
)

# Script shortcuts
SCRIPTMANAGER=$WRK/bin/ScriptManager-v0.15.jar
RESIZE=$WRK/bin/resize_png.py

# Inputs and outputs
GENOME=$WRK/data/hg38_files/hg38.fa
MOTIFS=$WRK/data/RefPT-Motif
OUTDIR=$WRK/Library

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

    # Expand 32bp
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 32 $BEDFILE -o $DIR/FourColor/$BED\_32bp.bed
    # Take off 1bp (since all motifs are odd-lengthed, we can take 1bp off to center 0bp label on motif center)
    awk '{OFS="\t"}{if($6=="-") print $1,$2+1,$3,$4,$5,$6; else print $1,$2,$3-1,$4,$5,$6}' $DIR/FourColor/$BED\_32bp.bed > $DIR/FourColor/$BED\_31bp.bed
    # Extract sequence
    java -jar $SCRIPTMANAGER sequence-analysis fasta-extract $GENOME $DIR/FourColor/$BED\_31bp.bed -o $DIR/FourColor/$BED\_31bp.fa
    # Generate four-color plot from sequence
    java -jar $SCRIPTMANAGER figure-generation four-color $DIR/FourColor/$BED\_31bp.fa -o $DIR/FourColor/$BED\_31bp.png
    # Resize PNG to standard px dimensions
    sips -z 600 200 $DIR/FourColor/$BED\_31bp.png --out $DIR/FourColor/$BED\_31bp_RESIZE.png
    # Count sites
   # NSITES=`wc -l $BEDFILE | awk '{print $1}'`

    # Add SVG label
    java -jar $SCRIPTMANAGER figure-generation label-heatmap $DIR/FourColor/$BED\_31bp_RESIZE.png \
        -l "-15" -m "0" -r "+15" -w 1 -f 20         -o $DIR/FourColor/$BED\_31bp.svg
done




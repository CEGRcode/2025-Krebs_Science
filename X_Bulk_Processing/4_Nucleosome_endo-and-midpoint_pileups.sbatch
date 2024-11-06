#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 00:15:00
#SBATCH -A open
#SBATCH -o logs/3_Nucleosome_endo-and-midpoint_pileups.log.out-%a
#SBATCH -e logs/3_Nucleosome_endo-and-midpoint_pileups.log.err-%a
#SBATCH --array 1-8

# Pileup midpoint (heatmaps) and composites (endo cut sites) for a custom
# combination of BAM x BED files (primarily nucleosome BAM). No scaling.
#  - data/BAM/K562_-_BI_rep1_hg19.bam

### CHANGE ME
WRK=/path/to/2024-Chen_Nature/X_Bulk_Processing
WRK=/storage/home/owl5022/scratch/2024-Chen_Nature/X_Bulk_Processing
METADATA=Nucleosome_pileups.txt
###

# Dependencies
# - java
# - samtools

set -exo
module load samtools

# Fill in placeholder constants with your directories
BAMDIR=../data/BAM

# Setup ScriptManager for job array
ORIGINAL_SCRIPTMANAGER=../bin/ScriptManager-v0.15.jar
SCRIPTMANAGER=../bin/ScriptManager-v0.15-$SLURM_ARRAY_TASK_ID.jar
cp $ORIGINAL_SCRIPTMANAGER $SCRIPTMANAGER

# Script shortcuts
COMPOSITE=../bin/sum_Col_CDT.pl

# Set up output directories
[ -d logs ] || mkdir logs
[ -d Library ] || mkdir Library

# Determine BAM file for the current job array index
BAMFILE=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $1}'`
BAM=`basename $BAMFILE ".bam"`
[ -f $BAMFILE.bai ] || samtools index $BAMFILE

# Determine BED file for the current job array index
BEDFILE=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $2}'`
BED=`basename $BEDFILE ".bed"`

echo "(${SLURM_ARRAY_TASK_ID}) ${BEDFILE} "x" ${BAMFILE} "

DIR=Library/$BED
[ -d $DIR ] || mkdir $DIR
[[ -d $DIR/CDT ]] || mkdir $DIR/CDT
[[ -d $DIR/Composites ]] || mkdir $DIR/Composites
[[ -d $DIR/PNG/Strand ]] || mkdir -p $DIR/PNG/Strand
[[ -d $DIR/SVG ]] || mkdir $DIR/SVG

echo "Run Nucleosome mid pileup"
BASE=$BAM\_$BED\_midpoint

# Pileup (mid)
java -jar $SCRIPTMANAGER read-analysis tag-pileup $BEDFILE $BAMFILE --cpu 4 -m --combined -o $DIR/Composites/$BASE.out -M $DIR/CDT/$BASE
# No scaling

# Two-color heatmap
java -jar $SCRIPTMANAGER figure-generation heatmap -p .95 --black $DIR/CDT/$BASE\_combined.cdt -o $DIR/PNG/$BASE\_combined_treeview.png

# Count sites
NSITES=`wc -l $BEDFILE | awk '{print $1-1}'`

# Label heatmap based on BED naming
if [[ "$BED" == *"1000bp"* ]]; then  
    java -jar "$SCRIPTMANAGER" figure-generation label-heatmap "$DIR/PNG/${BASE}_combined_treeview.png" \
        -l "-500" -m "0" -r "+500" -w 1 -f 20 \
        -o "$DIR/SVG/${BASE}_combined_treeview_label.svg"

elif [[ "$BED" == *"500bp"* ]]; then  
    java -jar "$SCRIPTMANAGER" figure-generation label-heatmap "$DIR/PNG/${BASE}_combined_treeview.png" \
        -l "-250" -m "0" -r "+250" -w 1 -f 20 \
        -o "$DIR/SVG/${BASE}_combined_treeview_label.svg"

elif [[ "$BED" == *"250bp"* ]]; then  
    java -jar "$SCRIPTMANAGER" figure-generation label-heatmap "$DIR/PNG/${BASE}_combined_treeview.png" \
        -l "0" -m "125" -r "+250" -w 1 -f 20 \
        -o "$DIR/SVG/${BASE}_combined_treeview_label.svg"
fi



echo "Run Nucleosome endo cut sites (read1 and read2) pileup"
BASE_read2=$BAM\_$BED\_read2
BASE_read1=$BAM\_$BED\_read1

# Pileup (endo)
java -jar $SCRIPTMANAGER read-analysis tag-pileup $BEDFILE $BAMFILE --cpu 4 -5 -1 -o $DIR/Composites/$BASE_read1.out
java -jar $SCRIPTMANAGER read-analysis tag-pileup $BEDFILE $BAMFILE --cpu 4 -5 -2 -o $DIR/Composites/$BASE_read2.out

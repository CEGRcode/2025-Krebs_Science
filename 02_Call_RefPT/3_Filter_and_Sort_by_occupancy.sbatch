#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 1:00:00
#SBATCH -A open
#SBATCH -o logs/3_Filter_and_Sort_by_occupancy.log.out-%a
#SBATCH -e logs/3_Filter_and_Sort_by_occupancy.log.err-%a
#SBATCH --array 1-8

# Sort filtered Motif BED Files by the first BX rep1 BAM file occupancy for each TF target (skip scrambled)

### CHANGE ME
WINDOW=100
WRK=/path/to/2024-Chen_Nature/02_Call_RefPT
###

# Dependencies
# - java
# - python
# - samtools

set -exo
module load anaconda3
module load bedtools
module load samtools
source activate bx

# Script shortcuts
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.14.jar

# Inputs and outputs
BLACKLIST=$WRK/../data/hg19_files/hg19_exclude.bed
MOTIF=$WRK/../data/RefPT-Motif

[ -d $MOTIF ] || mkdir $MOTIF

# Parse TF from PWM filename
PWM=`ls PWM/*.meme.txt | head -n $SLURM_ARRAY_TASK_ID | tail -1`;
TF=`basename $PWM ".meme.txt"`

# Skip if scrambled
[[ $TF == "scrambled"* ]] && exit

# Get BAM file
BAMFILE=$WRK/../data/BAM/K562_$TF\_BX_rep1_hg19.bam
[ -f $BAMFILE.bai ] || samtools index $BAMFILE
BAM=`basename $BAMFILE ".bam"`

## Sort Motifs by occupancy

# Expand WINDOW bp around motif
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c $WINDOW FIMO/$TF/$TF\_motif1_unsorted.bed -o FIMO/$TF/$TF\_MOTIF1_unsorted_$WINDOW\bp.bed

# Tag Pileup for sort...shift 8bp XL dist
java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -1 --combined FIMO/$TF/$TF\_MOTIF1_unsorted_$WINDOW\bp.bed $BAMFILE -M FIMO/$TF/$TF\_MOTIF1_unsorted_$WINDOW\bp_$BAM

# Sort BED
java -jar $SCRIPTMANAGER coordinate-manipulation sort-bed -c $WINDOW FIMO/$TF/$TF\_MOTIF1_unsorted_$WINDOW\bp.bed FIMO/$TF/$TF\_MOTIF1_unsorted_$WINDOW\bp_$BAM\_combined.cdt -o FIMO/$TF/$TF\_MOTIF1_sorted

# Heatmap of Sort
java -jar $SCRIPTMANAGER figure-generation heatmap --black -r 1 -l 2 -p .95 FIMO/$TF/$TF\_MOTIF1_sorted.cdt -o FIMO/$TF/$TF\_MOTIF1_sorted.png

# Clean-up
rm FIMO/$TF/$TF\_MOTIF1_unsorted_$WINDOW\bp_$BAM\_combined.cdt FIMO/$TF/$TF\_MOTIF1_unsorted_$WINDOW\bp.bed


## Calculate occupancy and sort by inflection point

[ -d FIMO/$TF/InflectionFilter ] || mkdir -p FIMO/$TF/InflectionFilter

# Sum occupancy
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum FIMO/$TF/$TF\_MOTIF1_sorted.cdt -o FIMO/$TF/InflectionFilter/

# Calculate and plot inflection point
python scripts/get_inflection_point.py -i FIMO/$TF/InflectionFilter/$TF\_MOTIF1_sorted_SCORES.out -o FIMO/$TF/InflectionFilter/$TF\_MOTIF1_sorted_SCORES_inflectionPT.png > FIMO/$TF/InflectionFilter/$TF\_MOTIF1_sorted_Threshold.txt

#ls FIMO/$TF/InflectionFilter/$TF\_MOTIF1_sorted_Threshold.txt
#THRESH=`grep "Second Derivative: " FIMO/$TF/InflectionFilter/$TF\_MOTIF1_sorted_Threshold.txt | awk '{print $3}'`
THRESH=`grep "KneeLocator: " FIMO/$TF/InflectionFilter/$TF\_MOTIF1_sorted_Threshold.txt | awk '{print $2}'`
echo $THRESH
awk -v THRESH=$THRESH '{FS="\t"}{OFS="\t"}{if($2>=THRESH) print}' FIMO/$TF/InflectionFilter/$TF\_MOTIF1_sorted_SCORES.out > FIMO/$TF/InflectionFilter/$TF\_MOTIF1_sorted_SCORES_Filter-$THRESH.out

# Pull top sites and save to ../data/Motif directory
NLINES=`wc -l FIMO/$TF/InflectionFilter/$TF\_MOTIF1_sorted_SCORES_Filter-$THRESH.out | awk '{print $1}'`
head -n $NLINES FIMO/$TF/$TF\_MOTIF1_sorted.bed >  $MOTIF/$TF\_Occupancy.bed

[ -d $MOTIF/1000bp ] || mkdir $MOTIF/1000bp
[ -d $MOTIF/1bp ] || mkdir $MOTIF/1bp

# Expand to 1000bp window and 1 bp Window
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/$TF\_Occupancy.bed -o $MOTIF/1000bp/$TF\_Occupancy_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1 $MOTIF/$TF\_Occupancy.bed -o $MOTIF/1bp/$TF\_Occupancy_1bp.bed
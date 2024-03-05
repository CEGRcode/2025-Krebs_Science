#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 00:10:00
#SBATCH -A open
#SBATCH -o logs/7_Motif_Sort_by_Nucleosome.log.out-%a
#SBATCH -e logs/7_Motif_Sort_by_Nucleosome.log.err-%a
#SBATCH --array 1-8

# Sort every TF by distance to the nearest nucleosome

### CHANGE ME
WRK=/path/to/2024-Chen_Nature/02_Call_RefPT
###

# Dependencies
# - bedtools (v2.31.0)
# - java

set -exo
module load bedtools

# Script shortcuts
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.14.jar

# Inputs and outputs
BLACKLIST=$WRK/../data/hg19_files/hg19_exclude.bed
MOTIF=$WRK/../data/RefPT-Motif
NUCLEOSOME=$WRK/../data/RefPT-Other/1bp/hg19_K562_Nucleosome_1bp_SORT-genomic_nostrand.bed

[ -d $MOTIF ] || mkdir $MOTIF

# Parse TF from PWM filename
PWM=`ls PWM/*.meme.txt | head -n $SLURM_ARRAY_TASK_ID | tail -1`;
TF=`basename $PWM ".meme.txt"`
BEDFILE=$MOTIF/1bp/$TF\_Occupancy_1bp.bed

# Change suffix name if scrambled
[[ $TF == "scrambled"* ]] && BEDFILE=$MOTIF/1bp/$TF\_top10k_1bp.bed

# Sort the motif BED file by genomic coord (for bedtools closest)
bedtools sort -i $BEDFILE > $TF\_genomic_sort.bed
# Report distance, wrt A's strand, use first occurance in case of tie
bedtools closest -d -D a -t first -a $TF\_genomic_sort.bed -b $NUCLEOSOME > $TF\_closest.bed
# Reformat BAM with genomic coord as ID and sort on distance
awk '{FS="\t"}{OFS="\t"}{print $1,$2,$3,$1"_"$2"_"$3,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}' $TF\_closest.bed \
  | sort -k13,13n > $MOTIF/1bp/$TF\_NucSort_1bp.bed                                                      # sort by distance and write

# Clean-up
rm $TF\_genomic_sort.bed $TF\_closest.bed

# | sort -k4,4 | awk '{FS="\t"}{OFS="\t"}{print $1,$2,$3,$1"_"$2"_"$3,$6,$7,$8,$9,$10,$11,$12,$13,$14}' \ # reformat BAM with genomic coord as ID

[ -d $MOTIF/1000bp ] || mkdir $MOTIF/1000bp

# Expand the BED files by 1000bp and save in the 1000bp directory
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/1bp/$TF\_NucSort_1bp.bed -o $MOTIF/1000bp/$TF\_NucSort_1000bp.bed

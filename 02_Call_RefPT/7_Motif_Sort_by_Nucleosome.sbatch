#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 00:10:00
#SBATCH -A open
#SBATCH -o logs/7_Motif_Sort_by_Nucleosome.log.out
#SBATCH -e logs/7_Motif_Sort_by_Nucleosome.log.err
#SBATCH --array 1-NSAMPLES

SLURM_ARRAY_TASK_ID=5

# Sort every TF by distance to the nearest nucleosome

### CHANGE ME
WRK=/path/to/2023-Chen_Benzonase-ChIP-exo/02_Call_RefPT
WRK=`pwd`
###

# Dependencies
# - bedtools

set -exo
# module load anaconda3
source activate bx

# Script shortcuts
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.14.jar

# Inputs and outputs
BLACKLIST=$WRK/../data/hg19_files/hg19_exclude.bed
MOTIF=$WRK/../data/RefPT-Motif
NUCLEOSOME=$WRK/../data/RefPT-Other/1bp/hg19_K562_Nucleosome_1bp_SORT-genomic.bed

[ -d $MOTIF ] || mkdir $MOTIF

# Parse TF from PWM filename
PWM=`ls PWM/*.meme.txt | head -n $SLURM_ARRAY_TASK_ID | tail -1`;
TF=$(basename "$PWM" .meme.txt)
BEDFILE=$MOTIF/1bp/$TF\_Occupancy_1bp.bed

# Skip if scrambled
[[ $TF == "scrambled"* ]] && BEDFILE=$MOTIF/1bp/$TF\_top10K_1bp.bed

# Sort and process the BED file and categorize results
bedtools sort -i $BEDFILE \
  | bedtools closest -a - -b $NUCLEOSOME -d -D a \
  | awk '{FS="\t"}{OFS="\t"}{print $1,$2,$3,$1"_"$2"_"$3,$6,$7,$8,$9,$10,$11,$12,$13,$14}' \
  | sort -k13,13n > $MOTIF/1bp/$TF\_NucSort_1bp.bed

[ -d $MOTIF/1000bp ] || mkdir $MOTIF/1000bp

# Expand the BED files by 1000bp and save in the 1000bp directory
java -jar "$SCRIPTMANAGER" coordinate-manipulation expand-bed -c 1000 $MOTIF/1bp/$TF\_NucSort_1bp.bed -o $MOTIF/1000bp/$TF\_NucSort_1000bp.bed

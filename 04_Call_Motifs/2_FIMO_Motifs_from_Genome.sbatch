#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=24gb
#SBATCH -t 5:00:00
#SBATCH -A open
#SBATCH -o logs/2_FIMO_Motifs_from_Genome.log.out-%a
#SBATCH -e logs/2_FIMO_Motifs_from_Genome.log.err-%a
#SBATCH --array 1-4

# FIMO the reference genome for each motif in the PWM directory

### CHANGE ME
EXCLUSION=100
WRK=/path/to/2024-Krebs_Science/04_Call_Motifs
WRK=/scratch/owl5022/2024-Krebs_Science/04_Call_Motifs
###

# Dependencies
# - bedtools
# - java
# - MEME suite (FIMO)

set -exo
module load bedtools
module load anaconda3
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Inputs and outputs
GENOME=../data/hg38_files/hg38.fa
BLACKLIST=../data/hg38_files/hg38-blacklist.bed

# Script shortcuts
ORIGINAL_SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15-$SLURM_ARRAY_TASK_ID.jar
cp $ORIGINAL_SCRIPTMANAGER $SCRIPTMANAGER

PWM_LIST=(
    PWM/CTCF_M1.meme.txt
    PWM/FOXA2_M1.meme.txt
    PWM/NFIA_M1.meme.txt
    PWM/ZKSCAN1_M1.meme.txt
)

# Define PWM file path based on SLURM_ARRAY_TASK_ID index
INDEX=$(($SLURM_ARRAY_TASK_ID-1))
PWM=${PWM_LIST[$INDEX]}

# Parse TF from PWM filename
TF=`basename $PWM "_M1.meme.txt"`

# Create output directories if they don't exist
[ -d logs ] || mkdir logs
[ -d FIMO ] || mkdir FIMO
[ -d FIMO/$TF ] || mkdir FIMO/$TF

echo "($SLURM_ARRAY_TASK_ID) $TF"

# Run FIMO to scan genome for motif occurrences (keep strand with max score)
fimo --verbosity 1 --thresh 1.0E-4 --max-strand --oc FIMO/$TF $PWM $GENOME

# Convert GFF to BED format
java -jar $SCRIPTMANAGER coordinate-manipulation gff-to-bed FIMO/$TF/fimo.gff -o FIMO/$TF/fimo.bed

# Assign a unique ID to the id column
awk '{OFS="\t"}{FS="\t"}{print $1,$2,$3,$1"_"$2"_"$3,$5,$6}' FIMO/$TF/fimo.bed > FIMO/$TF/fimo_uid.bed

# Filter out non-standard chromosomes
awk 'BEGIN{OFS="\t"}{if ($1 !~ /alt|random|chrUn/) print}' FIMO/$TF/fimo_uid.bed > FIMO/$TF/M1_unsorted_noproximity_unfiltered.bed

# Filter out motifs that overlap with blacklist regions
bedtools intersect -v -a FIMO/$TF/M1_unsorted_noproximity_unfiltered.bed -b $BLACKLIST > FIMO/$TF/M1_unsorted_noproximity.bed

# Apply proximity filter
java -jar $SCRIPTMANAGER peak-analysis filter-bed FIMO/$TF/M1_unsorted_noproximity.bed -o FIMO/$TF/M1_unsorted

# Rename motif-1 file
mv FIMO/$TF/M1_unsorted-FILTER.bed FIMO/$TF/$TF\_M1_unsorted.bed

# Clean-up
rm FIMO/$TF/M1_uid.bed
rm FIMO/$TF/M1_unsorted_noproximity_unfiltered.bed
rm FIMO/$TF/M1_unsorted_noproximity.bed
rm FIMO/$TF/M1_unsorted-CLUSTER.bed
#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 20:00:00
#SBATCH -A open
#SBATCH -o logs/2_FIMO_Motifs_from_Genome.log.out-%a
#SBATCH -e logs/2_FIMO_Motifs_from_Genome.log.err-%a
#SBATCH --array 1-5

# FIMO the reference genome for each motif in the PWM directory

### CHANGE ME
EXCLUSION=100
WRK=/storage/group/bfp2/default/hxc585_HainingChen/Fox_NFIA_CTCF/
CALL_RefPT_hg38=$WRK/02_Call_RefPT/
cd $CALL_RefPT_hg38


# Dependencies
# - bedtools
# - java
# - MEME suite (FIMO)

set -exo
module load anaconda3
source activate meme
source activate bioinfo

# Inputs and outputs
GENOME=$WRK/data/hg38_files/hg38.fa
BLACKLIST=$WRK/data/hg38_files/hg38-blacklist.bed
ORIGINAL_SCRIPTMANAGER=$WRK/bin/ScriptManager-v0.15.jar
SCRIPTMANAGER=$WRK/bin/ScriptManager-v0.15-$SLURM_ARRAY_TASK_ID.jar
cp $ORIGINAL_SCRIPTMANAGER $SCRIPTMANAGER

# Define PWM file path based on SLURM_ARRAY_TASK_ID index
PWM=$(ls PWM/{CTCF_M1.meme.txt,NFIA_M1.meme.txt,FOXA2_M1.meme.txt} | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

# Parse TF from PWM filename
TF=$(basename "$PWM" "_M1.meme.txt")


# Create output directories if they don't exist
[ -d logs ] || mkdir logs
[ -d FIMO ] || mkdir FIMO
[ -d FIMO/$TF ] || mkdir FIMO/$TF

echo "($SLURM_ARRAY_TASK_ID) $TF"

# Run FIMO to scan genome for motif occurrences
 fimo --verbosity 1 --thresh 1.0E-4 --oc FIMO/$TF $PWM $GENOME

grep $ALIAS FIMO/$TF/fimo.gff > FIMO/$TF/fimo.$TF.motif1.unsorted.noproximity.unfiltered.gff

# Convert GFF to BED format
java -jar $SCRIPTMANAGER coordinate-manipulation gff-to-bed FIMO/$TF/fimo.gff -o FIMO/$TF/motif1_unsorted_noproximity_unfiltered.bed

# Filter out motifs that overlap with blacklist regions
 tail -n +2 FIMO/$TF/motif1_unsorted_noproximity_unfiltered.bed | awk '{OFS="\t"}{FS="\t"}{print $1,$2,$3,$1"_"$2"_"$3,$5,$6}' - | bedtools intersect -v -a - -b $BLACKLIST | \
awk -v TF="${TF}" '{if ($6 == "+") print $0 > "FIMO/" TF "/motif1_unsorted_noproximity+.bed"; else print $0 > "FIMO/" TF "/motif1_unsorted_noproximity-.bed"}' 
bedtools intersect -v -a FIMO/$TF/motif1_unsorted_noproximity-.bed -b FIMO/$TF/motif1_unsorted_noproximity+.bed | cat FIMO/$TF/motif1_unsorted_noproximity+.bed - > FIMO/$TF/motif1_unsorted_noproximity.bed
rm  FIMO/$TF/motif1_unsorted_noproximity-.bed FIMO/$TF/motif1_unsorted_noproximity+.bed 

# Apply proximity filter
java -jar $SCRIPTMANAGER peak-analysis filter-bed  FIMO/$TF/motif1_unsorted_noproximity.bed -o FIMO/$TF/motif1_unsorted

# Rename motif-1 file and cleanup
mv FIMO/$TF/motif1_unsorted-FILTER.bed FIMO/$TF/$TF\_motif1_unsorted.bed
rm FIMO/$TF/motif1_unsorted_noproximity_unfiltered.bed FIMO/$TF/motif1_unsorted_noproximity.bed  FIMO/$TF/motif1_unsorted-CLUSTER.bed

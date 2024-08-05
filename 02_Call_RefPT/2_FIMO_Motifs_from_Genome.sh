#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 20:00:00
#SBATCH -A open
#SBATCH -o logs/2_FIMO_Motifs_from_Genome.log.out-%a
#SBATCH -e logs/2_FIMO_Motifs_from_Genome.log.err-%a
#SBATCH --array 1-4

# FIMO the reference genome for each motif in the PWM directory

### CHANGE ME
EXCLUSION=100
WRK=/storage/group/bfp2/default/hxc585_HainingChen/
02_Call_RefPT=$WRK/Fox_NFIA_CTCF/02_Call_RefPT/
cd 02_Call_RefPT
###

# Dependencies
# - bedtools
# - java
# - MEME suite (FIMO)

set -exo
module load anaconda3
#source activate meme
source activate bioinfo

# Inputs and outputs
GENOME=$WRK/2023_Chen_PIC3/hg19_files/hg19.fa
BLACKLIST=$WRK/2023_Chen_PIC3/hg19_files/hg19_exclude.bed
ORIGINAL_SCRIPTMANAGER=$WRK/2023_Chen_PIC3/bin/ScriptManager-v0.14.jar
SCRIPTMANAGER=$WRK/2023_Chen_PIC3/bin/ScriptManager-v0.14-$SLURM_ARRAY_TASK_ID.jar
cp $ORIGINAL_SCRIPTMANAGER $SCRIPTMANAGER

# Define PWM file path based on SLURM_ARRAY_TASK_ID index
PWM=`ls PWM/{CTCF.meme.txt,NFIA.meme.txt,scrambledCTCF.meme.txt,scrambledNFIA.meme.txt}  | head -n $SLURM_ARRAY_TASK_ID | tail -1`

# Parse TF from PWM filename
TF=`basename $PWM ".meme.txt"`

# Create output directories if they don't exist
[ -d logs ] || mkdir logs
[ -d FIMO ] || mkdir FIMO
[ -d FIMO/$TF ] || mkdir FIMO/$TF

echo "($SLURM_ARRAY_TASK_ID) $TF"

# Run FIMO to scan genome for motif occurrences
 fimo --verbosity 1 --thresh 1.0E-4 --oc FIMO/$TF $PWM $GENOME
# Extract motif-1 predictions using grep
#ALIAS='ID=1-MEME'
ALIAS='Alias=MEME-1;'
[[ $TF == "scrambled"* ]] && ALIAS='Alias=MEME;'
grep $ALIAS FIMO/$TF/fimo.gff > FIMO/$TF/fimo.$TF.motif1.unsorted.noproximity.unfiltered.gff

# Convert GFF to BED format
java -jar $SCRIPTMANAGER coordinate-manipulation gff-to-bed FIMO/$TF/fimo.$TF.motif1.unsorted.noproximity.unfiltered.gff -o FIMO/$TF/motif1_unsorted_noproximity_unfiltered.bed

# Filter out motifs that overlap with blacklist regions
awk '{OFS="\t"}{FS="\t"}{print $1,$2,$3,$1"_"$2"_"$3,$5,$6}' FIMO/$TF/motif1_unsorted_noproximity_unfiltered.bed | bedtools intersect -v -a - -b $BLACKLIST | \
awk -v TF="${TF}" '{if ($6 == "+") print $0 > "FIMO/" TF "/motif1_unsorted_noproximity+.bed"; else print $0 > "FIMO/" TF "/motif1_unsorted_noproximity-.bed"}' 
bedtools intersect -v -a FIMO/$TF/motif1_unsorted_noproximity-.bed -b FIMO/$TF/motif1_unsorted_noproximity+.bed | cat FIMO/$TF/motif1_unsorted_noproximity+.bed - > FIMO/$TF/motif1_unsorted_noproximity.bed
rm  FIMO/$TF/motif1_unsorted_noproximity-.bed FIMO/$TF/motif1_unsorted_noproximity+.bed 

# Apply proximity filter
java -jar $SCRIPTMANAGER peak-analysis filter-bed  FIMO/$TF/motif1_unsorted_noproximity.bed -o FIMO/$TF/motif1_unsorted

# Rename motif-1 file and cleanup
mv FIMO/$TF/motif1_unsorted-FILTER.bed FIMO/$TF/$TF\_motif1_unsorted.bed
rm FIMO/$TF/motif1_unsorted_noproximity_unfiltered.bed FIMO/$TF/motif1_unsorted_noproximity.bed  FIMO/$TF/motif1_unsorted-CLUSTER.bed

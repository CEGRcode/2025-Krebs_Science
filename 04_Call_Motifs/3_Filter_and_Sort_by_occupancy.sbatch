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
WRK=/path/to/2024-Krebs_Science/04_Call_Motifs
WRK=/scratch/owl5022/2024-Krebs_Science/04_Call_Motifs
###

# Dependencies
# - java
# - python
# - samtools

set -exo
module load anaconda3
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Inputs and outputs
BLACKLIST=../data/hg38_files/ENCFF356LFX_hg38_exclude.bed
MOTIF=../data/RefPT-Motif

# Script shortcuts
SCRIPTMANAGER=../bin/ScriptManager-v0.15.jar
INFLECTION=../bin/get_inflection_point.py

BAM_LIST=(
    ../data/BAM/K562_CTCF_BX_rep1_hg38.bam
    ../data/BAM/K562_FOXA1_BX_rep1_hg38.bam
    ../data/BAM/K562_FOXA2_BX_rep1_hg38.bam
    ../data/BAM/K562_NFIA_BX_rep1_hg38.bam
    ../data/BAM/HepG2_FOXA1_BX_rep1_hg38.bam
    ../data/BAM/HepG2_FOXA2_BX_rep1_hg38.bam
    ../data/BAM/K562_ZKSCAN1_BX_rep1_hg38.bam
    ../data/BAM/K562_HNF4A_BX_rep1_hg38.bam
)

# Define PWM file path based on SLURM_ARRAY_TASK_ID index
INDEX=$(($SLURM_ARRAY_TASK_ID-1))
BAMFILE=${BAM_LIST[$INDEX]}
BAM=`basename $BAMFILE ".bam"`

# Parse TF and Cell Line from BAM base filename
CELLLINE=`echo $BAM | awk -F"_" '{print $1}'`
TF=`echo $BAM | awk -F"_" '{print $2}'`
[[ $TF == "FOXA1" ]] && MTF=FOXA2 || MTF=$TF

TEMP=temp-3_Filter_and_Sort_by_occupancy
[ -d $TEMP ] || mkdir $TEMP

# [ -d $MOTIF/1bp ] || mkdir $MOTIF/1bp
# [ -d $MOTIF/1000bp ] || mkdir $MOTIF/1000bp

# Construct basename
BASE=$TEMP/${MTF}_${TF}-${CELLLINE}_M1_${WINDOW}bp

# =====Sort by occupancy=====

# Expand WINDOW
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c $WINDOW FIMO/$MTF/$MTF\_M1_unsorted.bed -o FIMO/$MTF/$MTF\_M1_unsorted_${WINDOW}bp.bed

# Tag Pileup
java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -1 -s 6 --combined FIMO/$MTF/$MTF\_M1_unsorted_${WINDOW}bp.bed $BAMFILE -M ${BASE}_unsorted

# Sort BED by Occupancy
java -jar $SCRIPTMANAGER coordinate-manipulation sort-bed -c $WINDOW FIMO/$MTF/$MTF\_M1_unsorted_${WINDOW}bp.bed ${BASE}_unsorted_combined.cdt -o ${BASE}_combined

# Visualize for reference
java -jar $SCRIPTMANAGER figure-generation heatmap --black -r 1 -l 2 -p .95 ${BASE}_combined.cdt -o ${BASE}_combined.png


# =====Call bound threshold=====

# Sum tags for Occupancy score
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum ${BASE}_combined.cdt -o ${BASE}_combined.out

# Add summed scores as 7th column to BED file
tail -n +2 ${BASE}_combined.out | cut -f 2 | paste ${BASE}_combined.bed - > ${BASE}_combined_7-Occupancy.bed

# Calculate inflection point
python $INFLECTION -i ${BASE}_combined.out -o ${BASE}_InflectionPT.png > ${BASE}_InflectionPT.txt

# Parse out threshold
THRESH=`grep "KneeLocator: " ${BASE}_InflectionPT.txt | awk '{print $2}'`
echo $THRESH

# Count sites above threshold
awk -v THRESH="$THRESH" 'BEGIN{FS="\t"}{if($2>=THRESH) print}' ${BASE}_combined.out > ${BASE}_AboveThreshold.out
NBOUND=`wc -l ${BASE}_AboveThreshold.out | awk '{print $1}'`

# Cut off BED at bound threshold
head -n $NBOUND ${BASE}_combined_7-Occupancy.bed > ${BASE}_7-Occupancy_BOUND.bed

# =====Expand RefPT=====

# Expand 1bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1 ${BASE}_7-Occupancy_BOUND.bed -o ${BASE}_7-Occupancy_BOUND_1bp.bed

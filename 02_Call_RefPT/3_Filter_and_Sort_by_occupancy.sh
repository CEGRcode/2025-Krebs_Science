#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 1:00:00
#SBATCH -A open
#SBATCH -o logs/3_Filter_and_Sort_by_occupancy.log.out-%a
#SBATCH -e logs/3_Filter_and_Sort_by_occupancy.log.err-%a
#SBATCH --array 1-5

# Sort filtered Motif BED Files by the first BX rep1 BAM file occupancy for each TF target (skip scrambled)

### CHANGE ME
WINDOW=100
# WRK=/storage/group/bfp2/default/hxc585_HainingChen/Fox_NFIA_CTCF/
WRK=/path/to/2024-Chen_Nature/02_Call_RefPT
WRK=/storage/home/owl5022/scratch/2024-Chen_Nature/02_Call_RefPT
###

# Dependencies
# - java
# - python
# - samtools

set -exo
module load anaconda3
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Inputs and outputs
BLACKLIST=../data/hg38_files/hg38-blacklist.bed
MOTIF=../data/RefPT-Motif

# Script shortcuts
SCRIPTMANAGER=../bin/ScriptManager-v0.15.jar

FIMO=FIMO
cd $FIMO || exit

[ -d FIMO ] || mkdir FIMO
[ -d $MOTIF/1bp ] || mkdir $MOTIF/1bp
[ -d $MOTIF/1000bp ] || mkdir $MOTIF/1000bp

# Parse TF from PWM filename
PWM=$(ls "$WRK/PWM/"*.meme.txt | head -n "$SLURM_ARRAY_TASK_ID" | tail -1)
TF=$(basename "$PWM" "_M1.meme.txt")

# TODO: Fix up filepaths below

# Get BAM file
# Define BAM file based on TF and sort by TF occupancy
if [[ $TF == "FOXA2" ]]; then
    HepG2FoxA1="$WRK/data/BAM/HepG2_FOXA1_BX_rep1_hg38.bam"
    HepG2FoxA2="$WRK/data/BAM/HepG2_FOXA2_BX_rep1_hg38.bam"
    K562FoxA1="$WRK/data/BAM/K562_FOXA1_BX_rep1_hg38.bam"
    K562FoxA2="$WRK/data/BAM/K562_FOXA2_BX_rep1_hg38.bam"

    java -jar "$SCRIPTMANAGER" coordinate-manipulation expand-bed -c "$WINDOW" "$TF/${TF}_motif1_unsorted.bed" -o "$TF/FOXA_MOTIF1_unsorted_${WINDOW}bp.bed"

    for FILE in "$HepG2FoxA1" "$HepG2FoxA2" "$K562FoxA1" "$K562FoxA2"; do
        BAM=$(basename "$FILE" ".bam")
        Cell=$(basename "$BAM" | cut -f1 -d "_")
        Factor=$(basename "$BAM" | cut -f2 -d "_")
        
        java -jar "$SCRIPTMANAGER" read-analysis tag-pileup -5 -1 -s 6 --combined "$TF/FOXA_MOTIF1_unsorted_${WINDOW}bp.bed" "$FILE" -M "$TF/FOXA_MOTIF1_unsorted_${WINDOW}bp_${BAM}"
        java -jar "$SCRIPTMANAGER" coordinate-manipulation sort-bed -c "$WINDOW" "$TF/FOXA_MOTIF1_unsorted_${WINDOW}bp.bed" "$TF/FOXA_MOTIF1_unsorted_${WINDOW}bp_${BAM}_combined.cdt" -o "$TF/FOXA_MOTIF1_${Cell}_${Factor}_sorted"
        java -jar "$SCRIPTMANAGER" figure-generation heatmap --black -r 1 -l 2 -p .95 "$TF/FOXA_MOTIF1_${Cell}_${Factor}_sorted.cdt" -o "$TF/FOXA_MOTIF1_${Cell}_${Factor}_sorted.png"
        
        rm "$TF/FOXA_MOTIF1_unsorted_${WINDOW}bp_${BAM}_combined.cdt"
        
        mkdir -p "$TF/${Cell}_${Factor}_InflectionFilter"
        
        java -jar "$SCRIPTMANAGER" read-analysis aggregate-data --sum "$TF/FOXA_MOTIF1_${Cell}_${Factor}_sorted.cdt" -o "$TF/${Cell}_${Factor}_InflectionFilter/"
        tail -n +2 "$TF/${Cell}_${Factor}_InflectionFilter/FOXA_MOTIF1_${Cell}_${Factor}_sorted_SCORES.out" | cut -f 2 | paste "$TF/FOXA_MOTIF1_${Cell}_${Factor}_sorted.bed" - > "$TF/FOXA_MOTIF1_${Cell}_${Factor}_sorted_SCORES.bed"
        
        python ../scripts/get_inflection_point.py -i "$TF/${Cell}_${Factor}_InflectionFilter/FOXA_MOTIF1_${Cell}_${Factor}_sorted_SCORES.out" -o "$TF/${Cell}_${Factor}_InflectionFilter/FOXA_MOTIF1_${Cell}_${Factor}_sorted_SCORES_inflectionPT.png" > "$TF/${Cell}_${Factor}_InflectionFilter/FOXA_MOTIF1_${Cell}_${Factor}_sorted_Threshold.txt"
        
        THRESH=$(grep "KneeLocator: " "$TF/${Cell}_${Factor}_InflectionFilter/FOXA_MOTIF1_${Cell}_${Factor}_sorted_Threshold.txt" | awk '{print $2}')
        echo "$THRESH"
        
        awk -v THRESH="$THRESH" '{FS="\t"; OFS="\t"; if($2>=THRESH) print}' "$TF/${Cell}_${Factor}_InflectionFilter/FOXA_MOTIF1_${Cell}_${Factor}_sorted_SCORES.out" > "$TF/${Cell}_${Factor}_InflectionFilter/FOXA_MOTIF1_${Cell}_${Factor}_sorted_SCORES_Filter-$THRESH.out"
        NLINES=$(wc -l < "$TF/${Cell}_${Factor}_InflectionFilter/FOXA_MOTIF1_${Cell}_${Factor}_sorted_SCORES_Filter-$THRESH.out")
        head -n "$NLINES" "$TF/FOXA_MOTIF1_${Cell}_${Factor}_sorted_SCORES.bed" > "$TF/FOXA_MOTIF1_${Cell}_${Factor}_Occupancy.bed"
        java -jar "$SCRIPTMANAGER" coordinate-manipulation expand-bed -c 1 "$TF/FOXA_MOTIF1_${Cell}_${Factor}_Occupancy.bed" -o "$TF/FOXA_MOTIF1_${Cell}_${Factor}_Occupancy_1bp.bed"
    done
    
    rm "$TF/${TF}_MOTIF1_unsorted_${WINDOW}bp.bed"
elif [[ $TF == "shuffle"* ]]; then
    #sort -k5,5nr "$TF/${TF}_motif1_unsorted.bed" | head -n 10000 > "$TF/${TF}_top10k.bed"
    java -jar "$SCRIPTMANAGER" coordinate-manipulation expand-bed -c 1 "$TF/${TF}_top10k.bed" -o "$TF/${TF}_top10k_1bp.bed"
else
    BAMFILE="$WRK/data/BAM/K562_${TF}_BX_rep1_hg38.bam"
    BAM=$(basename "$BAMFILE" ".bam")
    
    java -jar "$SCRIPTMANAGER" coordinate-manipulation expand-bed -c "$WINDOW" "$TF/${TF}_motif1_unsorted.bed" -o "$TF/${TF}_MOTIF1_unsorted_${WINDOW}bp.bed"
    java -jar "$SCRIPTMANAGER" read-analysis tag-pileup -5 -1 -s 6 --combined "$TF/${TF}_MOTIF1_unsorted_${WINDOW}bp.bed" "$BAMFILE" -M "$TF/${TF}_MOTIF1_unsorted_${WINDOW}bp_$BAM"
    java -jar "$SCRIPTMANAGER" coordinate-manipulation sort-bed -c "$WINDOW" "$TF/${TF}_MOTIF1_unsorted_${WINDOW}bp.bed" "$TF/${TF}_MOTIF1_unsorted_${WINDOW}bp_${BAM}_combined.cdt" -o "$TF/${TF}_MOTIF1_sorted"
    java -jar "$SCRIPTMANAGER" figure-generation heatmap --black -r 1 -l 2 -p .95 "$TF/${TF}_MOTIF1_sorted.cdt" -o "$TF/${TF}_MOTIF1_sorted.png"
    
    rm "$TF/${TF}_MOTIF1_unsorted_${WINDOW}bp_${BAM}_combined.cdt" "$TF/${TF}_MOTIF1_unsorted_${WINDOW}bp.bed"
    
    mkdir -p "$TF/InflectionFilter"
    
    java -jar "$SCRIPTMANAGER" read-analysis aggregate-data --sum "$TF/${TF}_MOTIF1_sorted.cdt" -o "$TF/InflectionFilter/"
    tail -n +2 "$TF/InflectionFilter/${TF}_MOTIF1_sorted_SCORES.out" | cut -f 2 | paste "$TF/${TF}_MOTIF1_sorted.bed" - > "$TF/${TF}_MOTIF1_sorted_SCORES.bed"
    
    python ../scripts/get_inflection_point.py -i "$TF/InflectionFilter/${TF}_MOTIF1_sorted_SCORES.out" -o "$TF/InflectionFilter/${TF}_MOTIF1_sorted_SCORES_inflectionPT.png" > "$TF/InflectionFilter/${TF}_MOTIF1_sorted_Threshold.txt"
    
    THRESH=$(grep "KneeLocator: " "$TF/InflectionFilter/${TF}_MOTIF1_sorted_Threshold.txt" | awk '{print $2}')
    echo "$THRESH"
    
    awk -v THRESH="$THRESH" '{FS="\t"; OFS="\t"; if($2>=THRESH) print}' "$TF/InflectionFilter/${TF}_MOTIF1_sorted_SCORES.out" > "$TF/InflectionFilter/${TF}_MOTIF1_sorted_SCORES_Filter-$THRESH.out"
    NLINES=$(wc -l < "$TF/InflectionFilter/${TF}_MOTIF1_sorted_SCORES_Filter-$THRESH.out")
    head -n "$NLINES" "$TF/${TF}_MOTIF1_sorted_SCORES.bed" > "$TF/${TF}_Occupancy.bed"
    java -jar "$SCRIPTMANAGER" coordinate-manipulation expand-bed -c 1 "$TF/${TF}_Occupancy.bed" -o "$TF/${TF}_Occupancy_1bp.bed"
fi

# Expand to 1000bp window and 1 bp Window
# Uncomment if needed
# java -jar "$SCRIPTMANAGER" coordinate-manipulation expand-bed -c 1000 "$MOTIF/${TF}_Occupancy.bed" -o "$MOTIF/1000bp/${TF}_Occupancy_1000bp.bed"
# java -jar "$SCRIPTMANAGER" coordinate-manipulation expand-bed -c 1 "$MOTIF/${TF}_Occupancy.bed" -o "$MOTIF/1bp/${TF}_Occupancy_1bp.bed"

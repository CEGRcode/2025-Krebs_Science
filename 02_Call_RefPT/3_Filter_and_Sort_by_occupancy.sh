#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 1:00:00
#SBATCH -A open
#SBATCH -o logs/3_Filter_and_Sort_by_occupancy.log.out-%a
#SBATCH -e logs/3_Filter_and_Sort_by_occupancy.log.err-%a
#SBATCH --array 1-3

# Sort filtered Motif BED Files by the first BX rep1 BAM file occupancy for each TF target (skip scrambled)

### CHANGE ME
WINDOW=100
WRK=/storage/group/bfp2/default/hxc585_HainingChen/
FIMO=$WRK/Fox_NFIA_CTCF/02_Call_RefPT/FIMO/
cd $FIMO
###

# Dependencies
# - java
# - python
# - samtools

set -exo
module load anaconda3
#source activate bioinfo
source activate virtualenv

# Script shortcuts
SCRIPTMANAGER=$WRK/2023_Chen_PIC3/bin/ScriptManager-v0.14.jar

# Inputs and outputs
BLACKLIST=$WRK/2023_Chen_PIC3/data/hg19_files/hg19_exclude.bed
MOTIF=$WRK/Fox_NFIA_CTCF/data/RefPT-Motif


# Parse TF from PWM filename
PWM=`ls $WRK/Fox_NFIA_CTCF/02_Call_RefPT/PWM/*.meme.txt | head -n $SLURM_ARRAY_TASK_ID | tail -1`
TF=`basename $PWM ".meme.txt"`

# Skip if scrambled
[[ $TF == "scrambled"* ]] && exit

# Get BAM file
# Define BAM file based on TF and sort by TF occupancy
if [[ $TF == "FoxA1" ]]; then
    HepG2FoxA1=$WRK/Fox_NFIA_CTCF/data/BAM/HepG2_FoxA1_BX_rep1_hg19.bam
    HepG2FoxA2=$WRK/Fox_NFIA_CTCF/data/BAM/HepG2_FoxA2_BX_rep3_hg19.bam
    # HepG2FoxA3=$WRK/Fox_NFIA_CTCF/data/BAM/HepG2_FoxA3_BX_rep3_hg19.bam 
    # K562FoxA1=$WRK/Fox_NFIA_CTCF/data/BAM/K562_FoxA1_BX_rep1_hg19.bam
    K562FoxA2=$WRK/Fox_NFIA_CTCF/data/BAM/K562_FoxA2_BX_rep1_hg19.bam
    K562FoxA3=$WRK/Fox_NFIA_CTCF/data/BAM/K562_FoxA3_BX_rep1_hg19.bam 
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c $WINDOW $TF/$TF\_motif1_unsorted.bed -o $TF/$TF\_MOTIF1_unsorted_$WINDOW\bp.bed
    for FILE in $HepG2FoxA2 $K562FoxA2 $K562FoxA3 ; do
        BAM=$(basename $FILE ".bam")
        Cell=$(basename "$BAM" | cut -f1 -d "_")
        Factor=$(basename "$BAM" | cut -f2 -d "_")
        java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -1 --combined $TF/$TF\_MOTIF1_unsorted_$WINDOW\bp.bed $FILE -M $TF/$TF\_MOTIF1_unsorted_$WINDOW\bp_$BAM
        java -jar $SCRIPTMANAGER coordinate-manipulation sort-bed -c $WINDOW $TF/$TF\_MOTIF1_unsorted_$WINDOW\bp.bed $TF/$TF\_MOTIF1_unsorted_$WINDOW\bp_$BAM\_combined.cdt -o $TF/$TF\_MOTIF1_${Cell}_${Factor}_sorted
        java -jar $SCRIPTMANAGER figure-generation heatmap --black -r 1 -l 2 -p .95 $TF/$TF\_MOTIF1_${Cell}_${Factor}_sorted.cdt -o $TF/$TF\_MOTIF1_${Cell}_${Factor}_sorted.png
        rm $TF/$TF\_MOTIF1_unsorted_$WINDOW\bp_$BAM\_combined.cdt
        [ -d $TF/${Cell}_${Factor}_InflectionFilter ] || mkdir -p $TF/${Cell}_${Factor}_InflectionFilter
        java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum $TF/$TF\_MOTIF1_${Cell}_${Factor}_sorted.cdt -o $TF/${Cell}_${Factor}_InflectionFilter/
        python ../scripts/get_inflection_point.py -i $TF/${Cell}_${Factor}_InflectionFilter/$TF\_MOTIF1_${Cell}_${Factor}_sorted_SCORES.out -o $TF/${Cell}_${Factor}_InflectionFilter/$TF\_MOTIF1_${Cell}_${Factor}_sorted_SCORES_inflectionPT.png > $TF/${Cell}_${Factor}_InflectionFilter/$TF\_MOTIF1_${Cell}_${Factor}_sorted_Threshold.txt
        THRESH=`grep "KneeLocator: " $TF/${Cell}_${Factor}_InflectionFilter/$TF\_MOTIF1_${Cell}_${Factor}_sorted_Threshold.txt | awk '{print $2}'`
        echo $THRESH
        awk -v THRESH=$THRESH '{FS="\t"}{OFS="\t"}{if($2>=THRESH) print}' $TF/${Cell}_${Factor}_InflectionFilter/$TF\_MOTIF1_${Cell}_${Factor}_sorted_SCORES.out > $TF/${Cell}_${Factor}_InflectionFilter/$TF\_MOTIF1_${Cell}_${Factor}_sorted_SCORES_Filter-$THRESH.out
        NLINES=`wc -l $TF/${Cell}_${Factor}_InflectionFilter/$TF\_MOTIF1_${Cell}_${Factor}_sorted_SCORES_Filter-$THRESH.out | awk '{print $1}'`
        head -n $NLINES $TF/$TF\_MOTIF1_${Cell}_${Factor}_sorted.bed > $TF/$TF\_${Cell}_${Factor}_Occupancy.bed
    done
    rm $TF/$TF\_MOTIF1_unsorted_$WINDOW\bp.bed
else
    BAMFILE=$WRK/Fox_NFIA_CTCF/data/BAM/K562_$TF\_BX_rep1_hg19.bam
    BAM=$(basename $BAMFILE ".bam")
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c $WINDOW $TF/$TF\_motif1_unsorted.bed -o $TF/$TF\_MOTIF1_unsorted_$WINDOW\bp.bed
    java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -1 --combined $TF/$TF\_MOTIF1_unsorted_$WINDOW\bp.bed $BAMFILE -M $TF/$TF\_MOTIF1_unsorted_$WINDOW\bp_$BAM
    java -jar $SCRIPTMANAGER coordinate-manipulation sort-bed -c $WINDOW $TF/$TF\_MOTIF1_unsorted_$WINDOW\bp.bed $TF/$TF\_MOTIF1_unsorted_$WINDOW\bp_$BAM\_combined.cdt -o $TF/$TF\_MOTIF1_sorted
    java -jar $SCRIPTMANAGER figure-generation heatmap --black -r 1 -l 2 -p .95 $TF/$TF\_MOTIF1_sorted.cdt -o $TF/$TF\_MOTIF1_sorted.png
    rm $TF/$TF\_MOTIF1_unsorted_$WINDOW\bp_$BAM\_combined.cdt $TF/$TF\_MOTIF1_unsorted_$WINDOW\bp.bed
    [ -d $TF/InflectionFilter ] || mkdir -p $TF/InflectionFilter
    java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum $TF/$TF\_MOTIF1_sorted.cdt -o $TF/InflectionFilter/
    python ../scripts/get_inflection_point.py -i $TF/InflectionFilter/$TF\_MOTIF1_sorted_SCORES.out -o $TF/InflectionFilter/$TF\_MOTIF1_sorted_SCORES_inflectionPT.png > $TF/InflectionFilter/$TF\_MOTIF1_sorted_Threshold.txt
    THRESH=`grep "KneeLocator: " $TF/InflectionFilter/$TF\_MOTIF1_sorted_Threshold.txt | awk '{print $2}'`
    echo $THRESH
    awk -v THRESH=$THRESH '{FS="\t"}{OFS="\t"}{if($2>=THRESH) print}' $TF/InflectionFilter/$TF\_MOTIF1_sorted_SCORES.out > $TF/InflectionFilter/$TF\_MOTIF1_sorted_SCORES_Filter-$THRESH.out
    NLINES=`wc -l $TF/InflectionFilter/$TF\_MOTIF1_sorted_SCORES_Filter-$THRESH.out | awk '{print $1}'`
    head -n $NLINES $TF/$TF\_MOTIF1_sorted.bed > $TF/$TF\_Occupancy.bed
fi


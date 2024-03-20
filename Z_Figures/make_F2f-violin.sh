#!/bin/bash

# Make violin occupancy plot (F2f)

### CHANGE ME
WRK=/path/to/2024-Chen_Nature/Z_Figures
WRK=/storage/home/owl5022/scratch/2024-Chen_Nature/Z_Figures
###

# Dependencies
# - java
# - samtools

set -exo
module load samtools

# Fill in placeholder constants with your directories
BAMDIR=$WRK/../data/BAM
MOTIF=$WRK/../data/RefPT-Motif

# Setup ScriptManager for job array
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.14.jar

# Script shortcuts
VIOLIN=$WRK/../bin/haining_occupancy_violin.py

# Set up output directories
[ -d logs ] || mkdir logs
[ -d F2/f ] || mkdir -p F2/f


# For each replicate...
for BAMFILE in "$BAMDIR/K562_NFIA_BX_rep1_hg19.bam" "$BAMDIR/K562_NFIA_BX_rep2_hg19.bam";
do
    BAM=`basename $BAMFILE ".bam"`

    [ -f $BAM\_violin_data.txt ] && rm $BAM\_violin_data.txt

    # Across each of the 3 NFIA BED groups
    for BEDFILE in "$MOTIF/NFIA_NucSort-DOWNSTREAM_1000bp.bed" "$MOTIF/NFIA_NucSort-OVERLAP_1000bp.bed" "$MOTIF/NFIA_NucSort-UPSTREAM_1000bp.bed";
    do
        BED=`basename $BEDFILE ".bed"`
        BASE=$BAM\_$BED\_read1
        
        # Pileup (exo)
        java -jar $SCRIPTMANAGER read-analysis tag-pileup $BEDFILE $BAMFILE --cpu 4 -5 -1 --combined -M $BASE
        java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum $BASE\_combined.cdt -o $BASE\_combined.txt

        # add group field and print sum to violin_data storage file
        awk -v FIELD=$BED '{OFS="\t"}{awk print FIELD,$2}' $BASE\_combined.txt >> $BAM\_violin_data.txt
    done
done

# # Plot occupancies
# python $VIOLIN -i K562_NFIA_BX_rep1_hg19_violin_data.txt -o F2/f/Violin_Rep1.png
# python $VIOLIN -i K562_NFIA_BX_rep2_hg19_violin_data.txt -o F2/f/Violin_Rep2.png

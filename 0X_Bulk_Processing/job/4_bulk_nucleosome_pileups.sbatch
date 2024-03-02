#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 00:15:00
#SBATCH -A open
#SBATCH -o logs/4_bulk_nucleosome_pileups.log.out
#SBATCH -e logs/4_bulk_nucleosome_pileups.log.err
#SBATCH --array 1-18

# Pileup midpoint (heatmaps) and composites (endo cut sites) for the 3 nucleosome BAM files:
#  - data/BAM/K562_-_BI_rep1_hg19.bam
#  - data/BAM/K562_H3K27ac_BX_rep1_hg19.bam
#  - data/BAM/K562_H3K4me3_BX_rep1_hg19.bam
# across the below hardcoded set of BED files. No scaling.

### CHANGE ME
WRK=/path/to/2023-Chen_Benzonase-ChIP-exo/0X_Bulk_Processing
###

BED_LIST=(
    "$WRK/../data/BED/TSS94_WDR5_sort_1000bp.bed"
    "$WRK/../data/BED/TSS_WDR5_same_1000bp.bed"
    "$WRK/../data/BED/K562_CoPRO-expressed_Gene-refSeqTSS_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/CTCF_Occupancy_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/CTCF_Occupancy-BOTTOM_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/CTCF_Occupancy-MIDDLE_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/CTCF_Occupancy-TOP_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/GABPA_Occupancy_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/SP1_Occupancy_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/NFIA_Occupancy_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/WDR5_Occupancy_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/WDR5_TSS_oppo_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/NFIA_NucSort_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/NFIA_NucSort-DOWNSTREAM_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/NFIA_NucSort-OVERLAP_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/NFIA_NucSort-UPSTREAM_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/scrambledCTCF_NucSort_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/scrambledNFIA_NucSort_1000bp.bed"
)

# Dependencies
# - java
# - samtools

set -exo
module load samtools

# Fill in placeholder constants with your directories
BAMDIR=$WRK/../data/BAM
OUTDIR=$WRK/Library

# Setup ScriptManager for job array
ORIGINAL_SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.14.jar
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.14-$SLURM_ARRAY_TASK_ID.jar
cp $ORIGINAL_SCRIPTMANAGER $SCRIPTMANAGER

# Set up output directories
[ -d $OUTDIR ] || mkdir $OUTDIR
[ -d logs ] || mkdir logs

# Determine BED file for the current job array index
BEDFILE=${BED_LIST[SLURM_ARRAY_TASK_ID]}
BED=`basename $BAMFILE ".bed"`

# Process for each BAM file in hardcoded list:
for BAMFILE in "$BAMDIR/K562_-_BI_rep1_hg19.bam" "$BAMDIR/K562_H3K27ac_BX_rep1_hg19.bam" "$BAMDIR/K562_H3K4me3_BX_rep1_hg19.bam";
do
	echo $BEDFILE "x" $BAMFILE

	# BAMFILE set-up
	BAM=`basename $BAMFILE ".bam"`
	[ -f $BAMFILE.bai ] || samtools index $BAMFILE

	DIR=$OUTDIR/$BED
	[ -d $DIR ] || mkdir $DIR
	[[ -d $DIR/CDT ]] || mkdir $DIR/CDT
	[[ -d $DIR/COMPOSITES ]] || mkdir $DIR/COMPOSITES
	[[ -d $DIR/PNG/Strand ]] || mkdir -p $DIR/PNG/Strand
	[[ -d $DIR/SVG ]] || mkdir $DIR/SVG

	echo "Run Nucleosome mid pileup"
	BASE=$BAM\_$BED\_midpoint

	# Pileup (mid)
    java -jar $SCRIPTMANAGER read-analysis tag-pileup $BEDFILE $BAMFILE --cpu 4 -m -a --combined -o $DIR/COMPOSITES/$BASE.out -M $DIR/CDT/$BASE

	# Two-color heatmap
	java -jar $SCRIPTMANAGER figure-generation heatmap -p .95 --black $DIR/CDT/$BASE\_combined.cdt -o $DIR/PNG/$BASE\_combined_treeview.png
        
    #Label heatmap
    java -jar $SCRIPTMANAGER figure-generation label-heatmap $DIR/PNG/$BASE\_combined_treeview.png -f 20 -l -500 -m 0 -r 500 -o $DIR/PNG/$BASE\_label.svg

	echo "Run Nucleosome endo cut sites (read2) pileup"
	BASE=$BAM\_$BED\_read2

	# Pileup (endo)
    java -jar $SCRIPTMANAGER read-analysis tag-pileup $BEDFILE $BAMFILE --cpu 4 -2 -o $DIR/COMPOSITES/$BASE.out
done
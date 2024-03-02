#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 00:15:00
#SBATCH -A open
#SBATCH -o logs/3_bulk_IgG_read1_pileups.log.out
#SBATCH -e logs/3_bulk_IgG_read1_pileups.log.err
#SBATCH --array 1-18

# Pileup read1 (exo cut sites) for the 2 IgG BAM files:
#  - data/BAM/K562_IgG_BX_merge.bam
#  - data/BAM/K562_IgG_XO_merge.bam
# across the below hardcoded set of BED files. No scaling.


### CHANGE ME
WRK=/path/to/2023-Chen_Benzonase-ChIP-exo/0X_Bulk_Processing
###

BED_LIST=(
    "$WRK/../data/BED/K562_CoPRO-expressed_Gene-refSeqTSS_1000bp.bed"
    "$WRK/../data/BED/TSS94_WDR5_sort_1000bp.bed"
    "$WRK/../data/BED/TSS_WDR5_same_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/CTCF_Occupancy_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/GABPA_Occupancy_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/SP1_Occupancy_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/NFIA_Occupancy_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/WDR5_Occupancy_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/WDR5_TSS_oppo_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/E4F1_Occupancy_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/MEIS2_Occupancy_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/NRF1_Occupancy_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/SRF1_Occupancy_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/USF1_Occupancy_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/USF2_Occupancy_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/YY1_Occupancy_1000bp.bed"
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

# Script shortcuts
# COMPOSITE=$WRK/bin/sum_Col_CDT.pl

# Set up output directories
[ -d $OUTDIR ] || mkdir $OUTDIR
[ -d logs ] || mkdir logs

# Determine BED file for the current job array index
BEDFILE=${BED_LIST[SLURM_ARRAY_TASK_ID]}
BED=`basename $BAMFILE ".bed"`

# Process for each BAM file in hardcoded list:
for BAMFILE in "$BAMDIR/K562_IgG_BX_merge.bam" "$BAMDIR/K562_IgG_XO_merge.bam";
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

	echo "Run IgG read1 pileup"
	BASE=$BAM\_$BED\_read1

	# Pileup (read1)
    java -jar $SCRIPTMANAGER read-analysis tag-pileup $BEDFILE $BAMFILE --cpu 4 -5 -1 \
		-o $DIR/COMPOSITES/$BASE.out -M $DIR/CDT/$BASE
	# no scaling

	# Two-color heatmap
	java -jar $SCRIPTMANAGER figure-generation heatmap -a 5 --blue $DIR/CDT/$BASE\_sense.cdt -o $DIR/PNG/$BASE\_sense_treeview.png
	java -jar $SCRIPTMANAGER figure-generation heatmap -a 5 --red $DIR/CDT/$BASE\_anti.cdt -o $DIR/PNG/$BASE\_anti_treeview.png
	java -jar $SCRIPTMANAGER figure-generation merge-heatmap $DIR/PNG/Strand/$BASE\_sense_treeview.png $DIR/PNG/Strand/$BASE\_anti_treeview.png -o $DIR/PNG/$BASE\_merge.png
 
    #Label heatmap
    java -jar $SCRIPTMANAGER figure-generation label-heatmap $DIR/PNG/$BASE\_merge.png \
		-f 20 -l -500 -m 0 -r 500 -o $DIR/PNG/$BASE\_label.svg
done
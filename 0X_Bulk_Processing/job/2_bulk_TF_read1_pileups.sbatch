#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 00:40:00
#SBATCH -A open
#SBATCH -o logs/2_bulk_TF_read1_pileups.log.out
#SBATCH -e logs/2_bulk_TF_read1_pileups.log.err
#SBATCH --array 1-7

# Generalized read 1 and read 2 TF pileups of BAM files against its
# corresponding motif file (sorted by Occupancy). With scaling.

### CHANGE ME
WRK=/path/to/2023-Chen_Benzonase-ChIP-exo/Z_Figures
###

# Dependencies
# - java
# - samtools

module load samtools
module load anaconda3

# Define the source files as TF_LIST
TF_LIST=(
	"CTCF"
	"E4F1"
	"MEIS2"
	"NRF1"
	"SRF1"
	"USF1"
	"USF2"
	"WDR5"
	"YY1"
)


# Fill in placeholder constants with your directories
BAMDIR=$WRK/../data/BAM
MOTIF=$WRK/../data/RefPT-Motif
OUTDIR=$WRK/SF1/b

# Setup ScriptManager for job array
ORIGINAL_SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.14.jar
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.14-$PBS_ARRAYID.jar
cp $ORIGINAL_SCRIPTMANAGER $SCRIPTMANAGER

# Script shortcuts
COMPOSITE=$WRK/../bin/sum_Col_CDT.pl

# Set up output directories
[ -d logs ] || mkdir logs
[ -d $OUTDIR ] || mkdir -p $OUTDIR

# Determine TF BED file for the current job array index
TF=${TF_LIST[SLURM_ARRAY_TASK_ID]}
BEDFILE=$MOTIF/1000bp/$TF\_Occupancy_1000bp.bed
BED=$TF\_Occupancy_1000bp

# Change TF for BAM filename construction (RBBP5 data aligned to WDR5 motif)
[[ "$TF" == "WDR5" ]] && TF="RBBP5"

# Process for each of 5 BAM files (TF_XO_rep1, TF_BX_rep1, TF_BX_rep2, IgG_XO_merge, IgG_BX_merge)
for BAMFILE in "$BAMDIR/K562_$TF\_XO_rep1_hg19.bam" "$BAMDIR/K562_$TF\_BX_rep1_hg19.bam" "$BAMDIR/K562_$TF\_BX_rep2_hg19.bam" "$BAMDIR/K562_IgG_XO_merge_hg19.bam" "$BAMDIR/K562_IgG_BX_merge_hg19.bam";
do
	echo $BEDFILE "x" $BAMFILE

	# BAMFILE set-up
	BAM=`basename $BAMFILE ".bam"`
	[ -f $BAMFILE.bai ] || samtools index $BAMFILE

	# Output set-up
	DIR=$OUTDIR/$BED
	[ -d $DIR ] || mkdir $DIR
	[[ -d $DIR/CDT ]] || mkdir $DIR/CDT
	[[ -d $DIR/COMPOSITES ]] || mkdir $DIR/COMPOSITES
	[[ -d $DIR/NormComposites ]] || mkdir $DIR/NormComposites
	[[ -d $DIR/PNG/Strand ]] || mkdir -p $DIR/PNG/Strand
	[[ -d $DIR/SVG ]] || mkdir $DIR/SVG

	# Get scaling factor
	FACTOR=`grep 'Scaling factor' $WRK/data/BAM/NormalizationFactors/$BAM\_NCISb_ScalingFactors.out | awk -F" " '{print $3}'`

	echo "Run generalized read1 pileup"
	BASE=$BAM\_$BED\_read1

	# Pileup (read 1)
	java -jar $SCRIPTMANAGER read-analysis tag-pileup $BEDFILE $BAMFILE --cpu 4 -5 -1 -o $DIR/COMPOSITES/$BASE.out -M $DIR/CDT/$BASE
	java -jar $SCRIPTMANAGER read-analysis scale-matrix $DIR/CDT/$BASE\_anti.cdt -s $FACTOR -o $DIR/CDT/$BASE\_anti_Normalized.cdt
	java -jar $SCRIPTMANAGER read-analysis scale-matrix $DIR/CDT/$BASE\_sense.cdt -s $FACTOR -o $DIR/CDT/$BASE\_sense_Normalized.cdt

	# Heatmaps
	java -jar $SCRIPTMANAGER figure-generation heatmap -a 5 --blue $DIR/CDT/$BASE\_sense_Normalized.cdt -o $DIR/PNG/Strand/$BASE\_sense_treeview.png
	java -jar $SCRIPTMANAGER figure-generation heatmap -a 5 --red $DIR/CDT/$BASE\_anti_Normalized.cdt -o $DIR/PNG/Strand/$BASE\_anti_treeview.png
	java -jar $SCRIPTMANAGER figure-generation merge-heatmap $DIR/PNG/Strand/$BASE\_sense_treeview.png $DIR/PNG/Strand/$BASE\_anti_treeview.png -o $DIR/PNG/$BASE\_merge.png

	# Label heatmap
	java -jar $SCRIPTMANAGER figure-generation label-heatmap $DIR/Heatmap/$BASE\_merge.png \
		-f 20 -l -500 -m 0 -r 500 -o $DIR/Heatmap/$BASE\_label.svg

	# Make stranded composites
	perl $COMPOSITE $DIR/CDT/$BASE\_anti_Normalized.cdt $DIR/CDT/$BASE\_ANTI
	perl $COMPOSITE $DIR/CDT/$BASE\_sense_Normalized.cdt $DIR/CDT/$BASE\_SENSE
	cat $DIR/CDT/$BASE\_ANTI <(tail -1 $DIR/CDT/$BASE\_SENSE) > $DIR/NormComposites/$BASE\_Normalized.out
	rm $DIR/CDT/$BASE\_SENSE $DIR/CDT/$BASE\_ANTI

	echo "Run generalized read2 pileup"
	BASE=$BAM\_$BED\_read2

	# Pileup (read 2)
	java -jar $SCRIPTMANAGER read-analysis tag-pileup $BEDFILE $BAMFILE --cpu 4 -5 -2 -o $DIR/COMPOSITES/$BASE.out -M $DIR/CDT/$BASE
	java -jar $SCRIPTMANAGER read-analysis scale-matrix $DIR/CDT/$BASE\_anti.cdt  -s $FACTOR -o $DIR/CDT/$BASE\_anti_Normalized.cdt
	java -jar $SCRIPTMANAGER read-analysis scale-matrix $DIR/CDT/$BASE\_sense.cdt -s $FACTOR -o $DIR/CDT/$BASE\_sense_Normalized.cdt

	# Make stranded composites
	perl $COMPOSITE $DIR/CDT/$BASE\_anti_Normalized.cdt $DIR/CDT/$BASE\_ANTI
	perl $COMPOSITE $DIR/CDT/$BASE\_sense_Normalized.cdt $DIR/CDT/$BASE\_SENSE
	cat $DIR/CDT/$BASE\_ANTI <(tail -1 $DIR/CDT/$BASE\_SENSE) > $DIR/NormComposites/$BASE\_Normalized.out
	rm $DIR/CDT/$BASE\_SENSE $DIR/CDT/$BASE\_ANTI

done
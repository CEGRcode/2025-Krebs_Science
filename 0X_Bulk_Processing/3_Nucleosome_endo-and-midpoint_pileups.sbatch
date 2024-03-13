#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 00:15:00
#SBATCH -A open
#SBATCH -o logs/3_Nucleosome_endo-and-midpoint_pileups.log.out-%a
#SBATCH -e logs/3_Nucleosome_endo-and-midpoint_pileups.log.err-%a
#SBATCH --array 1-48

# Pileup midpoint (heatmaps) and composites (endo cut sites) for a custom
# combination of BAM x BED files (primarily nucleosome BAM). No scaling.
#  - data/BAM/K562_-_BI_rep1_hg19.bam
#  - data/BAM/K562_H3K27ac_BX_rep1_hg19.bam
#  - data/BAM/K562_H3K4me3_BX_rep1_hg19.bam
# Read2 sense/anti flipped for plotting.

### CHANGE ME
WRK=/path/to/2024-Chen_Nature/0X_Bulk_Processing
WRK=/storage/home/owl5022/scratch/2024-Chen_Nature/0X_Bulk_Processing
METADATA=Nucleosome_pileups.txt
###

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
[ -d logs ] || mkdir logs
[ -d $OUTDIR ] || mkdir $OUTDIR

# Determine BAM file for the current job array index
BAMFILE=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $1}'`
BAM=`basename $BAMFILE ".bam"`
[ -f $BAMFILE.bai ] || samtools index $BAMFILE

# Determine BED file for the current job array index
BEDFILE=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $2}'`
BED=`basename $BEDFILE ".bed"`

echo "(${SLURM_ARRAY_TASK_ID}) ${BEDFILE} "x" ${BAMFILE} "

DIR=$OUTDIR/$BED
[ -d $DIR ] || mkdir $DIR
[[ -d $DIR/CDT ]] || mkdir $DIR/CDT
[[ -d $DIR/Composites ]] || mkdir $DIR/Composites
[[ -d $DIR/PNG/Strand ]] || mkdir -p $DIR/PNG/Strand
[[ -d $DIR/SVG ]] || mkdir $DIR/SVG

echo "Run Nucleosome mid pileup"
BASE=$BAM\_$BED\_midpoint

# Pileup (mid)
java -jar $SCRIPTMANAGER read-analysis tag-pileup $BEDFILE $BAMFILE --cpu 4 -m --combined -o $DIR/Composites/$BASE.out -M $DIR/CDT/$BASE
# No scaling

# Two-color heatmap
java -jar $SCRIPTMANAGER figure-generation heatmap -p .95 --black $DIR/CDT/$BASE\_combined.cdt -o $DIR/PNG/$BASE\_combined_treeview.png

# Count sites
NSITES=`wc -l $BEDFILE | awk '{print $1-1}'`

# Label heatmap
java -jar $SCRIPTMANAGER figure-generation label-heatmap $DIR/PNG/$BASE\_combined_treeview.png \
	-l "-500" -m "0" -r "+500" -w 1 -f 20 \
	-x $BED -y "$BED occurences (${NSITES} sites)" \
	-o $DIR/SVG/$BASE\_combined_treeview_label.svg

echo "Run Nucleosome endo cut sites (read2) pileup"
BASE=$BAM\_$BED\_read2

# Pileup (endo)
java -jar $SCRIPTMANAGER read-analysis tag-pileup $BEDFILE $BAMFILE --cpu 4 -5 -2 -M $DIR/CDT/$BASE

# Make stranded composites (need to swap sense/anti labels for locus plotter if read2)
perl $COMPOSITE $DIR/CDT/$BASE\_anti_Normalized.cdt $DIR/CDT/$BASE\_ANTI
perl $COMPOSITE $DIR/CDT/$BASE\_sense_Normalized.cdt $DIR/CDT/$BASE\_SENSE
head -1 $DIR/CDT/$BASE\_ANTI > $DIR/Composites/$BASE\_Normalized.out
paste <(echo $BASE\_Normalized_sense_RESTRANDED.cdt ) <(tail -1 $DIR/CDT/$BASE\_ANTI | cut -f2-10010 ) >> $DIR/Composites/$BASE\_Normalized.out
paste <(echo $BASE\_Normalized_anti_RESTRANDED.cdt ) <(tail -1 $DIR/CDT/$BASE\_SENSE | cut -f2-10010 ) >> $DIR/Composites/$BASE\_Normalized.out
rm $DIR/CDT/$BASE\_SENSE $DIR/CDT/$BASE\_ANTI
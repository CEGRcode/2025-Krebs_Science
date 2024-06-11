#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=24gb
#SBATCH -t 00:40:00
#SBATCH -A open
#SBATCH -o logs/2_exo_read_pileups.log.out-%a
#SBATCH -e logs/2_exo_read_pileups.log.err-%a
#SBATCH --array 1-100

# Pileup read1 (exo cut sites) for a custom combinations of BAM x BED files.
# Heatmaps included and all with scaling. Read2 sense/anti flipped for plotting.

### CHANGE ME
WRK=/path/to/2024-Chen_Nature/0X_Bulk_Processing
WRK=/storage/home/owl5022/scratch/2024-Chen_Nature/0X_Bulk_Processing
METADATA=Read1_pileups.txt
###

# Dependencies
# - java
# - perl
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
COMPOSITE=$WRK/../bin/sum_Col_CDT.pl

# Set up output directories
[ -d logs ] || mkdir logs
[ -d $OUTDIR ] || mkdir $OUTDIR

# Determine BAM file for the current job array index and process each line of METADATA 
BAMFILE=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $1}'`
BAM=`basename $BAMFILE ".bam"`
[ -f $BAMFILE.bai ] || samtools index $BAMFILE
	while read line; do
  		echo "$line"
	done <$METADATA

# Determine BED file for the current job array index
BEDFILE=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $2}'`
BED= basename $BEDFILE ".bam"`
	
	while read line; do
  		echo "$line"
	done <$METADATA

echo "(${SLURM_ARRAY_TASK_ID}) ${BEDFILE} "x" ${BAMFILE} "

DIR=$OUTDIR/$BED
[ -d $DIR ] || mkdir $DIR
[[ -d $DIR/CDT ]] || mkdir $DIR/CDT
[[ -d $DIR/Composites ]] || mkdir $DIR/Composites
[[ -d $DIR/PNG/Strand ]] || mkdir -p $DIR/PNG/Strand
[[ -d $DIR/SVG ]] || mkdir $DIR/SVG

echo "Run Custom read1 pileup"
BASE=$BAM\_$BED\_read1

# Pileup (read 1)
java -jar $SCRIPTMANAGER read-analysis tag-pileup $BEDFILE $BAMFILE --cpu 4 -5 -1 -o $DIR/Composites/$BASE.out -M $DIR/CDT/$BASE
java -jar $SCRIPTMANAGER read-analysis scale-matrix $DIR/CDT/$BASE\_anti.cdt  -s $FACTOR -o $DIR/CDT/$BASE\_anti_Normalized.cdt
java -jar $SCRIPTMANAGER read-analysis scale-matrix $DIR/CDT/$BASE\_sense.cdt -s $FACTOR -o $DIR/CDT/$BASE\_sense_Normalized.cdt

# Two-color heatmap & merge
java -jar $SCRIPTMANAGER figure-generation heatmap -a 5 --blue $DIR/CDT/$BASE\_sense_Normalized.cdt -o $DIR/PNG/Strand/$BASE\_sense_Normalized_treeview.png
java -jar $SCRIPTMANAGER figure-generation heatmap -a 5 --red  $DIR/CDT/$BASE\_anti_Normalized.cdt  -o $DIR/PNG/Strand/$BASE\_anti_Normalized_treeview.png
java -jar $SCRIPTMANAGER figure-generation merge-heatmap $DIR/PNG/Strand/$BASE\_sense_Normalized_treeview.png $DIR/PNG/Strand/$BASE\_anti_Normalized_treeview.png -o $DIR/PNG/$BASE\_Normalized_merge.png

# Count sites
NSITES=`wc -l $BEDFILE | awk '{print $1-1}'`

# Label heatmap
java -jar $SCRIPTMANAGER figure-generation label-heatmap $DIR/PNG/$BASE\_Normalized_merge.png \
	-l "-500" -m "0" -r "+500" -w 1 -f 20 \
	-x $BED -y "$BED occurences (${NSITES} sites)" \
	-o $DIR/SVG/$BASE\_Normalized_merge_label.svg

# Make stranded composites
perl $COMPOSITE $DIR/CDT/$BASE\_anti_Normalized.cdt $DIR/CDT/$BASE\_ANTI
perl $COMPOSITE $DIR/CDT/$BASE\_sense_Normalized.cdt $DIR/CDT/$BASE\_SENSE
cat $DIR/CDT/$BASE\_ANTI <(tail -1 $DIR/CDT/$BASE\_SENSE) > $DIR/Composites/$BASE\_Normalized.out
rm $DIR/CDT/$BASE\_SENSE $DIR/CDT/$BASE\_ANTI

echo "Run Custom read2 pileup"
BASE=$BAM\_$BED\_read2

# Pileup (read 2)
java -jar $SCRIPTMANAGER read-analysis tag-pileup $BEDFILE $BAMFILE --cpu 4 -5 -2 -o $DIR/Composites/$BASE.out -M $DIR/CDT/$BASE
java -jar $SCRIPTMANAGER read-analysis scale-matrix $DIR/CDT/$BASE\_anti.cdt  -s $FACTOR -o $DIR/CDT/$BASE\_anti_Normalized.cdt
java -jar $SCRIPTMANAGER read-analysis scale-matrix $DIR/CDT/$BASE\_sense.cdt -s $FACTOR -o $DIR/CDT/$BASE\_sense_Normalized.cdt

# Make stranded composites (need to swap sense/anti labels for locus plotter if read2)
perl $COMPOSITE $DIR/CDT/$BASE\_anti_Normalized.cdt $DIR/CDT/$BASE\_ANTI
perl $COMPOSITE $DIR/CDT/$BASE\_sense_Normalized.cdt $DIR/CDT/$BASE\_SENSE
head -1 $DIR/CDT/$BASE\_ANTI > $DIR/Composites/$BASE\_Normalized.out
paste <(echo $BASE\_Normalized_sense_RESTRANDED.cdt ) <(tail -1 $DIR/CDT/$BASE\_ANTI | cut -f2-10010 ) >> $DIR/Composites/$BASE\_Normalized.out
paste <(echo $BASE\_Normalized_anti_RESTRANDED.cdt ) <(tail -1 $DIR/CDT/$BASE\_SENSE | cut -f2-10010 ) >> $DIR/Composites/$BASE\_Normalized.out
rm $DIR/CDT/$BASE\_SENSE $DIR/CDT/$BASE\_ANTI
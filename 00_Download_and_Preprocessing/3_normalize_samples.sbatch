#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 00:40:00
#SBATCH -A open
#SBATCH -o logs/3_normalize_samples.log.out-%a
#SBATCH -e logs/3_normalize_samples.log.err-%a
#SBATCH --array 1-20

# Calculate normalization factors and get fragment distribution for every *.bam file in a directory

### CHANGE ME
WRK=/path/to/2024-Krebs_Science/00_Download_and_Preprocessing
WRK=/storage/home/owl5022/scratch/2024-Krebs_Science/00_Download_and_Preprocessing
###

# Dependencies
# - java
# - samtools

set -exo
module load samtools

# Fill in placeholder constants with your directories
BAMDIR=$WRK/../data/BAM
NDIR=$BAMDIR/NormalizationFactors
FDIR=$BAMDIR/FragmentDistributions
BLACKLIST=$WRK/../data/hg38_files/hg38-blacklist.bed

# Setup ScriptManager for job array
ORIGINAL_SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15-$SLURM_ARRAY_TASK_ID.jar
cp $ORIGINAL_SCRIPTMANAGER $SCRIPTMANAGER

# Script shortcuts
# (none)

# Set up output directories
[ -d logs ] || mkdir logs
[ -d $NDIR ] || mkdir $NDIR
[ -d $FDIR ] || mkdir $FDIR

# Determine BAM file for the current job array index
BAMFILE=`ls $BAMDIR/*hg38.bam | head -n $SLURM_ARRAY_TASK_ID | tail -1`
BAM=`basename $BAMFILE ".bam"`
[ -f $BAMFILE.bai ] || samtools index $BAMFILE

echo "Calculate Total Tag normalization factors w/ blacklist"

# Calculate Total Tag normalization factor
java -jar -Djava.awt.headless=true $SCRIPTMANAGER read-analysis scaling-factor $BAMFILE -f $BLACKLIST --total-tag -o $NDIR/$BAM\_TotalTag

STRAIN=`echo $BAM | cut -d "_" -f 1`
TARGET=`echo $BAM | cut -d "_" -f 2`
[ "$TARGET" == "IgG" ] && exit
[ "$STRAIN" == "BNase-seq" ] && exit
[ "$STRAIN" == "DNase-seq" ] && exit
[ "$STRAIN" == "MNase-seq" ] && exit
[ "$STRAIN" == "MNase-ChIP" ] && exit
[ "$STRAIN" == "CUTRUN" ] && exit

CONTROL=$BAMDIR/${STRAIN}_IgG_BX_merge_hg38.bam
echo "Calculate NCIS normalization factors w/ blacklist and control: ${CONTROL}"

# Calculate NCIS normalization factor
java -jar -Djava.awt.headless=true $SCRIPTMANAGER read-analysis scaling-factor $BAMFILE -f $BLACKLIST --ncis -c $CONTROL -w 500 -o $NDIR/$BAM\_NCISb

# Calculate fragment size distribution
java -jar -Djava.awt.headless=true $SCRIPTMANAGER bam-statistics pe-stat -x 500 -d $BAMFILE -o $FDIR/$BAM

 
rm $SCRIPTMANAGER

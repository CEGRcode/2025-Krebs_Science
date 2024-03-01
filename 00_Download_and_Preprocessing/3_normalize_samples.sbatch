#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 00:10:00
#SBATCH -A open
#SBATCH -o logs/3_normalize_samples.log.out
#SBATCH -e logs/3_normalize_samples.log.err
#SBATCH --array 1-NSAMPLES

# Calculate normalization factors and get fragment distribution for every *.bam file in a directory

### CHANGE ME
WRK=/path/to/2023-Chen_Benzonase-ChIP-exo/00_Download_and_Preprocessing
###

# Dependencies
# - java
# - samtools

module load samtools

# Fill in placeholder constants with your directories
BAMDIR=$WRK/../data/BAM
BEDDIR=$WRK/../data/BED
NDIR=$WRK/../data/BAM/NormalizationFactors
FDIR=$WRK/../data/BAM/FragmentDistributions
BLACKLIST=$WRK/../data/GENOME/hg19_exclude.bed

# Setup ScriptManager for job array
ORIGINAL_SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.14.jar
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.14-$SLURM_ARRAY_TASK_ID.jar
cp $ORIGINAL_SCRIPTMANAGER $SCRIPTMANAGER

# Script shortcuts
# (none)

# Set up output directories
[ -d logs ] || mkdir logs
[ -d $NDIR ] || mkdir $NDIR
[ -d $FDIR ] || mkdir $FDIR

# Determine BAM file for the current job array index
BAMFILE=`ls $BAMDIR/*.bam | head -n $SLURM_ARRAY_TASK_ID | tail -1`;
BAM=`basename $BAMFILE ".bam"`
[ -f $BAMFILE.bai ] || samtools index $BAMFILE

# Parse out sample info
CLINE=`echo $BAM | awk -F"_" '{print $1}'`
ASSAY=`echo $BAM | awk -F"_" '{print $3}'`

# Build control filename based on sample info
CONTROL=$BAMDIR/$CLINE\_IgG_$ASSAY\_merge_hg19.bam

# Calculate fragment size distribution
java -jar $SCRIPTMANAGER bam-statistics pe-stat -x 500 -s $BAMFILE -o $FDIR/$BAM

# Use different normalization method depending on target/assay
echo "Calculate Total Tag normalization factors"
java -jar $SCRIPTMANAGER read-analysis scaling-factor $BAMFILE -o $NDIR/$BAM\_TotalTag
[ "$ASSAY" == "BI" ] && exit
echo "Calculate classic TF NCIS normalization factors w/ blacklist"
java -jar $SCRIPTMANAGER read-analysis scaling-factor $BAMFILE --ncis -c $CONTROL -f $BLACKLIST -w 100 -o $NDIR/$BAM\_NCISb  # Do we want to add a custom window?

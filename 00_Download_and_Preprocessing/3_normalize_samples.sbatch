#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 00:10:00
#SBATCH -A open
#SBATCH -o logs/3_normalize_samples.log.out-%a
#SBATCH -e logs/3_normalize_samples.log.err-%a
#SBATCH --array 1-20

# Calculate normalization factors and get fragment distribution for every *.bam file in a directory

### CHANGE ME
WRK=/storage/group/bfp2/default/hxc585_HainingChen/Fox_NFIA_CTCF/
###

# Dependencies
# - java
# - samtools

module load samtools

# Fill in placeholder constants with your directories
BAMDIR=$WRK/data/hg38_BAM
NDIR=$WRK/data/NormalizationFactors
FDIR=$WRK/data/FragmentDistributions
BLACKLIST=$WRK/data/hg38_files/hg38-blacklist.bed

# Setup ScriptManager for job array
ORIGINAL_SCRIPTMANAGER="$WRK/bin/ScriptManager-v0.15.jar"
SCRIPTMANAGER=$WRK/bin/ScriptManager-v0.15-$SLURM_ARRAY_TASK_ID.jar
cp $ORIGINAL_SCRIPTMANAGER $SCRIPTMANAGER

# Script shortcuts
# (none)

# Set up output directories
[ -d logs ] || mkdir logs
[ -d $NDIR ] || mkdir $NDIR
[ -d $FDIR ] || mkdir $FDIR

# Determine BAM file for the current job array index
BAMFILE=`ls $BAMDIR/*hg38.bam | head -n $SLURM_ARRAY_TASK_ID | tail -1`
BAMFILENEME=`basename $BAMFILE ".bam"`
[ -f $BAMFILENEME.bai ] || samtools index $BAMFILE
STRAIN=`basename $BAMFILE ".bam" | cut -d "_" -f 1`
Target=`basename $BAMFILE ".bam" | cut -d "_" -f 2`
ASSAY=`basename $BAMFILE ".bam" | cut -d "_" -f 3`

[ "$ASSAY" == "BI" ] && exit
echo "Calculate classic TF NCIS normalization factors w/ blacklist"
java -jar $SCRIPTMANAGER read-analysis scaling-factor $BAMFILE -f $BLACKLIST --ncis -c $BAMDIR/${STRAIN}_IgG_*_hg38.bam -w 500 -o $NDIR/$BAMFILENEME\_NCISb


# Calculate fragment size distribution
java -jar $SCRIPTMANAGER bam-statistics pe-stat -x 500 -s $BAMFILE -o $FDIR/$BAMFILENEME

 
rm $SCRIPTMANAGER
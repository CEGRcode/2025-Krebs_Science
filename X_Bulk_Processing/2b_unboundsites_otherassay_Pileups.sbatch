#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=24gb
#SBATCH -t 00:40:00
#SBATCH -A open
#SBATCH -o logs/4_Five_Read1_Pileups.log.out-%a
#SBATCH -e logs/4_Five_Read1_Pileups.log.err-%a
#SBATCH --array 1-19

# Pileup heatmaps and composites (5' read1, exo cut sites) for a custom
# combination of BAM x BED files (primarily nucleosome BAM). Configuration
# file specifies insert size filters, contrast threshold, and normalization
# strategy.

### CHANGE ME
WRK=/storage/group/bfp2/default/hxc585_HainingChen/2025_Chen_TF-Nuc/X_Bulk_Processing
#WRK=/ocean/projects/see180003p/owlang/2024-Krebs_Science/X_Bulk_Processing
#WRK=/storage/group/bfp2/default/owl5022-OliviaLang/2024-Krebs_Science/X_Bulk_Processing
#WRK=/scratch/owl5022/2024-Krebs_Science/X_Bulk_Processing
METADATA=$WRK/Unboundsites_otherassay.txt
THREADS=4
###

# Dependencies
# - java
# - samtools

set -exo
module load samtools

# Fill in placeholder constants with your directories
BAMDIR=$WRK/../data/BAM

# Setup ScriptManager for job array
ORIGINAL_SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15-$SLURM_ARRAY_TASK_ID.jar
cp $ORIGINAL_SCRIPTMANAGER $SCRIPTMANAGER

# Set up output directories
[ -d logs ] || mkdir logs
[ -d Library ] || mkdir Library

# Determine BAM file for the current job array index
MIN=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $1}'`
MAX=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $2}'
TP_PARAMS="--cpu ${THREADS} -5 -a"
FILTER=""
[[ $MIN =~ "NaN" ]] || FILTER=${FILTER}-MIN${MIN}
[[ $MIN =~ "NaN" ]] || TP_PARAMS="${TP_PARAMS} -n ${MIN}"
[[ $MAX =~ "NaN" ]] || FILTER=${FILTER}-MAX${MAX}
[[ $MAX =~ "NaN" ]] || TP_PARAMS="${TP_PARAMS} -x ${MAX}"

BAMFILE=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $3}'`
BAM=`basename $BAMFILE ".bam"`
[ -f $BAMFILE.bai ] || samtools index $BAMFILE

# Determine BED file for the current job array index
BEDFILE=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $4}'`
BED=`basename $BEDFILE ".bed"`

# ===============================================================================================================================

BASE=${BAM}_${BED}_5both${FILTER}

# Setup output directories
DIR=Library/$BED
[ -d $DIR ] || mkdir $DIR
[[ -d $DIR/Composites ]] || mkdir $DIR/Composites

# ===============================================================================================================================

echo "Run Custom read1 pileup"

# Pileup (read 1)
java -jar $SCRIPTMANAGER read-analysis tag-pileup $BEDFILE $BAMFILE ${TP_PARAMS} -o $DIR/Composites/$BASE.out
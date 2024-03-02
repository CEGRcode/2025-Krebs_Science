#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l pmem=14gb
#PBS -l walltime=00:40:00
#PBS -A open
#PBS -o logs/bulk_cross.log.out
#PBS -e logs/bulk_cross.log.err
#PBS -t 1-5

# Pileup read1 and read2 (exo & endo cut sites) for custom combinations of BAM x BED files
# that are filtered by fragment read size (>100bp). Just composites with scaling.

BAM_FILES=(
    "$WRK/../data/BAM/K562_NFIA_BX_rep1_hg19.bam"
    "$WRK/../data/BAM/K562_NFIA_BX_rep1_hg19.bam"
    "$WRK/../data/BAM/K562_NFIA_BX_rep1_hg19.bam"
    "$WRK/../data/BAM/K562_NFIA_BX_rep1_hg19.bam"
#4
    "$WRK/../data/BAM/K562_CTCF_BX_rep1_hg19.bam"
    "$WRK/../data/BAM/K562_CTCF_BX_rep1_hg19.bam"
    "$WRK/../data/BAM/K562_CTCF_BX_rep1_hg19.bam"
    "$WRK/../data/BAM/K562_CTCF_BX_rep1_hg19.bam"
    "$WRK/../data/BAM/K562_CTCF_BX_rep1_hg19.bam"
#9
    "$WRK/../data/BAM/K562_CTCF_BX_rep2_hg19.bam"
    "$WRK/../data/BAM/K562_CTCF_BX_rep2_hg19.bam"
    "$WRK/../data/BAM/K562_CTCF_BX_rep2_hg19.bam"
    "$WRK/../data/BAM/K562_CTCF_BX_rep2_hg19.bam"
    "$WRK/../data/BAM/K562_CTCF_BX_rep2_hg19.bam"
)

BED_FILES=(
    "$WRK/../data/RefPT-Motif/1000bp/NFIA_NucSort_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/NFIA_NucSort-UPSTREAM_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/NFIA_NucSort-OVERLAP_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/NFIA_NucSort-DOWNSTREAM_1000bp.bed"
#4
    "$WRK/../data/RefPT-Motif/1000bp/CTCF_Occupancy_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/CTCF_Occupancy_TOP_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/CTCF_Occupancy_MIDDLE_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/CTCF_Occupancy_BOTTOM_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/scrambledCTCF_NucSort_1000bp.bed"
#9
    "$WRK/../data/RefPT-Motif/1000bp/CTCF_Occupancy_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/CTCF_Occupancy_TOP_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/CTCF_Occupancy_MIDDLE_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/CTCF_Occupancy_BOTTOM_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/scrambledCTCF_NucSort_1000bp.bed"
)

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
ORIGINAL_SCRIPTMANAGER=$WRK/bin/ScriptManager-v0.14.jar
SCRIPTMANAGER=$WRK/bin/ScriptManager-v0.14-$PBS_ARRAYID.jar
cp $ORIGINAL_SCRIPTMANAGER $SCRIPTMANAGER

# Script shortcuts
COMPOSITE=$WRK/bin/sum_Col_CDT.pl

# Set up output directories
[ -d logs ] || mkdir logs
[ -d $OUTDIR ] || mkdir $OUTDIR

# Determine BED file for the current job array index
BEDFILE=${BED_LIST[SLURM_ARRAY_TASK_ID]}
BED=`basename $BAMFILE ".bed"`

# Determine BAM file for the current job array index
BAMFILE=${BAM_LIST[SLURM_ARRAY_TASK_ID]}
BAM=`basename $BAMFILE ".bam"`
[ -f $BAMFILE.bai ] || samtools index $BAMFILE
  
echo $BEDFILE "x" $BAMFILE

DIR=$OUTDIR/$BED
[ -d $DIR ] || mkdir $DIR
[[ -d $DIR/CDT ]] || mkdir $DIR/CDT
[[ -d $DIR/COMPOSITES ]] || mkdir $DIR/COMPOSITES
[[ -d $DIR/NormComposites ]] || mkdir $DIR/NormComposites

FACTOR=`grep 'Scaling factor' $WRK/data/NormalizationFactors/$BAM\_NCISb_ScalingFactors.out | awk -F" " '{print $3}'`

### READ 1

echo "Run Custom read1 pileup (filter fragment>100bp)"
BASE=$BAM\_$BED\_read1_MIN100

# Pileup (read1_MIN)
java -jar $SCRIPTMANAGER read-analysis tag-pileup $BEDFILE $BAMFILE --cpu 4 -5 -1 -n 100 -o $DIR/COMPOSITES/$BASE.out -M $DIR/CDT/$BASE
java -jar $SCRIPTMANAGER read-analysis scale-matrix $DIR/CDT/$BASE\_anti.cdt  -s $FACTOR -o $DIR/CDT/$BASE\_anti_Normalized.cdt
java -jar $SCRIPTMANAGER read-analysis scale-matrix $DIR/CDT/$BASE\_sense.cdt -s $FACTOR -o $DIR/CDT/$BASE\_sense_Normalized.cdt

# Make stranded composites
perl $COMPOSITE $DIR/CDT/$BASE\_anti_Normalized.cdt $DIR/CDT/$BASE\_ANTI
perl $COMPOSITE $DIR/CDT/$BASE\_sense_Normalized.cdt $DIR/CDT/$BASE\_SENSE
cat $DIR/CDT/$BASE\_ANTI <(tail -1 $DIR/CDT/$BASE\_SENSE) > $DIR/NormComposites/$BASE\_Normalized.out
rm $DIR/CDT/$BASE\_SENSE $DIR/CDT/$BASE\_ANTI

### READ 2

echo "Run Custom read2 pileup (filter fragment>100bp)"
BASE=$BAM\_$BED\_read2_MIN100

# Pileup (read2_MIN)
java -jar $SCRIPTMANAGER read-analysis tag-pileup $BEDFILE $BAMFILE --cpu 4 -5 -2 -n 100 -o $DIR/COMPOSITES/$BASE.out -M $DIR/CDT/$BASE
java -jar $SCRIPTMANAGER read-analysis scale-matrix $DIR/CDT/$BASE\_anti.cdt  -s $FACTOR -o $DIR/CDT/$BASE\_anti_Normalized.cdt
java -jar $SCRIPTMANAGER read-analysis scale-matrix $DIR/CDT/$BASE\_sense.cdt -s $FACTOR -o $DIR/CDT/$BASE\_sense_Normalized.cdt

# Make stranded composites
perl $COMPOSITE $DIR/CDT/$BASE\_anti_Normalized.cdt $DIR/CDT/$BASE\_ANTI
perl $COMPOSITE $DIR/CDT/$BASE\_sense_Normalized.cdt $DIR/CDT/$BASE\_SENSE
cat $DIR/CDT/$BASE\_ANTI <(tail -1 $DIR/CDT/$BASE\_SENSE) > $DIR/NormComposites/$BASE\_Normalized.out
rm $DIR/CDT/$BASE\_SENSE $DIR/CDT/$BASE\_ANTI
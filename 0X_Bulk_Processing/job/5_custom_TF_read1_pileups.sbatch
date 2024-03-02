#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 00:40:00
#SBATCH -A open
#SBATCH -o logs/5_custom_TF_read1_pileups.log.out
#SBATCH -e logs/5_custom_TF_read1_pileups.log.err
#SBATCH --array 1-3

# Pileup read1 (exo cut sites) for a custom combinations of BAM x BED files.
# Heatmaps included and all with scaling.

BAM_FILES=(
    "$WRK/../data/BAM/K562_NFIA_BX_rep1_hg19.bam"
    "$WRK/../data/BAM/K562_NFIA_BX_rep1_hg19.bam"
    "$WRK/../data/BAM/K562_NFIA_BX_rep1_hg19.bam"
    "$WRK/../data/BAM/K562_NFIA_BX_rep1_hg19.bam"
    "$WRK/../data/BAM/K562_NFIA_BX_rep1_hg19.bam"
    "$WRK/../data/BAM/K562_NFIA_BX_rep1_hg19.bam" 
    "$WRK/../data/BAM/K562_NFIA_BX_rep1_hg19.bam"
    "$WRK/../data/BAM/K562_NFIA_BX_rep1_hg19.bam"
    "$WRK/../data/BAM/K562_NFIA_BX_rep1_hg19.bam"
#9
    "$WRK/../data/BAM/K562_NFIA_BX_rep2_hg19.bam"
    "$WRK/../data/BAM/K562_NFIA_BX_rep2_hg19.bam"
    "$WRK/../data/BAM/K562_NFIA_BX_rep2_hg19.bam"
    "$WRK/../data/BAM/K562_NFIA_BX_rep2_hg19.bam"
    "$WRK/../data/BAM/K562_NFIA_BX_rep2_hg19.bam"
    "$WRK/../data/BAM/K562_NFIA_BX_rep2_hg19.bam"
    "$WRK/../data/BAM/K562_NFIA_BX_rep2_hg19.bam"
    "$WRK/../data/BAM/K562_NFIA_BX_rep2_hg19.bam"
    "$WRK/../data/BAM/K562_NFIA_BX_rep2_hg19.bam"
#9
    "$WRK/../data/BAM/K562_SP1_BX_rep1_hg19.bam"
    "$WRK/../data/BAM/K562_SP1_BX_rep2_hg19.bam"
    "$WRK/../data/BAM/K562_EP300_BX_rep1_hg19.bam"
    "$WRK/../data/BAM/K562_EP300_BX_rep2_hg19.bam"
    "$WRK/../data/BAM/K562_Pol2_BX_rep1_hg19.bam"
    "$WRK/../data/BAM/K562_Pol2_BX_rep2_hg19.bam"
#6
    "$WRK/../data/BAM/K562_GABPA_BX_rep1_hg19.bam"
    "$WRK/../data/BAM/K562_GABPA_BX_rep2_hg19.bam"
    "$WRK/../data/BAM/K562_EP300_BX_rep1_hg19.bam"
    "$WRK/../data/BAM/K562_EP300_BX_rep2_hg19.bam"
    "$WRK/../data/BAM/K562_Pol2_BX_rep1_hg19.bam"
    "$WRK/../data/BAM/K562_Pol2_BX_rep2_hg19.bam"
#6
    "$WRK/../data/BAM/K562_GABPA_BX_rep1_hg19.bam"
    "$WRK/../data//K562_GABPA_BX_rep2_hg19.bam"
    "$WRK/../data/BAM/K562_EP300_BX_rep1_hg19.bam"
    "$WRK/../data/BAM/K562_EP300_BX_rep2_hg19.bam"
    "$WRK/../data/BAM/K562_Pol2_BX_rep1_hg19.bam"
    "$WRK/../data/BAM/K562_Pol2_BX_rep2_hg19.bam"
    "$WRK/../data/.BAM/K562_SP1_BX_rep1_hg19.bam"
    "$WRK/../data//K562_SP1_BX_rep2_hg19.bam"
#8
)

BED_FILES=(
    "$WRK/../data/RefPT-Motif/1000bp/NFIA_Occupancy_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/scrambledNFIA_top10K_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/NFIA_NucSort_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/NFIA_NucSort-DOWNSTREAM_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/NFIA_NucSort-OVERLAP_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/NFIA_NucSort-UPSTREAM_1000bp.bed"
    "$WRK/../data/RefPT-Motif/100bp/NFIA_NucSort-DOWNSTREAM_100bp.bed"
    "$WRK/../data/RefPT-Motif/100bp/NFIA_NucSort-OVERLAP_100bp.bed"
    "$WRK/../data/RefPT-Motif/100bp/NFIA_NucSort-UPSTREAM_100bp.bed"
#
    "$WRK/../data/RefPT-Motif/1000bp/NFIA_Occupancy_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/scrambledNFIA_top10K_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/NFIA_NucSort_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/NFIA_NucSort-DOWNSTREAM_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/NFIA_NucSort-OVERLAP_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/NFIA_NucSort-UPSTREAM_1000bp.bed"
    "$WRK/../data/RefPT-Motif/100bp/NFIA_NucSort-DOWNSTREAM_100bp.bed"
    "$WRK/../data/RefPT-Motif/100bp/NFIA_NucSort-OVERLAP_100bp.bed"
    "$WRK/../data/RefPT-Motif/100bp/NFIA_NucSort-UPSTREAM_100bp.bed"
#
    "$WRK/../data/RefPT-Motif/1000bp/GABPA_Occupancy_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/GABPA_Occupancy_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/GABPA_Occupancy_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/GABPA_Occupancy_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/GABPA_Occupancy_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/GABPA_Occupancy_1000bp.bed"
#
    "$WRK/../data/RefPT-Motif/1000bp/SP1_Occupancy_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/SP1_Occupancy_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/SP1_Occupancy_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/SP1_Occupancy_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/SP1_Occupancy_1000bp.bed"
    "$WRK/../data/RefPT-Motif/1000bp/SP1_Occupancy_1000bp.bed"
#
    "$WRK/BED/K562_CoPRO-expressed_Gene-refSeqTSS_1000bp.bed"
    "$WRK/BED/K562_CoPRO-expressed_Gene-refSeqTSS_1000bp.bed"
    "$WRK/BED/K562_CoPRO-expressed_Gene-refSeqTSS_1000bp.bed"
    "$WRK/BED/K562_CoPRO-expressed_Gene-refSeqTSS_1000bp.bed"
    "$WRK/BED/K562_CoPRO-expressed_Gene-refSeqTSS_1000bp.bed"
    "$WRK/BED/K562_CoPRO-expressed_Gene-refSeqTSS_1000bp.bed"
    "$WRK/BED/K562_CoPRO-expressed_Gene-refSeqTSS_1000bp.bed"
    "$WRK/BED/K562_CoPRO-expressed_Gene-refSeqTSS_1000bp.bed"
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
[[ -d $DIR/PNG/Strand ]] || mkdir -p $DIR/PNG/Strand
[[ -d $DIR/SVG ]] || mkdir $DIR/SVG

# Get scaling factor
FACTOR=`grep 'Scaling factor' $WRK/data/NormalizationFactors/$BAM\_NCISb_ScalingFactors.out | awk -F" " '{print $3}'`

echo "Run Custom read1 pileup"
BASE=$BAM\_$BED\_read1

# Pileup (read 1)
java -jar $SCRIPTMANAGER read-analysis tag-pileup $BEDFILE $BAMFILE --cpu 4 -5 -1 -o $DIR/COMPOSITES/$BASE.out -M $DIR/CDT/$BASE
java -jar $SCRIPTMANAGER read-analysis scale-matrix $DIR/CDT/$BASE\_anti.cdt  -s $FACTOR -o $DIR/CDT/$BASE\_anti_Normalized.cdt
java -jar $SCRIPTMANAGER read-analysis scale-matrix $DIR/CDT/$BASE\_sense.cdt -s $FACTOR -o $DIR/CDT/$BASE\_sense_Normalized.cdt

# Heatmaps
java -jar $SCRIPTMANAGER figure-generation heatmap -a 5 --blue $DIR/CDT/$BASE\_sense_Normalized.cdt -o $DIR/PNG/Strand/$BASE\_sense_treeview.png
java -jar $SCRIPTMANAGER figure-generation heatmap -a 5 --red  $DIR/CDT/$BASE\_anti_Normalized.cdt  -o $DIR/PNG/Strand/$BASE\_anti_treeview.png
java -jar $SCRIPTMANAGER figure-generation merge-heatmap $DIR/PNG/Strand/$BASE\_sense_treeview.png $DIR/PNG/Strand/$BASE\_anti_treeview.png -o $DIR/PNG/$BASE\_merge.png
	
#Label heatmap
java -jar $SCRIPTMANAGER figure-generation label-heatmap $DIR/PNG/$BASE\_merge.png \
	-f 20 -l -500 -m 0 -r 500 -o $DIR/SVG/$BASE\_label.svg

# Make stranded composites
perl $COMPOSITE $DIR/CDT/$BASE\_anti_Normalized.cdt $DIR/CDT/$BASE\_ANTI
perl $COMPOSITE $DIR/CDT/$BASE\_sense_Normalized.cdt $DIR/CDT/$BASE\_SENSE
cat $DIR/CDT/$BASE\_ANTI <(tail -1 $DIR/CDT/$BASE\_SENSE) > $DIR/NormComposites/$BASE\_Normalized.out
rm $DIR/CDT/$BASE\_SENSE $DIR/CDT/$BASE\_ANTI


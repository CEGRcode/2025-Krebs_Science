#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 00:15:00
#SBATCH -A open
#SBATCH -o logs/4_midpoint_pileups.log.out-%a
#SBATCH -e logs/4_midpoint_pileups.log.err-%a
#SBATCH --array 1-4

# Pileup midpoint (heatmaps) and composites (endo cut sites) for a custom
# combination of BAM x BED files (primarily nucleosome BAM). No scaling.
#  - data/BAM/BNase-seq_50U-10min_merge_hg38.bam

### CHANGE ME
WRK=/path/to/2024-Chen_Nature/X_Bulk_Processing
WRK=/storage/home/owl5022/scratch/2024-Chen_Nature/X_Bulk_Processing
METADATA=Nucleosome_pileups.txt
###

# Dependencies
# - java
# - samtools

set -exo
module load samtools

# Fill in placeholder constants with your directories
BAMDIR=../data/BAM

# Setup ScriptManager for job array
ORIGINAL_SCRIPTMANAGER=../bin/ScriptManager-v0.15.jar
SCRIPTMANAGER=../bin/ScriptManager-v0.15-$SLURM_ARRAY_TASK_ID.jar
cp $ORIGINAL_SCRIPTMANAGER $SCRIPTMANAGER

# Script shortcuts
COMPOSITE=../bin/sum_Col_CDT.pl

# Set up output directories
[ -d logs ] || mkdir logs
[ -d Library ] || mkdir Library

# Determine rgb heatmap color
RGB=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $1}'`

# Determine insert filters
MIN=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $2}'`
MAX=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $3}'`

# Determine heatmap contrast threshold
CONTRAST_THRESH=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $4}'`

# Determine scaling method
NORMALIZATION=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $5}'`

# Determine BAM file for the current job array index
BAMFILE=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $6}'`
BAM=`basename $BAMFILE ".bam"`
[ -f $BAMFILE.bai ] || samtools index $BAMFILE

# Determine BED file for the current job array index
BEDFILE=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $7}'`
BED=`basename $BEDFILE ".bed"`

# Parse min/max params
TP_PARAMS=""
FILTER=""
[[ $MIN =~ "NaN" ]] || FILTER=-MIN${MIN}
[[ $MIN =~ "NaN" ]] || TP_PARAMS="-n $MIN"
[[ $MAX =~ "NaN" ]] || FILTER=${FILTER}-MAX${MAX}
[[ $MAX =~ "NaN" ]] || TP_PARAMS=$TP_PARAMS" -x $MAX"

# Parse heatmap params
HM_PARAMS="-p 0.95"
[[ $CONTRAST_THRESH =~ "NaN" ]] || HM_PARAMS="-a $CONTRAST_THRESH"
[[ $RGB =~ "NaN" ]] || HM_PARAMS="$HM_PARAMS --color $RGB"

# Parse normalization params and get scaling factor
FACTOR=1
if [ "$NORMALIZATION" == "NCIS" ]; then
    NFFILE=$BAMDIR/NormalizationFactors/$BAM\_NCISb_ScalingFactors.out
    [[ -f $NFFILE ]] || exit "Missing Normalization factor file (NCIS): $NFFILE"
    FACTOR=`grep 'Scaling factor' $NFFILE | awk -F" " '{print $3}'`
elif [ "$NORMALIZATION" == "TotalTag" ]; then
    NFFILE=$BAMDIR/NormalizationFactors/$BAM\_TotalTag_ScalingFactors.out
    [[ -f $NFFILE ]] || exit "Missing Normalization factor file (TotalTag): $NFFILE"
    FACTOR=`grep 'Scaling factor' $NFFILE | awk -F" " '{print $3}'`
else
    NORMALIZATION="Raw"
fi

echo "(${SLURM_ARRAY_TASK_ID}) ${BEDFILE} x ${BAMFILE} x ${MIN} x ${MAX} x ${NORMALIZATION}"
BASE=${BAM}_${BED}_midpoint${FILTER}

# Setup output directories
DIR=Library/$BED
[ -d $DIR ] || mkdir $DIR
[[ -d $DIR/CDT ]] || mkdir $DIR/CDT
[[ -d $DIR/Composites ]] || mkdir $DIR/Composites
[[ -d $DIR/PNG/Strand ]] || mkdir -p $DIR/PNG/Strand
[[ -d $DIR/SVG ]] || mkdir $DIR/SVG


# ===============================================================================================================================

echo "Run Custom midpoint pileup"

# Pileup (mid)
java -jar $SCRIPTMANAGER read-analysis tag-pileup $BEDFILE $BAMFILE --cpu 4 -m --combined -o $DIR/Composites/${BASE}_combined.out -M $DIR/CDT/$BASE

# Scale matrix
java -jar $SCRIPTMANAGER read-analysis scale-matrix $DIR/Composites/${BASE}_combined.out -l 1 -s $FACTOR -o $DIR/Composites/${BASE}_${NORMALIZATION}_combined.out
java -jar $SCRIPTMANAGER read-analysis scale-matrix $DIR/CDT/${BASE}_combined.cdt    -s $FACTOR -o $DIR/CDT/${BASE}_${NORMALIZATION}_combined.cdt

BASE=${BASE}_${NORMALIZATION}

# Two-color heatmap
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation heatmap -p .95 --black $DIR/CDT/${BASE}_combined.cdt -o $DIR/PNG/${BASE}_combined.png

# Count sites
NSITES=`wc -l $BEDFILE | awk '{print $1-1}'`

# Label heatmap based on BED naming
if [[ "$BED" == *"1000bp"* ]]; then  
    java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation label-heatmap $DIR/PNG/${BASE}_combined.png \
        -l "-500" -m "0" -r "+500" -w 1 -f 20 \
        -o "$DIR/SVG/${BASE}_combined.svg"

elif [[ "$BED" == *"500bp"* ]]; then  
    java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation label-heatmap $DIR/PNG/${BASE}_combined.png \
        -l "-250" -m "0" -r "+250" -w 1 -f 20 \
        -o "$DIR/SVG/${BASE}_combined.svg"

elif [[ "$BED" == *"250bp"* ]]; then  
    java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation label-heatmap $DIR/PNG/${BASE}_combined.png \
        -l "0" -m "125" -r "+250" -w 1 -f 20 \
        -o "$DIR/SVG/${BASE}_combined.svg"
fi

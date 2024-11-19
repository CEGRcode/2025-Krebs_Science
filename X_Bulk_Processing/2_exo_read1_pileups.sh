#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=24gb
#SBATCH -t 00:40:00
#SBATCH -A open
#SBATCH -o logs/2_exo_read1_pileups.log.out-%a
#SBATCH -e logs/2_exo_read1_pileups.log.err-%a
#SBATCH --array 1-37

# Pileup read1 (exo cut sites) for a custom combinations of BAM x BED files.
# Heatmaps included and all with scaling. 

### CHANGE ME
WRK=/path/to/2024-Krebs_Science/X_Bulk_Processing
WRK=/storage/home/owl5022/scratch/2024-Krebs_Science/X_Bulk_Processing
METADATA=Read1_pileups.txt
###

# Dependencies
# - java
# - perl
# - samtools

set -exo
module load samtools

# Fill in placeholder constants with your directories
BAMDIR=../data/BAM

# Setup ScriptManager for job array
#SCRIPTMANAGER=$WRK/bin/ScriptManager-v0.14.jar
ORIGINAL_SCRIPTMANAGER=../bin/ScriptManager-v0.15.jar
SCRIPTMANAGER=../bin/ScriptManager-v0.15-$SLURM_ARRAY_TASK_ID.jar
cp $ORIGINAL_SCRIPTMANAGER $SCRIPTMANAGER

# Script shortcuts
COMPOSITE=../bin/sum_Col_CDT.pl

# Set up output directories
[ -d logs ] || mkdir logs
[ -d Library ] || mkdir Library

# Determine insert filters
MIN=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $1}'`
MAX=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $2}'`

# Determine heatmap contrast threshold
CONTRAST_THRESH=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $3}'`

# Determine scaling method
NORMALIZATION=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $4}'`

# Determine BAM file for the current job array index
BAMFILE=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $5}'`
BAM=`basename $BAMFILE ".bam"`
[ -f $BAMFILE.bai ] || samtools index $BAMFILE

# Determine BED file for the current job array index
BEDFILE=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $6}'`
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
BASE=${BAM}_${BED}_5read1${FILTER}

# Setup output directories
DIR=Library/$BED
[ -d $DIR ] || mkdir $DIR
[[ -d $DIR/CDT ]] || mkdir $DIR/CDT
[[ -d $DIR/Composites ]] || mkdir $DIR/Composites
[[ -d $DIR/PNG/Strand ]] || mkdir -p $DIR/PNG/Strand
[[ -d $DIR/SVG ]] || mkdir $DIR/SVG

# ===============================================================================================================================

echo "Run Custom read1 pileup"

# Pileup (read 1)
java -jar $SCRIPTMANAGER read-analysis tag-pileup $BEDFILE $BAMFILE $TP_PARAMS --cpu 4 -5 -1 -o $DIR/Composites/$BASE.out -M $DIR/CDT/$BASE

# Scale matrix
java -jar $SCRIPTMANAGER read-analysis scale-matrix $DIR/Composites/$BASE.out  -s $FACTOR -o $DIR/Composites/${BASE}_${NORMALIZATION}.out -l 1
java -jar $SCRIPTMANAGER read-analysis scale-matrix $DIR/CDT/${BASE}_anti.cdt  -s $FACTOR -o $DIR/CDT/${BASE}_${NORMALIZATION}_anti.cdt
java -jar $SCRIPTMANAGER read-analysis scale-matrix $DIR/CDT/${BASE}_sense.cdt -s $FACTOR -o $DIR/CDT/${BASE}_${NORMALIZATION}_sense.cdt

BASE=${BASE}_${NORMALIZATION}

# Two-color heatmap & merge
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation heatmap $HM_PARAMS --blue $DIR/CDT/${BASE}_sense.cdt -o $DIR/PNG/Strand/${BASE}_sense_treeview.png
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation heatmap $HM_PARAMS --red  $DIR/CDT/${BASE}_anti.cdt  -o $DIR/PNG/Strand/${BASE}_anti_treeview.png
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation merge-heatmap $DIR/PNG/Strand/${BASE}_sense_treeview.png $DIR/PNG/Strand/${BASE}_anti_treeview.png -o $DIR/PNG/${BASE}_merge.png

# Count sites
NSITES=`wc -l $BEDFILE | awk '{print $1-1}'`

# Label heatmap based on BED naming
if [[ "$BED" == *"1000bp"* ]]; then
    java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation label-heatmap $DIR/PNG/${BASE}_merge.png \
        -l "-500" -m "0" -r "+500" -w 1 -f 20 \
        -o "$DIR/SVG/${BASE}_merge_label.svg"

elif [[ "$BED" == *"500bp"* ]]; then
    java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation label-heatmap $DIR/PNG/${BASE}_merge.png \
        -l "-250" -m "0" -r "+250" -w 1 -f 20 \
        -o "$DIR/SVG/${BASE}_merge_label.svg"

elif [[ "$BED" == *"250bp"* ]]; then
    java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation label-heatmap $DIR/PNG/${BASE}_merge.png \
        -l "0" -m "125" -r "+250" -w 1 -f 20 \
        -o "$DIR/SVG/${BASE}_merge_label.svg"
fi

# Make stranded composites
perl $COMPOSITE $DIR/CDT/$BASE\_anti.cdt $DIR/CDT/$BASE\_ANTI
perl $COMPOSITE $DIR/CDT/$BASE\_sense.cdt $DIR/CDT/$BASE\_SENSE
cat $DIR/CDT/$BASE\_ANTI <(tail -1 $DIR/CDT/$BASE\_SENSE) > $DIR/Composites/$BASE.out
rm $DIR/CDT/$BASE\_SENSE $DIR/CDT/$BASE\_ANTI

#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=24gb
#SBATCH -t 00:40:00
#SBATCH -A open
#SBATCH -o logs/5_Five_Read2_Pileups.log.out-%a
#SBATCH -e logs/5_Five_Read2_Pileups.log.err-%a
#SBATCH --array 1-26

# Pileup heatmaps and composites (5' read2, endo cut sites) for a custom
# combination of BAM x BED files (primarily nucleosome BAM). Configuration
# file specifies insert size filters, contrast threshold, RGB heatmap color,
# and normalization strategy.

### CHANGE ME
WRK=/path/to/2024-Krebs_Science/X_Bulk_Processing
WRK=/ocean/projects/see180003p/owlang/2024-Krebs_Science/X_Bulk_Processing
WRK=/storage/group/bfp2/default/owl5022-OliviaLang/2024-Krebs_Science/X_Bulk_Processing
WRK=/scratch/owl5022/2024-Krebs_Science/X_Bulk_Processing
METADATA=Five_Read2_Pileups.txt
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

# Determine insert filters for the current job array index
MIN=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $1}'`
MAX=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $2}'`

# Determine heatmap contrast threshold for the current job array index
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

# ===============================================================================================================================

# Count sites
NSITES=`wc -l $BEDFILE | awk '{print $1}'`

# Set window labels
RLABEL="R"
LLABEL="L"
WINDOW=`echo $BED | rev | cut -d"_" -f1 | rev`
[ $WINDOW == "2000bp" ] && LLABEL="-1" && RLABEL="+1"
[ $WINDOW == "1000bp" ] && LLABEL="-500" && RLABEL="+500"
[ $WINDOW == "500bp" ] && LLABEL="-250" && RLABEL="+250"
[ $WINDOW == "250bp" ] && LLABEL="-125" && RLABEL="+125"
[ $WINDOW == "200bp" ] && LLABEL="-100" && RLABEL="+100"

# Parse min/max params
TP_PARAMS="--cpu $THREADS -5 -2"
FILTER=""
[[ $MIN =~ "NaN" ]] || FILTER=${FILTER}-MIN${MIN}
[[ $MIN =~ "NaN" ]] || TP_PARAMS="${TP_PARAMS} -n ${MIN}"
[[ $MAX =~ "NaN" ]] || FILTER=${FILTER}-MAX${MAX}
[[ $MAX =~ "NaN" ]] || TP_PARAMS="${TP_PARAMS} -x ${MAX}"

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
BASE=${BAM}_${BED}_5read2${FILTER}

# Setup output directories
DIR=Library/$BED
[ -d $DIR ] || mkdir $DIR
[[ -d $DIR/CDT ]] || mkdir $DIR/CDT
[[ -d $DIR/Composites ]] || mkdir $DIR/Composites
[[ -d $DIR/PNG/Strand ]] || mkdir -p $DIR/PNG/Strand
[[ -d $DIR/SVG ]] || mkdir $DIR/SVG

# ===============================================================================================================================

echo "Run Custom read2 pileup"

# Pileup (read 1)
java -jar $SCRIPTMANAGER read-analysis tag-pileup $BEDFILE $BAMFILE $TP_PARAMS -o $DIR/Composites/$BASE.out -M $DIR/CDT/$BASE

# Scale matrix
java -jar $SCRIPTMANAGER read-analysis scale-matrix $DIR/Composites/$BASE.out  -s $FACTOR -o $DIR/Composites/${BASE}_${NORMALIZATION}.out -l 1
java -jar $SCRIPTMANAGER read-analysis scale-matrix $DIR/CDT/${BASE}_anti.cdt  -s $FACTOR -o $DIR/CDT/${BASE}_${NORMALIZATION}_anti.cdt
java -jar $SCRIPTMANAGER read-analysis scale-matrix $DIR/CDT/${BASE}_sense.cdt -s $FACTOR -o $DIR/CDT/${BASE}_${NORMALIZATION}_sense.cdt

BASE=${BASE}_${NORMALIZATION}

# Two-color heatmap & merge
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation heatmap $HM_PARAMS -c FF00FF $DIR/CDT/${BASE}_sense.cdt -o $DIR/PNG/Strand/${BASE}_sense_treeview.png
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation heatmap $HM_PARAMS -c 00BFFF $DIR/CDT/${BASE}_anti.cdt  -o $DIR/PNG/Strand/${BASE}_anti_treeview.png
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation merge-heatmap $DIR/PNG/Strand/${BASE}_sense_treeview.png $DIR/PNG/Strand/${BASE}_anti_treeview.png -o $DIR/PNG/${BASE}_merge.png

# Label heatmap
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation label-heatmap $DIR/PNG/${BASE}_merge.png \
	-l $LLABEL -m "0" -r $RLABEL -w 1 -f 20 \
	-x $BED -y "$BED occurences (${NSITES} sites)" \
	-o $DIR/SVG/${BASE}_merge_label.svg

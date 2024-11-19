#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 00:15:00
#SBATCH -A open
#SBATCH -o logs/6_Conservation-phylo_pileups.log.out-%a
#SBATCH -e logs/6_Conservation-phylo_pileups.log.err-%a
#SBATCH --array 1

### CHANGE ME
WRK=/path/to/2024-Krebs_Science/X_Bulk_Processing
WRK=/storage/home/owl5022/scratch/2024-Krebs_Science/X_Bulk_Processing
###

# Dependencies
# - java
# - perl
# - python
# - pyBigWig

set -exo
module load anaconda
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Fill in placeholder constants with your directories
MOTIF=../data/RefPT-Motif
CDIR=../data/Conservation-SNP

# Script shortcuts
SCRIPTMANAGER=../bin/ScriptManager-v0.15.jar
COMPOSITE=../bin/sum_Col_CDT.pl
WIGTOBG=../bin/convert_wig_to_bedgraph.py
PILEUPBG=../bin/pileup_BedGraph_on_RefPT.py
PILEUPBW=../bin/pileup_BigWig_on_RefPT.py
SPILEUPBW=../bin/pileup_BigWig_on_RefPT_stranded.py
SUMMAT=../bin/sum_each_CDT.py

# Set up output directories
[ -d logs ] || mkdir logs
[ -d Library ] || mkdir Library

# Define set of BED files to pileup on
LIST=(
    "$MOTIF/500bp/NFIA_SORT-Occupancy_500bp.bed"
)

# Select BED file corresponding to this job index
INDEX=$(($SLURM_ARRAY_TASK_ID-1))
BEDFILE=${LIST[$INDEX]}
BED=`basename $BEDFILE ".bed"`

DIR=Library/$BED
[ -d $DIR ] || mkdir $DIR
[ -d $DIR/CDT ] || mkdir $DIR/CDT
[ -d $DIR/Composites ] || mkdir $DIR/Composites
[ -d $DIR/PNG ] || mkdir $DIR/PNG
[ -d $DIR/SVG ] || mkdir $DIR/SVG

# Pileup conservation scores
CONSERVATION=$CDIR/hg38.phyloP30way.bw
CONS=`basename $CONSERVATION ".bw"`

echo $BED x $CONS

# Pileup BigWig
python $PILEUPBW -i $CONSERVATION -r $BEDFILE -o $DIR/CDT/${CONS}_$BED.cdt
# Make composite
perl $COMPOSITE $DIR/CDT/${CONS}_$BED.cdt $DIR/Composites/${CONS}_$BED.out

# Count sites
NSITES=`wc -l $BEDFILE | awk '{print $1-1}'`

# Two-color heatmap
java -jar $SCRIPTMANAGER figure-generation heatmap -p 0.95 --black $DIR/CDT/${CONS}_$BED.cdt -o $DIR/PNG/${CONS}_$BED.png
java -jar $SCRIPTMANAGER figure-generation label-heatmap $DIR/PNG/${CONS}_$BED.png \
	-l "-250" -m "0" -r "+250" -w 1 -f 20 \
	-x $BED -y "$BED occurences (${NSITES} sites)" \
	-o $DIR/SVG/${CONS}_$BED.svg

#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 00:15:00
#SBATCH -A open
#SBATCH -o logs/8_dbSnp153_pileups.log.out-%a
#SBATCH -e logs/8_dbSnp153_pileups.log.err-%a
#SBATCH --array 1

### CHANGE ME
WRK=/path/to/2024-Krebs_Science/X_Bulk_Processing
WRK=/scratch/owl5022/2024-Krebs_Science/X_Bulk_Processing
###

# Dependencies
# - java
# - perl
# - python
# - samtools

set -exo
module load anaconda
module load samtools
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

# Count sites
DIR=Library/$BED
[ -d $DIR ] || mkdir $DIR
[ -d $DIR/CDT ] || mkdir $DIR/CDT
[ -d $DIR/Composites ] || mkdir $DIR/Composites

# Pileup SNPs
DBSNP=$CDIR/dbSnp153_snv.bw
SNP=`basename $DBSNP ".bw"`

echo $BED x $SNP

# Pileup BigWig
python $PILEUPBW -i $DBSNP -r $BEDFILE -o $DIR/CDT/${SNP}_$BED.cdt
# Make composite
perl $COMPOSITE $DIR/CDT/${SNP}_$BED.cdt $DIR/Composites/${SNP}_$BED.out

# Two-color heatmap
java -jar $SCRIPTMANAGER figure-generation heatmap -p 0.95 --black $DIR/CDT/${SNP}_$BED.cdt -o $DIR/PNG/${SNP}_$BED.png
java -jar $SCRIPTMANAGER figure-generation label-heatmap $DIR/PNG/${SNP}_$BED.png \
	-l "-500" -m "0" -r "+500" -w 1 -f 20 \
	-x $BED -y "$BED occurences (${NSITES} sites)" \
	-o $DIR/SVG/${SNP}_$BED.svg

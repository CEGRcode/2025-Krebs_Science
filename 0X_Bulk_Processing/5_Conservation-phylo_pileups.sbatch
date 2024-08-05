#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 02:00:00
#SBATCH -A open
#SBATCH -o logs/5_Conservation-phylo_pileups.log.out-%a
#SBATCH -e logs/5_Conservation-phylo_pileups.log.err-%a
#SBATCH --array 1-4

# Pileup read1 and read2 (exo & endo cut sites) for custom combinations of BAM x BED files
# that are filtered by fragment read size (>100bp). Just composites with scaling.

### CHANGE ME
WRK=/path/to/2024-Chen_Nature
OUTDIR=$WRK/Library
#WRK=/storage/home/owl5022/scratch/2024-Chen_Nature/0X_Bulk_Processing
###

# Dependencies
# - java
# - perl
# - python

set -exo
module load anaconda
source activate bx

# Fill in placeholder constants with your directories
MOTIF=$WRK/data/RefPT-Motif/1000bp
CDIR=$WRK/data/Conservation-SNP


# Setup ScriptManager for job array
ORIGINAL_SCRIPTMANAGER=$WRK/bin/ScriptManager-v0.14.jar
SCRIPTMANAGER=$WRK/bin/ScriptManager-v0.14-$SLURM_ARRAY_TASK_ID.jar
cp $ORIGINAL_SCRIPTMANAGER $SCRIPTMANAGER

# Script shortcuts
SCRIPTMANAGER=$WRK/bin/scriptmanager-v0.14-dev.jar
COMPOSITE=$WRK/bin/sum_Col_CDT.pl
# WIGTOBG=$WRK/../bin/convert_wig_to_bedgraph.py
# PILEUPBG=$WRK/../bin/pileup_BedGraph_on_RefPT.py
PILEUPBW=$WRK/bin/pileup_BigWig_on_RefPT.py
# SPILEUPBW=$WRK/../bin/pileup_BigWig_on_RefPT_stranded.py
# SUMMAT=$WRK/../bin/sum_each_CDT.py

# Set up output directories
[ -d logs ] || mkdir logs
[ -d $OUTDIR ] || mkdir $OUTDIR

# Define set of BED files to pileup on
BED_LIST=(
	"$MOTIF/CTCF_Occupancy_1000bp.bed"
	"$MOTIF/NFIA_engageNucdown_NucSort_1000bp.bed"
	"$MOTIF/scrambledCTCF_NucSort_1000bp.bed"
	"$MOTIF/scrambledNFIA_NucSort_1000bp.bed"
)
# Determine BED file for the current job array index
INDEX=$(($SLURM_ARRAY_TASK_ID-1))
BEDFILE=${BED_LIST[$INDEX]}
BED=`basename $BEDFILE ".bed"`

[ -d $OUTDIR/$BED ] || mkdir $OUTDIR/$BED
[ -d $OUTDIR/$BED/CDT ] || mkdir $OUTDIR/$BED/CDT
[ -d $OUTDIR/$BED/Composites ] || mkdir $OUTDIR/$BED/Composites
[ -d $OUTDIR/$BED/PNG ] || mkdir $OUTDIR/$BED/PNG
[ -d $OUTDIR/$BED/SVG ] || mkdir $OUTDIR/$BED/SVG

# Pileup conservation scores
for CONSERVATION in "$CDIR/placentalMammals.phyloP46way.bw" "$CDIR/primates.phyloP46way.bw" "$CDIR/vertebrate.phyloP46way.bw";
do
	CONS=`basename $CONSERVATION ".bw"`

	echo $BED x $CONS

	# Pileup BigWig
	python $PILEUPBW -i $CONSERVATION -r $BEDFILE -o $OUTDIR/$BED/CDT/$CONS\_$BED.cdt
	# Make composite
	perl $COMPOSITE $OUTDIR/$BED/CDT/$CONS\_$BED.cdt $OUTDIR/$BED/Composites/$CONS\_$BED.out

	# Count sites
	NSITES=`wc -l $BEDFILE | awk '{print $1-1}'`

	# Two-color heatmap
	java -jar $SCRIPTMANAGER figure-generation heatmap -p 0.95 --black $OUTDIR/$BED/CDT/$CONS\_$BED.cdt -o $OUTDIR/$BED/PNG/$CONS\_$BED.png
	java -jar $SCRIPTMANAGER figure-generation label-heatmap $OUTDIR/$BED/PNG/$CONS\_$BED.png \
		-l "-500" -m "0" -r "+500" -w 1 -f 20 \
		#-x $BED -y "$BED occurences (${NSITES} sites)" \
		-o $OUTDIR/$BED/SVG/$CONS\_$BED.svg
done


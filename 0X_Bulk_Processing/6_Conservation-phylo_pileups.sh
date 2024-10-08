
### CHANGE ME
WRK=/storage/group/bfp2/default/hxc585_HainingChen/Fox_NFIA_CTCF/
OUTDIR=$WRK/Library
###

# Dependencies
# - java
# - perl
# - python

set -exo
module load anaconda
#source activate bx

# Fill in placeholder constants with your directories
MOTIF=$WRK/data/RefPT-Motif/1000bp
CDIR=$WRK/data/Conservation-SNP
#pip install pyBigWig   

# Script shortcuts
SCRIPTMANAGER=$WRK/bin/scriptmanager-v0.15-dev.jar
COMPOSITE=$WRK/bin/sum_Col_CDT.pl
WIGTOBG=$WRK/bin/convert_wig_to_bedgraph.py
PILEUPBG=$WRK/bin/pileup_BedGraph_on_RefPT.py
PILEUPBW=$WRK/bin/pileup_BigWig_on_RefPT.py
SPILEUPBW=$WRK/bin/pileup_BigWig_on_RefPT_stranded.py
SUMMAT=$WRK/bin/sum_each_CDT.py

# Set up output directories
[ -d logs ] || mkdir logs
[ -d $OUTDIR ] || mkdir $OUTDIR

# Define set of BED files to pileup on
BEDFILE="$MOTIF/NFIA_downNuc_500bp.bed"
BED=`basename $BEDFILE ".bed"`

[ -d $OUTDIR/$BED ] || mkdir $OUTDIR/$BED
[ -d $OUTDIR/$BED/CDT ] || mkdir $OUTDIR/$BED/CDT
[ -d $OUTDIR/$BED/Composites ] || mkdir $OUTDIR/$BED/Composites
[ -d $OUTDIR/$BED/PNG ] || mkdir $OUTDIR/$BED/PNG
[ -d $OUTDIR/$BED/SVG ] || mkdir $OUTDIR/$BED/SVG

# Pileup conservation scores
for CONSERVATION in "$CDIR/hg38.phyloP30way.wig" ;
do
	CONS=`basename $CONSERVATION ".wig"`

	echo $BED x $CONS

	# Pileup BigWig
	python $PILEUPBW -i $CONSERVATION -r $BEDFILE -o $OUTDIR/$BED/CDT/$CONS\_$BED.cdt
	# Make composite
	perl $COMPOSITE $OUTDIR/$BED/CDT/$CONS\_$BED.cdt $OUTDIR/$BED/Composites/$CONS\_$BED.out

	# Count sites
	#NSITES=`wc -l $BEDFILE | awk '{print $1-1}'`

	# Two-color heatmap
	#java -jar $SCRIPTMANAGER figure-generation heatmap -p 0.95 --black $OUTDIR/$BED/CDT/$CONS\_$BED.cdt -o $OUTDIR/$BED/PNG/$CONS\_$BED.png
	#java -jar $SCRIPTMANAGER figure-generation label-heatmap $OUTDIR/$BED/PNG/$CONS\_$BED.png \
		#-l "-250" -m "0" -r "+250" -w 1 -f 20 \
		#-x $BED -y "$BED occurences (${NSITES} sites)" \
		#-o $OUTDIR/$BED/SVG/$CONS\_$BED.svg
done


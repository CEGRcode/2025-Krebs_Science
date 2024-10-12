
### CHANGE ME
WRK=/Path/to/Title
OUTDIR=$WRK/Library
###

# Dependencies
# - java
# - perl
# - python
# - samtools

set -exo
module load anaconda
module load bedtools
module load samtools
#source activate bx

# Fill in placeholder constants with your directories
MOTIF=$WRK/data/RefPT-Motif/1000bp
CDIR=$WRK/data/Conservation-SNP


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

# Count sites

[ -d $OUTDIR/$BED ] || mkdir $OUTDIR/$BED
[ -d $OUTDIR/$BED/CDT ] || mkdir $OUTDIR/$BED/CDT
[ -d $OUTDIR/$BED/Composites ] || mkdir $OUTDIR/$BED/Composites
#[ -d $OUTDIR/$BED/PNG/Strand ] || mkdir -p $OUTDIR/$BED/PNG/Strand
#[ -d $OUTDIR/$BED/SVG ] || mkdir $OUTDIR/$BED/SVG

# Pileup SNPs
for DBSNP in $CDIR/dbSnp153_snv.bw;
do
	SNP=`basename $DBSNP ".bw"`
	[[ "$SNP" == "dbSnp153" ]] && continue

	echo $BED x $SNP

	python $PILEUPBW -i $DBSNP -r $BEDFILE -o $OUTDIR/$BED/CDT/$SNP\_$BED.cdt
	perl $COMPOSITE $OUTDIR/$BED/CDT/$SNP\_$BED.cdt $OUTDIR/$BED/Composites/$SNP\_$BED.out

	# Two-color heatmap
	#java -jar $SCRIPTMANAGER figure-generation heatmap -p 0.95 --black $OUTDIR/$BED/CDT/$SNP\_$BED.cdt -o $OUTDIR/$BED/PNG/$SNP\_$BED.png
	#java -jar $SCRIPTMANAGER figure-generation label-heatmap $OUTDIR/$BED/PNG/$SNP\_$BED.png \
		#-l "-500" -m "0" -r "+500" -w 1 -f 20 \
		#-x $BED -y "$BED occurences (${NSITES} sites)" \
		#-o $OUTDIR/$BED/SVG/$SNP\_$BED.svg

	#python $SPILEUPBW -i $DBSNP -r $BEDFILE -o $OUTDIR/$BED/CDT/$SNP\_$BED
done

#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 02:00:00
#SBATCH -A open
#SBATCH -o logs/6_conservation_and_SNP_pileups.log.out-%a
#SBATCH -e logs/6_conservation_and_SNP_pileups.log.err-%a
#SBATCH --array 1-12

# Pileup read1 and read2 (exo & endo cut sites) for custom combinations of BAM x BED files
# that are filtered by fragment read size (>100bp). Just composites with scaling.

### CHANGE ME
WRK=/path/to/2024-Chen_Nature/0X_Bulk_Processing
WRK=/storage/home/owl5022/scratch/2024-Chen_Nature/0X_Bulk_Processing
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
source activate bx

# Fill in placeholder constants with your directories
MOTIF=$WRK/../data/RefPT-Motif/1000bp
CDIR=$WRK/../data/Conservation-SNP
OUTDIR=$WRK/Library

# Setup ScriptManager for job array
ORIGINAL_SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.14.jar
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.14-$SLURM_ARRAY_TASK_ID.jar
cp $ORIGINAL_SCRIPTMANAGER $SCRIPTMANAGER

# Script shortcuts
SCRIPTMANAGER=$WRK/../bin/scriptmanager-v0.14-dev.jar
COMPOSITE=$WRK/../bin/sum_Col_CDT.pl
WIGTOBG=$WRK/../bin/convert_wig_to_bedgraph.py
PILEUPBG=$WRK/../bin/pileup_BedGraph_on_RefPT.py
PILEUPBW=$WRK/../bin/pileup_BigWig_on_RefPT.py
SPILEUPBW=$WRK/../bin/pileup_BigWig_on_RefPT_stranded.py
SUMMAT=$WRK/../bin/sum_each_CDT.py


# Set up output directories
[ -d logs ] || mkdir logs
[ -d $OUTDIR ] || mkdir $OUTDIR

# Define set of BED files to pileup on
BED_LIST=(
	"$MOTIF/CTCF_Occupancy_1000bp.bed"
	"$MOTIF/GABPA_Occupancy_1000bp.bed"
	"$MOTIF/NFIA_Occupancy_1000bp.bed"
	"$MOTIF/SP1_Occupancy_1000bp.bed"
	"$MOTIF/WDR5_Occupancy_1000bp.bed"
	"$MOTIF/scrambledCTCF_NucSort_1000bp.bed"
	"$MOTIF/scrambledNFIA_NucSort_1000bp.bed"
	"$MOTIF/scrambledWDR5_NucSort_1000bp.bed"
	"$MOTIF/NFIA_Occupancy-TOP_1000bp.bed"
	"$MOTIF/NFIA_Occupancy-MIDDLE_1000bp.bed"
	"$MOTIF/NFIA_Occupancy-BOTTOM_1000bp.bed"
)

# Determine BED file for the current job array index
INDEX=$(($SLURM_ARRAY_TASK_ID-1))
BEDFILE=${BED_LIST[$INDEX]}
BED=`basename $BEDFILE ".bed"`
[ $OUTDIR/$BED/CDT ] || mkdir -p $OUTDIR/$BED/CDT


# Pileup SNPs
#for DBSNP in "$CDIR/mini/dbSnp153_delins.bed" "$CDIR/mini/dbSnp153_ins.bed" "$CDIR/mini/dbSnp153_mnv.bed" "$CDIR/mini/dbSnp153_refA.bed" "$CDIR/mini/dbSnp153_refC.bed" "$CDIR/mini/dbSnp153_refG.bed" "$CDIR/mini/dbSnp153_refT.bed" "$CDIR/mini/dbSnp153_snv_alt-A.bed" "$CDIR/mini/dbSnp153_snv_alt-C.bed" "$CDIR/mini/dbSnp153_snv_alt-G.bed" "$CDIR/mini/dbSnp153_snv_alt-T.bed" "$CDIR/mini/dbSnp153_snv_type-commonAll.bed" "$CDIR/mini/dbSnp153_snv_type-commonSome.bed" "$CDIR/mini/dbSnp153_snv_type-diffMajor.bed" "$CDIR/mini/dbSnp153_snv_type-overlapDiffClass.bed" "$CDIR/mini/dbSnp153_snv_type-overlapSameClass.bed" "$CDIR/mini/dbSnp153_snv_type-rareAll.bed" "$CDIR/mini/dbSnp153_snv_type-rareSome.bed" "$CDIR/mini/dbSnp153_snv_type-rareAll.bed" "$CDIR/mini/dbSnp153_snv_type-refIsMinor.bed" "$CDIR/mini/dbSnp153_snv_type-refIsRare.bed" "$CDIR/mini/dbSnp153_snv_type-revStrand" "$CDIR/mini/dbSnp153_snv_type-clinvar";

for DBSNP in $CDIR/mini/dbSnp153_snv_*.bw;
do
	SNP=`basename $DBSNP ".bw"`
	[[ "$SNP" == "dbSnp153" ]] && continue

	echo $BED x $SNP

	python $PILEUPBW -i $DBSNP -r $BEDFILE -o $OUTDIR/$BED/CDT/$SNP\_$BED.cdt
	python $SPILEUPBW -i $DBSNP -r $BEDFILE -o $OUTDIR/$BED/CDT/$SNP\_$BED
	perl $COMPOSITE $OUTDIR/$BED/CDT/$SNP\_$BED.cdt $OUTDIR/$BED/Composites/$SNP\_$BED.out
	gzip $OUTDIR/$BED/CDT/$SNP\_$BED.cdt

	java -jar $SCRIPTMANAGER figure-generation heatmap -p 0.95 --black $OUTDIR/$BED/CDT/$SNP\_$BED.cdt.gz -o $OUTDIR/$BED/PNG/$SNP\_$BED.png
	java -jar $SCRIPTMANAGER figure-generation label-heatmap $OUTDIR/$BED/PNG/$SNP\_$BED.png \
		-f 20 -l -500 -m 0 -r 500 -o $OUTDIR/$BED/SVG/$SNP\_$BED.svg

done

# Merge alt with strand correction
python $SUMMAT -1 $OUTDIR/$BED/CDT/dbSnp153_snv_alt-A_$BED\_same.cdt -2 $OUTDIR/$BED/CDT/dbSnp153_snv_alt-T_$BED\_opposite.cdt -o $OUTDIR/$BED/CDT/dbSnp153_snv_alt-A_$BED\_merged.cdt
python $SUMMAT -1 $OUTDIR/$BED/CDT/dbSnp153_snv_alt-T_$BED\_same.cdt -2 $OUTDIR/$BED/CDT/dbSnp153_snv_alt-A_$BED\_opposite.cdt -o $OUTDIR/$BED/CDT/dbSnp153_snv_alt-T_$BED\_merged.cdt
python $SUMMAT -1 $OUTDIR/$BED/CDT/dbSnp153_snv_alt-G_$BED\_same.cdt -2 $OUTDIR/$BED/CDT/dbSnp153_snv_alt-C_$BED\_opposite.cdt -o $OUTDIR/$BED/CDT/dbSnp153_snv_alt-G_$BED\_merged.cdt
python $SUMMAT -1 $OUTDIR/$BED/CDT/dbSnp153_snv_alt-C_$BED\_same.cdt -2 $OUTDIR/$BED/CDT/dbSnp153_snv_alt-G_$BED\_opposite.cdt -o $OUTDIR/$BED/CDT/dbSnp153_snv_alt-C_$BED\_merged.cdt

# Merge ref with strand correction
python $SUMMAT -1 $OUTDIR/$BED/CDT/dbSnp153_snv_ref-A_$BED\_same.cdt -2 $OUTDIR/$BED/CDT/dbSnp153_snv_ref-T_$BED\_opposite.cdt -o $OUTDIR/$BED/CDT/dbSnp153_snv_ref-A_$BED\_merged.cdt
python $SUMMAT -1 $OUTDIR/$BED/CDT/dbSnp153_snv_ref-T_$BED\_same.cdt -2 $OUTDIR/$BED/CDT/dbSnp153_snv_ref-A_$BED\_opposite.cdt -o $OUTDIR/$BED/CDT/dbSnp153_snv_ref-T_$BED\_merged.cdt
python $SUMMAT -1 $OUTDIR/$BED/CDT/dbSnp153_snv_ref-G_$BED\_same.cdt -2 $OUTDIR/$BED/CDT/dbSnp153_snv_ref-C_$BED\_opposite.cdt -o $OUTDIR/$BED/CDT/dbSnp153_snv_ref-G_$BED\_merged.cdt
python $SUMMAT -1 $OUTDIR/$BED/CDT/dbSnp153_snv_ref-C_$BED\_same.cdt -2 $OUTDIR/$BED/CDT/dbSnp153_snv_ref-G_$BED\_opposite.cdt -o $OUTDIR/$BED/CDT/dbSnp153_snv_ref-C_$BED\_merged.cdt

for CDT in $OUTDIR/$BED/CDT/dbSnp153_snv_*_merged.cdt:
do
	BASE=`basename $CDT ".cdt"`
	java -jar $SCRIPTMANAGER figure-generation heatmap -p 0.95 --black $CDT -o $OUTDIR/$BED/PNG/$BASE.png
	java -jar $SCRIPTMANAGER figure-generation label-heatmap $OUTDIR/$BED/PNG/$BASE.png \
		-f 20 -l -500 -m 0 -r 500 -o $OUTDIR/$BED/SVG/$BASE.svg
done

perl $COMPOSITE $OUTDIR/$BED/CDT/dbSnp153_snv_alt-A_$BED\_merged.cdt $OUTDIR/$BED/Composites/dbSnp153_snv_alt-A_$BED\_merged.out
perl $COMPOSITE $OUTDIR/$BED/CDT/dbSnp153_snv_alt-T_$BED\_merged.cdt $OUTDIR/$BED/Composites/dbSnp153_snv_alt-T_$BED\_merged.out
perl $COMPOSITE $OUTDIR/$BED/CDT/dbSnp153_snv_alt-C_$BED\_merged.cdt $OUTDIR/$BED/Composites/dbSnp153_snv_alt-C_$BED\_merged.out
perl $COMPOSITE $OUTDIR/$BED/CDT/dbSnp153_snv_alt-G_$BED\_merged.cdt $OUTDIR/$BED/Composites/dbSnp153_snv_alt-G_$BED\_merged.out

perl $COMPOSITE $OUTDIR/$BED/CDT/dbSnp153_snv_ref-A_$BED\_merged.cdt $OUTDIR/$BED/Composites/dbSnp153_snv_ref-A_$BED\_merged.out
perl $COMPOSITE $OUTDIR/$BED/CDT/dbSnp153_snv_ref-T_$BED\_merged.cdt $OUTDIR/$BED/Composites/dbSnp153_snv_ref-T_$BED\_merged.out
perl $COMPOSITE $OUTDIR/$BED/CDT/dbSnp153_snv_ref-G_$BED\_merged.cdt $OUTDIR/$BED/Composites/dbSnp153_snv_ref-G_$BED\_merged.out
perl $COMPOSITE $OUTDIR/$BED/CDT/dbSnp153_snv_ref-C_$BED\_merged.cdt $OUTDIR/$BED/Composites/dbSnp153_snv_ref-C_$BED\_merged.out

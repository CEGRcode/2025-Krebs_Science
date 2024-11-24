#!/bin/bash

# Make PE insert size histogram of BNase-seq, MNase-seq ([21 U],[304 U]), and DNase-seq
# (1c) BNase-seq
# (1e) MNase-seq
# (1f) DNase-seq
# see 04_Figures/Fig_1c.sh
# see 04_Figures/Fig_1e_f.sh
# (F1b and S1) Use peak-align to "pileup" CpG island annotations on TSS-centered RefPT
# (F1b) Tag Pileup deep ENCODE MNase-seq signal on TSS-centered RefPT

### CHANGE ME
WRK=/path/to/2024-Krebs_Science/Z_Figures
WRK=/ocean/projects/see180003p/owlang/2024-Krebs_Science/Z_Figures
WRK=/storage/home/owl5022/scratch/2024-Krebs_Science/Z_Figures
THREADS=4
###

# Dependencies
# - java
# - pandas
# - python
# - seaborn

set -exo
module load anaconda
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Inputs and outputs
LIBRARY=$WRK/../X_Bulk_Processing/Library
BAMDIR=$WRK/../data/BAM

# Script shortcuts
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar
COMPOSITE=$WRK/../bin/sum_Col_CDT.pl
BOXPLOT=$WRK/../bin/make_box_plot.py


[ -d F7 ] || mkdir F7

# ===============================================================================================================================

[ -d F7/a ] || mkdir F7/a

# Heatmaps
BED=TSS_GROUP-Expressed_SORT-CpG_2000bp
cp $LIBRARY/$BED/SVG/BNase-seq_50U-10min_merge_hg38_${BED}_midpoint_TotalTag_combined.svg F7/a/
cp $LIBRARY/$BED/SVG/DNase-seq_ENCFF518XTC_rep1_hg38_${BED}_midpoint_TotalTag_combined.svg F7/a/

# =====Make CpG heatmaps=====

REFPT=$WRK/../data/RefPT-Krebs/2000bp/TSS_GROUP-Expressed_SORT-CpG_2000bp.bed
BED=`basename $REFPT ".bed"`
CPGBED=$WRK/../data/RefPT-Other/CpGIslands.bed
BASE=CpG-Islands_$BED

# Pileup CpG islands
java -jar $SCRIPTMANAGER peak-analysis peak-align-ref $CPGBED $REFPT -o F7/a/$BASE

# Two-color heatmap
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation heatmap -a 1 --blue F7/a/$BASE\_combined.cdt -o F7/a/$BASE\_treeview.png

# Count sites
NSITES=`wc -l $REFPT | awk '{print $1-1}'`

# Label heatmap
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation label-heatmap F7/a/$BASE\_treeview.png \
	-l "-1" -m "0" -r "+1" -w 2 -f 18 \
	-x "Distance from TSS (kb)" -y "${NSITES} CoPRO determined TSSs sorted by CpG island length" \
	-o F7/a/$BASE\_treeview.svg

# Copy same file into supplement
[ -d S1 ] || mkdir S1
cp F7/a/$BASE\_treeview.svg S1/$BASE\_treeview.svg


# =====Pileup MNase heatmap (TSS RefPT w/ 80bp shift)=====

BAMFILE=$WRK/../data/BAM/MNase-seq_ENCODE_merge_hg38.bam
BAM=`basename $BAMFILE ".bam"`
BASE=${BAM}_${BED}_TotalTag

# Pileup SE MNase data (shift 80bp)
java -jar $SCRIPTMANAGER read-analysis tag-pileup $REFPT $BAMFILE \
	-5 -1 --shift 80 --combined --cpu $THREADS -z \
	-o F7/a/$BASE\_combined.out -M F7/a/$BASE

# Two-color heatmap
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation heatmap --black -p 0.95 F7/a/$BASE\_combined.cdt.gz -o F7/a/$BASE\_combined.png

# Label heatmap
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation label-heatmap F7/a/$BASE\_combined.png \
	-l "-1" -m "0" -r "+1" -w 2 -f 18 \
	-x "Distance from TSS (kb)" -y "${NSITES} CoPRO determined TSSs sorted by CpG island length" \
	-o F7/a/$BASE\_combined.svg


# =====Pileup MNase/BNase for violins=====

# Set BAM basenames
MNASE_BAM=MNase-seq_ENCODE_merge_hg38
BNASE_BAM=BNase-seq_50U-10min_merge_hg38

# Parse scaling factors
MNASE_FACTOR=`grep 'factor' $BAMDIR/NormalizationFactors/${MNASE_BAM}_TotalTag_ScalingFactors.out | awk '{print $3}'`
BNASE_FACTOR=`grep 'factor' $BAMDIR/NormalizationFactors/${BNASE_BAM}_TotalTag_ScalingFactors.out | awk '{print $3}'`

DATAFILE=F7/a/violin_data.txt
[ -f $DATAFILE ] && rm $DATAFILE

for GROUP in "NoOverlap" "Overlap";
do
	BEDFILE=../data/RefPT-Krebs/PlusOneDyad_SORT-DistToExpressedTSS_GROUP-CpGIsland-$GROUP.bed

	# Tag pileup BNase-seq and MNase-seq
	java -jar $SCRIPTMANAGER read-analysis tag-pileup -1 -5 --combined --shift 80 $BEDFILE $BAMDIR/$MNASE_BAM.bam -M F7/a/MNase_$GROUP\_5read1-SHIFT80
	java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -p --combined            $BEDFILE $BAMDIR/$BNASE_BAM.bam -M F7/a/BNase_$GROUP\_midpoint
	# java -jar $SCRIPTMANAGER read-analysis tag-pileup -1 -5 --combined --shift 80 $BEDFILE $BNASE_BAMFILE -M F7/a/BNase_$GROUP\_5read1-SHIFT80

	# Aggregate data into violin data file
	java -jar $SCRIPTMANAGER read-analysis aggregate-data -m F7/a/MNase_$GROUP\_5read1-SHIFT80_combined.cdt -o F7/a/MNase_$GROUP.out
	java -jar $SCRIPTMANAGER read-analysis aggregate-data -m F7/a/BNase_$GROUP\_midpoint_combined.cdt       -o F7/a/BNase_$GROUP.out

	# Scale data
	java -jar $SCRIPTMANAGER read-analysis scale-matrix -l 1 -s $MNASE_FACTOR F7/a/MNase_$GROUP.out -o F7/a/MNase_${GROUP}_SCALE.out
	java -jar $SCRIPTMANAGER read-analysis scale-matrix -l 1 -s $BNASE_FACTOR F7/a/BNase_$GROUP.out -o F7/a/BNase_${GROUP}_SCALE.out

	# Append aggregated data into a merged violin data file
	sed '1d' F7/a/MNase_${GROUP}_SCALE.out | awk -v GROUP=$GROUP 'BEGIN{OFS="\t";FS="\t"}{print $2,"MNase-"GROUP}' >> $DATAFILE
	sed '1d' F7/a/BNase_${GROUP}_SCALE.out | awk -v GROUP=$GROUP 'BEGIN{OFS="\t";FS="\t"}{print $2,"BNase-"GROUP}' >> $DATAFILE
done

# Plot violin data (skip zeros)
python $BOXPLOT -i <(awk 'BEGIN{OFS="\t";FS="\t"}{if ($1!=0) print }' $DATAFILE) \
	-o F7/a/violin_data.svg \
	--xlabel "+1 dyad group" --ylabel "Tag Occupancy" \
	--width 4 --height 4 --ylog10 --preset1 \
	--title "BNase v MNase in CpGIsland overlap v non-overlap regions"

# Stat zeros (all MNase)
awk '{if ($1==0) print $2}' F7/a/violin_data.txt | sort | uniq -c
# 28 MNase-NoOverlap
# 1026 MNase-Overlap

# ===============================================================================================================================

[ -d F7/c ] || mkdir F7/c


# Composites
BED=PlusOneDyad_SORT-Expression_2000bp
cp $LIBRARY/$BED/Composites/BNase-ChIP_H2A_merge_hg38_$BED\_midpoint-MAX80_TotalTag_combined.out F7/c/BNase-ChIP-H2A.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H2AZ_merge_hg38_$BED\_midpoint-MAX80_TotalTag_combined.out F7/c/BNase-ChIP-H2AZ.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3_merge_hg38_$BED\_midpoint-MAX80_TotalTag_combined.out F7/c/BNase-ChIP-H3.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3K4me3_merge_hg38_$BED\_midpoint-MAX80_TotalTag_combined.out F7/c/BNase-ChIP-H3K4me3.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3K9ac_merge_hg38_$BED\_midpoint-MAX80_TotalTag_combined.out F7/c/BNase-ChIP-H3K9ac.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3K27ac_merge_hg38_$BED\_midpoint-MAX80_TotalTag_combined.out F7/c/BNase-ChIP-H3K27ac.out
cp $LIBRARY/$BED/Composites/MNase-ChIP_H3K4me3_merge_hg38_$BED\_midpoint-MAX80_TotalTag_combined.out F7/c/MNase-ChIP-H3K4me3.out


# ===============================================================================================================================

[ -d F7/d ] || mkdir F7/d

# Heatmaps
BED=PlusOneDyad_SORT-pHN-dHN_400bp
cp $LIBRARY/$BED/SVG/BNase-ChIP_H3K4me3_merge_hg38_$BED\_midpoint-MAX80_TotalTag_combined.svg F7/d
cp $LIBRARY/$BED/SVG/BNase-ChIP_H3K9ac_merge_hg38_$BED\_midpoint-MAX80_TotalTag_combined.svg F7/d
cp $LIBRARY/$BED/SVG/BNase-ChIP_H3K27ac_merge_hg38_$BED\_midpoint-MAX80_TotalTag_combined.svg F7/d
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation label-heatmap $LIBRARY/$BED/PNG/Strand/CoPRO_Capped_merge_hg38_$BED\_5read1_anti_treeview.png \
	-l -200 -m 0 -r +200 -w 2 -f 18 \
	-x $BED -y "$BED occurences (NSITES sites)" \
	-o F7/d/CoPRO_Capped_merge_hg38_$BED\_5read1_anti_treeview_label.svg

# Composites
BED=PlusOneDyad_SORT-pHN-dHN_GROUP-TOP-2500_400bp
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3K4me3_merge_hg38_$BED\_midpoint-MAX80_TotalTag_combined.out F7/d/Top-H3K4me3.out
cp $LIBRARY/$BED/Composites/CoPRO_Capped_merge_hg38_$BED\_5read1.out F7/d/Top-CoPRO.out

BED=PlusOneDyad_SORT-pHN-dHN_GROUP-BOTTOM-2500_400bp
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3K4me3_merge_hg38_$BED\_midpoint-MAX80_TotalTag_combined.out F7/d/Bottom-H3K4me3.out
cp $LIBRARY/$BED/Composites/CoPRO_Capped_merge_hg38_$BED\_5read1.out F7/d/Bottom-CoPRO.out
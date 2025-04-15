#!/bin/bash

# data/RefPT-JASPAR-nonK562
#   |--<TFNAME>_<JASPARID>_Unbound.bed
#   |--<TFNAME>_<JASPARID>_K562-specific-Unbound.bed

### CHANGE ME
WRK=/storage/group/bfp2/default/hxc585_HainingChen/2025_Chen_TF-Nuc/03_Call_JASPAR
#WRK=/ocean/projects/see180003p/owlang/2024-Krebs_Science/03_Call_JASPAR
#WRK=/scratch/owl5022/2024-Krebs_Science/03_Call_JASPAR

# Dependencies
# - bedtools
# - java

#set -exo
module load bedtools
module load anaconda3
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Inputs and outputs
GENOME=$WRK/../data/hg38_files/hg38.fa
BAMFILE=$WRK/../data/BAM/BNase-seq_50U-10min_merge_hg38.bam
BLACKLIST=$WRK/../data/hg38_files/ENCFF356LFX_hg38_exclude.bed
MOTIF=$WRK/NonK562_narrowPeaks
ODIR=$WRK/../data/RefPT-JASPAR-nonK562
# Script shortcuts
DEDUP=$WRK/../bin/dedup_coord_by_ID.py
RATIO=$WRK/../bin/calculate_BED_ScoreRatio.pl
UPDATES=$WRK/../bin/update_BED_score_with_TAB_score.pl
ORIGINAL_SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar

# Parse BED 6 + 3 narroPeak files into summits (1bp)

for file in $WRK/NonK562_narrowPeaks/*.bed.gz ; do
	filename=`basename $file ".bed.gz"`
	gzip -dc $file > $WRK/NonK562_narrowPeaks/${filename}.bed
	awk '{OFS="\t"}{FS="\t"}{print $1,$2+$10,$2+$10+1,"P-"$1"_"$2"_"$3,$5,$6}' $WRK/NonK562_narrowPeaks/${filename}.bed >  $WRK/NonK562_narrowPeaks/${filename}_summit.bed
	java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $WRK/NonK562_narrowPeaks/${filename}_summit.bed -o $WRK/NonK562_narrowPeaks/${filename}_1000bp.bed
	rm $WRK/NonK562_narrowPeaks/${filename}_summit.bed
done


## take unbound TF sites
cat $WRK/NonK562_narrowPeaks/CTCF_*_1000bp.bed | bedtools intersect -v -a $WRK/FIMO/CTCF_MA1929.1/filtered.bed -b - | bedtools intersect -v -a - -b $WRK/Intersect/CTCF_MA1929.1/BoundMotifs_SORT-TFnucRatio.bed > $ODIR/CTCF_MA1929.1_Unbound.bed
cat $WRK/NonK562_narrowPeaks/ZKSCAN1_*_1000bp.bed | bedtools intersect -v -a $WRK/FIMO/ZKSCAN1_MA1585.1/filtered.bed -b - | bedtools intersect -v -a - -b $WRK/Intersect/ZKSCAN1_MA1585.1/BoundMotifs_SORT-TFnucRatio.bed > $ODIR/ZKSCAN1_MA1585.1_Unbound.bed
cat $WRK/NonK562_narrowPeaks/ATF2_*_1000bp.bed | bedtools intersect -v -a $WRK/FIMO/ATF2_MA1632.2/filtered.bed -b - | bedtools intersect -v -a - -b $WRK/Intersect/ATF2_MA1632.2/BoundMotifs_SORT-TFnucRatio.bed > $ODIR/ATF2_MA1632.2_Unbound.bed
cat $WRK/NonK562_narrowPeaks/CREM_*_1000bp.bed | bedtools intersect -v -a $WRK/FIMO/CREM_MA0609.3/filtered.bed -b - | bedtools intersect -v -a - -b $WRK/Intersect/CREM_MA0609.3/BoundMotifs_SORT-TFnucRatio.bed > $ODIR/CREM_MA0609.3_Unbound.bed
cat $WRK/NonK562_narrowPeaks/CREM_*_1000bp.bed | bedtools intersect -v -a $WRK/FIMO/CREM_MA0609.3/filtered.bed -b - | bedtools intersect -v -a - -b $WRK/Intersect/CREM_MA0609.3/BoundMotifs_SORT-TFnucRatio.bed > $ODIR/CREM_MA0609.3_Unbound.bed
cat $WRK/NonK562_narrowPeaks/ZNF263_*_1000bp.bed  | bedtools intersect -v -a $WRK/FIMO/ZNF263_MA0528.1/filtered.bed -b - | bedtools intersect -v -a - -b $WRK/Intersect/ZNF263_MA0528.1/BoundMotifs_SORT-TFnucRatio.bed > $ODIR/ZNF263_MA0528.1_Unbound.bed

## take unbound K562 specific TF sites

cat $WRK/NonK562_narrowPeaks/CTCF_*_1000bp.bed  $WRK/Intersect/CTCF_MA1929.1/BoundMotifs_SORT-TFnucRatio.bed | bedtools intersect -u -a $WRK/FIMO/CTCF_MA1929.1/filtered.bed -b -  | bedtools intersect -v -a - -b $WRK/Intersect/CTCF_MA1929.1/BoundMotifs_SORT-TFnucRatio.bed > $ODIR/CTCF_MA1929.1_K562-specific-Unbound.bed
cat $WRK/NonK562_narrowPeaks/ZKSCAN1_*_1000bp.bed $WRK/Intersect/ZKSCAN1_MA1585.1/BoundMotifs_SORT-TFnucRatio.bed | bedtools intersect -u -a $WRK/FIMO/ZKSCAN1_MA1585.1/filtered.bed -b - | bedtools intersect -v -a - -b $WRK/Intersect/ZKSCAN1_MA1585.1/BoundMotifs_SORT-TFnucRatio.bed > $ODIR/ZKSCAN1_MA1585.1_K562-specific-Unbound.bed
cat $WRK/NonK562_narrowPeaks/ATF2_*_1000bp.bed $WRK/Intersect/ATF2_MA1632.2/BoundMotifs_SORT-TFnucRatio.bed | bedtools intersect -u -a $WRK/FIMO/ATF2_MA1632.2/filtered.bed -b - | bedtools intersect -v -a - -b $WRK/Intersect/ATF2_MA1632.2/BoundMotifs_SORT-TFnucRatio.bed > $ODIR/ATF2_MA1632.2_K562-specific-Unbound.bed
cat $WRK/NonK562_narrowPeaks/CREM_*_1000bp.bed $WRK/Intersect/CREM_MA0609.3/BoundMotifs_SORT-TFnucRatio.bed | bedtools intersect -u -a $WRK/FIMO/CREM_MA0609.3/filtered.bed -b - | bedtools intersect -v -a - -b $WRK/Intersect/CREM_MA0609.3/BoundMotifs_SORT-TFnucRatio.bed > $ODIR/CREM_MA0609.3_K562-specific-Unbound.bed
cat $WRK/NonK562_narrowPeaks/CREM_*_1000bp.bed $WRK/Intersect/CREM_MA0609.3/BoundMotifs_SORT-TFnucRatio.bed | bedtools intersect -u -a $WRK/FIMO/CREM_MA0609.3/filtered.bed -b - | bedtools intersect -v -a - -b $WRK/Intersect/CREM_MA0609.3/BoundMotifs_SORT-TFnucRatio.bed > $ODIR/CREM_MA0609.3_K562-specific-Unbound.bed
cat $WRK/NonK562_narrowPeaks/ZNF263_*_1000bp.bed $WRK/Intersect/ZNF263_MA0528.1/BoundMotifs_SORT-TFnucRatio.bed | bedtools intersect -u -a $WRK/FIMO/ZNF263_MA0528.1/filtered.bed -b - | bedtools intersect -v -a - -b $WRK/Intersect/ZNF263_MA0528.1/BoundMotifs_SORT-TFnucRatio.bed > $ODIR/ZNF263_MA0528.1_K562-specific-Unbound.bed

mkdir -p $WRK/../data/RefPT-JASPAR-nonK562/1000bp
for file in $WRK/../data/RefPT-JASPAR-nonK562/*.bed ; do
	filename=`basename $file ".bed"
	java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $file -o "$WRK/../data/RefPT-JASPAR-nonK562/1000bp/${filename}_1000bp.bed"
done


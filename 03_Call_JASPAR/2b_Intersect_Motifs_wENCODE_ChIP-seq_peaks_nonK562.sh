#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 01:00:00
#SBATCH -A open
#SBATCH -o logs/2b_Intersect_Motifs_wENCODE_ChIP-seq_peaks_nonK562.log.out-%a
#SBATCH -e logs/2b_Intersect_Motifs_wENCODE_ChIP-seq_peaks_nonK562.log.err-%a

# data/RefPT-JASPAR-nonK562
#   |--<TFNAME>_<JASPARID>_Unbound.bed
#   |--<TFNAME>_<JASPARID>_K562-specific-Unbound.bed

# Dependencies
# - bedtools
# - java

set -exo
source activate bx

# Inputs and outputs
GENOME=../data/hg38_files/hg38.fa
BAMFILE=../data/BAM/BNase-seq_50U-10min_merge_hg38.bam
BLACKLIST=../data/hg38_files/ENCFF356LFX_hg38_exclude.bed
MOTIF=NonK562_narrowPeaks
ODIR=../data/RefPT-JASPAR-nonK562
# Script shortcuts
DEDUP=../bin/dedup_coord_by_ID.py
RATIO=../bin/calculate_BED_ScoreRatio.pl
UPDATES=../bin/update_BED_score_with_TAB_score.pl
ORIGINAL_SCRIPTMANAGER=../bin/ScriptManager-v0.15.jar
SCRIPTMANAGER=../bin/ScriptManager-v0.15.jar

[[ -d $ODIR/1000bp ]] || mkdir -p $ODIR/1000bp

# Parse BED 6 + 3 narrowPeak files into summits (1bp)

for file in NonK562_narrowPeaks/*.bed.gz ; do
	filename=`basename $file ".bed.gz"`
	gzip -dc $file > NonK562_narrowPeaks/${filename}.bed
	awk '{OFS="\t"}{FS="\t"}{print $1,$2+$10,$2+$10+1,"P-"$1"_"$2"_"$3,$5,$6}' NonK562_narrowPeaks/${filename}.bed >  NonK562_narrowPeaks/${filename}_summit.bed
	java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 NonK562_narrowPeaks/${filename}_summit.bed -o NonK562_narrowPeaks/${filename}_1000bp.bed
	rm NonK562_narrowPeaks/${filename}_summit.bed
done


## take unbound TF sites
cat NonK562_narrowPeaks/CTCF_*_1000bp.bed | bedtools intersect -v -a FIMO/CTCF_MA1929.1/filtered.bed -b - | bedtools intersect -v -a - -b Intersect/CTCF_MA1929.1/BoundMotifs_SORT-TFnucRatio.bed > $ODIR/CTCF_MA1929.1_Unbound.bed
cat NonK562_narrowPeaks/ZKSCAN1_*_1000bp.bed | bedtools intersect -v -a FIMO/ZKSCAN1_MA1585.1/filtered.bed -b - | bedtools intersect -v -a - -b Intersect/ZKSCAN1_MA1585.1/BoundMotifs_SORT-TFnucRatio.bed > $ODIR/ZKSCAN1_MA1585.1_Unbound.bed
cat NonK562_narrowPeaks/ATF2_*_1000bp.bed | bedtools intersect -v -a FIMO/ATF2_MA1632.2/filtered.bed -b - | bedtools intersect -v -a - -b Intersect/ATF2_MA1632.2/BoundMotifs_SORT-TFnucRatio.bed > $ODIR/ATF2_MA1632.2_Unbound.bed
cat NonK562_narrowPeaks/CREM_*_1000bp.bed | bedtools intersect -v -a FIMO/CREM_MA0609.3/filtered.bed -b - | bedtools intersect -v -a - -b Intersect/CREM_MA0609.3/BoundMotifs_SORT-TFnucRatio.bed > $ODIR/CREM_MA0609.3_Unbound.bed
cat NonK562_narrowPeaks/CREM_*_1000bp.bed | bedtools intersect -v -a FIMO/CREM_MA0609.3/filtered.bed -b - | bedtools intersect -v -a - -b Intersect/CREM_MA0609.3/BoundMotifs_SORT-TFnucRatio.bed > $ODIR/CREM_MA0609.3_Unbound.bed
cat NonK562_narrowPeaks/ZNF263_*_1000bp.bed  | bedtools intersect -v -a FIMO/ZNF263_MA0528.1/filtered.bed -b - | bedtools intersect -v -a - -b Intersect/ZNF263_MA0528.1/BoundMotifs_SORT-TFnucRatio.bed > $ODIR/ZNF263_MA0528.1_Unbound.bed

## take unbound K562 specific TF sites

cat NonK562_narrowPeaks/CTCF_*_1000bp.bed  Intersect/CTCF_MA1929.1/BoundMotifs_SORT-TFnucRatio.bed | bedtools intersect -u -a FIMO/CTCF_MA1929.1/filtered.bed -b -  | bedtools intersect -v -a - -b Intersect/CTCF_MA1929.1/BoundMotifs_SORT-TFnucRatio.bed > $ODIR/CTCF_MA1929.1_K562-specific-Unbound.bed
cat NonK562_narrowPeaks/ZKSCAN1_*_1000bp.bed Intersect/ZKSCAN1_MA1585.1/BoundMotifs_SORT-TFnucRatio.bed | bedtools intersect -u -a FIMO/ZKSCAN1_MA1585.1/filtered.bed -b - | bedtools intersect -v -a - -b Intersect/ZKSCAN1_MA1585.1/BoundMotifs_SORT-TFnucRatio.bed > $ODIR/ZKSCAN1_MA1585.1_K562-specific-Unbound.bed
cat NonK562_narrowPeaks/ATF2_*_1000bp.bed Intersect/ATF2_MA1632.2/BoundMotifs_SORT-TFnucRatio.bed | bedtools intersect -u -a FIMO/ATF2_MA1632.2/filtered.bed -b - | bedtools intersect -v -a - -b Intersect/ATF2_MA1632.2/BoundMotifs_SORT-TFnucRatio.bed > $ODIR/ATF2_MA1632.2_K562-specific-Unbound.bed
cat NonK562_narrowPeaks/CREM_*_1000bp.bed Intersect/CREM_MA0609.3/BoundMotifs_SORT-TFnucRatio.bed | bedtools intersect -u -a FIMO/CREM_MA0609.3/filtered.bed -b - | bedtools intersect -v -a - -b Intersect/CREM_MA0609.3/BoundMotifs_SORT-TFnucRatio.bed > $ODIR/CREM_MA0609.3_K562-specific-Unbound.bed
cat NonK562_narrowPeaks/CREM_*_1000bp.bed Intersect/CREM_MA0609.3/BoundMotifs_SORT-TFnucRatio.bed | bedtools intersect -u -a FIMO/CREM_MA0609.3/filtered.bed -b - | bedtools intersect -v -a - -b Intersect/CREM_MA0609.3/BoundMotifs_SORT-TFnucRatio.bed > $ODIR/CREM_MA0609.3_K562-specific-Unbound.bed
cat NonK562_narrowPeaks/ZNF263_*_1000bp.bed Intersect/ZNF263_MA0528.1/BoundMotifs_SORT-TFnucRatio.bed | bedtools intersect -u -a FIMO/ZNF263_MA0528.1/filtered.bed -b - | bedtools intersect -v -a - -b Intersect/ZNF263_MA0528.1/BoundMotifs_SORT-TFnucRatio.bed > $ODIR/ZNF263_MA0528.1_K562-specific-Unbound.bed

for file in $ODIR/*.bed ; do
	filename=`basename $file ".bed"
	java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $file -o $ODIR/1000bp/${filename}_1000bp.bed
done


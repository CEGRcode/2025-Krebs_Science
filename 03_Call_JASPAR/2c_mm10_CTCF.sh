#!/bin/bash

# data/RefPT-JASPAR-nonK562
#   |--<TFNAME>_<JASPARID>_Unbound.bed
#   |--<TFNAME>_<JASPARID>_K562-specific-Unbound.bed

### CHANGE ME
WRK=/storage/group/bfp2/default/hxc585_HainingChen/2025_Chen_TF-Nuc/03_Call_JASPAR
#WRK=/ocean/projects/see180003p/owlang/2024-Krebs_Science/03_Call_JASPAR
#WRK=/scratch/owl5022/2024-Krebs_Science/03_Call_JASPAR

# - java

#set -exo
module load bedtools
module load anaconda3
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Inputs and outputs
GENOME=$WRK/../data/mm10_files/mm10.fa
BLACKLIST=$WRK/../data/mm10_files/mm10_blacklist.bed
# Script shortcuts
DEDUP=$WRK/../bin/dedup_coord_by_ID.py
RATIO=$WRK/../bin/calculate_BED_ScoreRatio.pl
UPDATES=$WRK/../bin/update_BED_score_with_TAB_score.pl
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar
OUTDIR=$WRK/FIMO/CTCF_MA1929.1_mm10
mkdir -p $OUTDIR
## FIMO ctcf motif agaist mm10 genome
fimo --verbosity 1 --thresh 1.0E-4 --max-strand --oc $OUTDIR $WRK/../data/JASPAR/CTCF_MA1929.1.meme $GENOME
java -jar $SCRIPTMANAGER coordinate-manipulation gff-to-bed $OUTDIR/fimo.gff -o $OUTDIR/fimo_unformatted_unfiltered.bed
awk '{OFS="\t"}{FS="\t"}{print $1,$2,$3,"M-"$1"_"$2"_"$3,$5,$6}' $OUTDIR/fimo_unformatted_unfiltered.bed \
	| sort -k4 -rnk5 > $OUTDIR/fimo_unfiltered.bed

python $DEDUP -i $OUTDIR/fimo_unfiltered.bed -o $OUTDIR/fimo.bed

java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $OUTDIR/fimo.bed -o $OUTDIR/fimo_1000bp.bed

bedtools intersect -v -a <(awk '{if($2>0) print}' $OUTDIR/fimo_1000bp.bed) -b $BLACKLIST > $OUTDIR/filtered.bed

gunzip -c $WRK/NonK562_narrowPeaks/CTCF_J1_ENCFF533APC.bed.gz > $OUTDIR/ENCFF533APC.bed
#shuffle bedfiles
shuf $OUTDIR/ENCFF533APC.bed > $OUTDIR/ENCFF533APC_shuffled.bed
shuf $OUTDIR/filtered.bed > $OUTDIR/MA1929_1_final_1000bp_shuffled.bed
#expand befiles
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000  $OUTDIR/ENCFF533APC_shuffled.bed -o $OUTDIR/ENCFF533APC_shuffled_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 20 $OUTDIR/MA1929_1_final_1000bp_shuffled.bed -o $OUTDIR/MA1929_1_final_1000bp_shuffled_20bp.bed
# Intersect peaks with motifs - filter to keep overlap - move ENCODE ChIP value (signal value) to score col - sort by ID, then score
bedtools intersect -loj -a $OUTDIR/MA1929_1_final_1000bp_shuffled_20bp.bed -b $OUTDIR/ENCFF533APC_shuffled_1000bp.bed | awk '{OFS="\t"}{FS="\t"}{if($8>0) print $1,$2,$3,$4,$13,$6}' | sort -rnk4,5 > $OUTDIR/MA1929_1_final_1000bp_intersected_wDUP.bed
bedtools intersect -loj -a $OUTDIR/MA1929_1_final_1000bp_shuffled_20bp.bed -b $OUTDIR/ENCFF533APC_shuffled_1000bp.bed | awk '{OFS="\t"}{FS="\t"}{if($8==-1) print $1,$2,$3,$4,$13,$6}' > $OUTDIR/MA1929_1_final_1000bp_NOTintersected_wDUP.bed
#Deduplicate bound motifs by keeping first instance (larger ENCODE score based on previous command sort)
python $DEDUP -i $OUTDIR/MA1929_1_final_1000bp_intersected_wDUP.bed -o $OUTDIR/MA1929_1_final_1000bp_intersected.bed
#Deduplicate of unbound motifs does  NOT work as each sites seems to have its own unique 4th column
#get number of rows from intersected bedfile
#expand intersected bedfile
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 164 MA1929_1_final_1000bp_intersected.bed -o $OUTDIR/MA1929_1_final_1000bp_intersected_164bp.bed
#do size selection so that nucleosomal fragments (126-185 bp) are pooled from DNAse_FLASH BAM file
samtools view -h $WRK/data/BAM/MPE-seq_10min_merge_rep1_mm10.bam | awk 'BEGIN {OFS=t} /^@/ || ($9 >= 126 && $9 <= 185) || ($9 <= -126 && $9 >= -185)' | samtools view -b -o $OUTDIR/MPE-seq_10min_merge_rep1_mm10_126to185bp.bam
java -jar $WRK/../bin/picard.jar BuildBamIndex --INPUT $OUTDIR/MPE-seq_10min_merge_rep1_mm10_126to185bp.bam --OUTPUT $OUTDIR/MPE-seq_10min_merge_rep1_mm10_126to185bp.bam.bai

#do initial tag-pileUp (output is input directory). Settings: midpoint(m) OR 5 prime end (-5) with read 1 (-1), Gizp output cdt (z), No smoothing (N), required proper PEs (p), load blacklist **total tag option (-t) removed**
java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -M $OUTDIR/MPE-seq_10min_merge_rep1_mm10_126to185bp_MA1929_1_final_1000bp_intersected_164bp -N -p --cpu 4 --blacklist-filter $BLACKLIST -o $OUTDIR/MPE-seq_10min_merge_rep1_mm10_126to185bp_MA1929_1_final_1000bp_intersected_164bp_midpoint.out $OUTDIR/MA1929_1_final_1000bp_intersected_164bp.bed $OUTDIR/MPE-seq_10min_merge_rep1_mm10_126to185bp.bam
#sum the number of tags by each row
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -o $OUTDIR/MPE-seq_10min_merge_rep1_mm10_126to185bp_MA1929_1_final_1000bp_intersected_164bp_combined_sum.tsv  $OUTDIR/MPE-seq_10min_merge_rep1_mm10_126to185bp_MA1929_1_final_1000bp_intersected_164bp_combined.cdt
#remove header from CDT1_sum file
cat $OUTDIR/MPE-seq_10min_merge_rep1_mm10_126to185bp_MA1929_1_final_1000bp_intersected_164bp_combined_sum.tsv | sed '1d' > $OUTDIR/MPE-seq_10min_merge_rep1_mm10_126to185bp_MA1929_1_final_1000bp_intersected_164bp_combined_sum_noHeader.tsv
#paste bedfile and CTD1_sum_noHeader, make sure all rows match first, avoid any rows with 0 in TF signal (column 5) or nucleosome occupancy (column 8) then divide encode TF signal to nucleosome occupancy (ratio) in column 7 and sort
paste $OUTDIR/MA1929_1_final_1000bp_intersected_164bp.bed $OUTDIR/MPE-seq_10min_merge_rep1_mm10_126to185bp_MA1929_1_final_1000bp_intersected_164bp_combined_sum_noHeader.tsv | awk '{OFS="\t"}{FS="\t"}{if ($12=$7 && $5!=0 && $8!=0) print $1,$2,$3,$4,$5,$6,($5/$8)}' | sort -k7,7n > $OUTDIR/MA1929_1_final_1000bp_intersected_ratio.TSV
#get number of rows from intersected bedfile
cat $OUTDIR/MA1929_1_final_1000bp_intersected_ratio.TSV | wc -l | awk '{printf "%.f\n", $1 * 0.25}' > $OUTDIR/MA1929_1_final_1000bp_intersected_rowsNumber.tab
cat $OUTDIR/MA1929_1_final_1000bp_intersected_ratio.TSV | wc -l | awk '{printf "%.f\n", $1 * 0.5}' > $OUTDIR/MA1929_1_final_1000bp_intersected_rowsNumber2.tab
cat $OUTDIR/MA1929_1_final_1000bp_intersected_ratio.TSV | wc -l | awk '{printf "%.f\n", $1 * 0.75}' > $OUTDIR/MA1929_1_final_1000bp_intersected_rowsNumber3.tab
#take sorted sites and split into quartiles
#take above motif dedup bedfile that is sorted by TSV.
cat $OUTDIR/MA1929_1_final_1000bp_intersected_ratio.TSV | head -$(cat $OUTDIR/MA1929_1_final_1000bp_intersected_rowsNumber.tab)  | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > $WRK/../data/RefPT-JASPAR-nonK562/MA1929_1_mm10_intersected_MPE-seq10min_164bp_category1.bed
cat $OUTDIR/MA1929_1_final_1000bp_intersected_ratio.TSV | head -$(cat $OUTDIR/MA1929_1_final_1000bp_intersected_rowsNumber2.tab)  | tail -$(cat $OUTDIR/MA1929_1_final_1000bp_intersected_rowsNumber.tab) | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > $WRK/../data/RefPT-JASPAR-nonK562/MA1929_1_mm10_intersected_MPE-seq10min_164bp_category2.bed
cat $OUTDIR/MA1929_1_final_1000bp_intersected_ratio.TSV | head -$(cat $OUTDIR/MA1929_1_final_1000bp_intersected_rowsNumber3.tab)  | tail -$(cat $OUTDIR/MA1929_1_final_1000bp_intersected_rowsNumber.tab) | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > $WRK/../data/RefPT-JASPAR-nonK562/MA1929_1_mm10_intersected_MPE-seq10min_164bp_category3.bed
cat $OUTDIR/MA1929_1_final_1000bp_intersected_ratio.TSV | tail -$(cat $OUTDIR/MA1929_1_final_1000bp_intersected_rowsNumber.tab) | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > $WRK/../data/RefPT-JASPAR-nonK562/MA1929_1_mm10_intersected_MPE-seq10min_164bp_category4.bed

#expand bedfiles
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $WRK/../data/RefPT-JASPAR-nonK562/MA1929_1_mm10_intersected_MPE-seq10min_164bp_category1.bed $WRK/../data/RefPT-JASPAR-nonK562/MA1929_1_mm10_intersected_MPE-seq10min_164bp_category1_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $WRK/../data/RefPT-JASPAR-nonK562/MA1929_1_mm10_intersected_MPE-seq10min_164bp_category2.bed $WRK/../data/RefPT-JASPAR-nonK562/MA1929_1_mm10_intersected_MPE-seq10min_164bp_category2_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $WRK/../data/RefPT-JASPAR-nonK562/MA1929_1_mm10_intersected_MPE-seq10min_164bp_category3.bed $WRK/../data/RefPT-JASPAR-nonK562/MA1929_1_mm10_intersected_MPE-seq10min_164bp_category3_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $WRK/../data/RefPT-JASPAR-nonK562/MA1929_1_mm10_intersected_MPE-seq10min_164bp_category4.bed $WRK/../data/RefPT-JASPAR-nonK562/MA1929_1_mm10_intersected_MPE-seq10min_164bp_category4_1000bp.bed



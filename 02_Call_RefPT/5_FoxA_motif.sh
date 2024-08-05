#!/bin/bash

# get 

### CHANGE ME
WRK=/storage/group/bfp2/default/hxc585_HainingChen/
MOTIF=$WRK/Fox_NFIA_CTCF/data/RefPT-Motif
SCRIPTMANAGER=$WRK/2023_Chen_PIC3/bin/ScriptManager-v0.14.jar
GENOME=$WRK/2023_Chen_PIC3/data/hg19_files/hg19.fa
Nuc="$WRK/2024-Chen_Benzonase-ChIP-exo/03_RfMotif_sort_Nuc/K562_nuc_uHex_uTetra_1bp_sort_nostrand.bed"
# merge HepG2 FoxA1 sites
cat $MOTIF/1bp/FoxA1_HepG2_*_1bp.bed | bedtools sort -i | uniq | awk '{OFS="\t"}{FS="\t"}{print $1,$2,$3,$1"_"$2"_"$3,$5,$6,"HepG2"}' > $MOTIF/1bp/FoxA1_HepG2.bed
cat $MOTIF/1bp/FoxA1_K562_*_1bp.bed | bedtools sort -i | uniq | awk '{OFS="\t"}{FS="\t"}{print $1,$2,$3,$1"_"$2"_"$3,$5,$6,"K562"}' > $MOTIF/1bp/FoxA1_K562.bed

# test overlap
bedtools intersect -v -a $MOTIF/1bp/FoxA1_HepG2.bed -b $MOTIF/1bp/FoxA1_K562.bed | awk '{OFS="\t"}{FS="\t"}{print $1,$2,$3,$1"_"$2"_"$3,$5,$6,"HepG2_uniq"}'  > $MOTIF/1bp/FoxA1_uniq_HepG2.bed

wc -l $MOTIF/1bp/FoxA1_K562.bed
wc -l $MOTIF/1bp/FoxA1_uniq_HepG2.bed
1675 /storage/group/bfp2/default/hxc585_HainingChen//Fox_NFIA_CTCF/data/RefPT-Motif/1bp/FoxA1_K562.bed
2939 /storage/group/bfp2/default/hxc585_HainingChen//Fox_NFIA_CTCF/data/RefPT-Motif/1bp/FoxA1_uniq_HepG2.bed

bedtools sort -i $MOTIF/1bp/FoxA1_K562.bed | uniq | bedtools closest -a - -b $Nuc -d -D a | sort -k14,14n > $MOTIF/1bp/FoxA1_K562_NucSort.bed
bedtools sort -i $MOTIF/1bp/FoxA1_uniq_HepG2.bed | uniq | bedtools closest -a - -b $Nuc -d -D a | sort -k14,14n > $MOTIF/1bp/FoxA1_uniq_HepG2_NucSort.bed

cat $MOTIF/1bp/FoxA1_K562_NucSort.bed $MOTIF/1bp/FoxA1_uniq_HepG2_NucSort.bed > $MOTIF/1bp/FoxA1_all.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/1bp/FoxA1_all.bed -o $MOTIF/1000bp/FoxA1_all_1000bp.bed

awk -v motif="$MOTIF" '{
    if (($14 >= -73) && ($14 <= 73)) {
        print $0 > motif "/1bp/FoxA1_K562_NucSort-OVERLAP.bed"
    } else {
        print $0 > motif "/1bp/FoxA1_K562_NucSort-NFR.bed"
    }
}' "$MOTIF/1bp/FoxA1_K562_NucSort.bed"

wc -l $MOTIF/1bp/FoxA1_*_NucSort-OVERLAP.bed
wc -l $MOTIF/1bp/FoxA1_*_NucSort-NFR.bed

702 /storage/group/bfp2/default/hxc585_HainingChen//Fox_NFIA_CTCF/data/RefPT-Motif/1bp/FoxA1_K562_NucSort-OVERLAP.bed
1372 /storage/group/bfp2/default/hxc585_HainingChen//Fox_NFIA_CTCF/data/RefPT-Motif/1bp/FoxA1_uniq_K562_NucSort-NFR.bed

java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/1bp/FoxA1_K562_NucSort-OVERLAP.bed -o $MOTIF/1000bp/FoxA1_K562_NucSort-OVERLAP_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/1bp/FoxA1_K562_NucSort-NFR.bed -o $MOTIF/1000bp/FoxA1_K562_NucSort-NFR_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/1bp/FoxA1_K562_NucSort.bed -o $MOTIF/1000bp/FoxA1_K562_NucSort_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/1bp/FoxA1_uniq_HepG2_NucSort.bed -o $MOTIF/1000bp/FoxA1_uniq_HepG2_NucSort_1000bp.bed

exit

mkdir -p FoxA1
HepG2FoxA1=$WRK/Fox_NFIA_CTCF/data/BAM/HepG2_FoxA1_BX_rep1_hg19.bam
HepG2FoxA2=$WRK/Fox_NFIA_CTCF/data/BAM/HepG2_FoxA2_BX_rep3_hg19.bam
# HepG2FoxA3=$WRK/Fox_NFIA_CTCF/data/BAM/HepG2_FoxA3_BX_rep3_hg19.bam 
# K562FoxA1=$WRK/Fox_NFIA_CTCF/data/BAM/K562_FoxA1_BX_rep1_hg19.bam
K562FoxA2=$WRK/Fox_NFIA_CTCF/data/BAM/K562_FoxA2_BX_rep1_hg19.bam
K562FoxA3=$WRK/Fox_NFIA_CTCF/data/BAM/K562_FoxA3_BX_rep1_hg19.bam
K562IgG=$WRK/Fox_NFIA_CTCF/data/BAM/K562_IgG_BX_merge.bam
Input=$WRK/Fox_NFIA_CTCF/data/BAM/K562_-_BI_rep1_hg19.bam

for file in $HepG2FoxA1 $HepG2FoxA2 $K562FoxA2 $K562FoxA3 $K562IgG ; do
  BAM=$(basename "$file" ".bam")
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup $MOTIF/1000bp/FoxA1_all_1000bp.bed "$file" --cpu 4 -M FoxA1/${BAM}_FoxA1_all_1000bp_read1
  java -jar "$SCRIPTMANAGER" figure-generation heatmap -p .95 --blue FoxA1/${BAM}_FoxA1_all_1000bp_read1_sense.cdt -o FoxA1/${BAM}_FoxA1_all_1000bp_read1_sense_treeview.png
  java -jar "$SCRIPTMANAGER" figure-generation heatmap -p .95 --red FoxA1/${BAM}_FoxA1_all_1000bp_read1_anti.cdt -o FoxA1/${BAM}_FoxA1_all_1000bp_read1_anti_treeview.png
  java -jar "$SCRIPTMANAGER" figure-generation merge-heatmap FoxA1/${BAM}_FoxA1_all_1000bp_read1_sense_treeview.png FoxA1/${BAM}_FoxA1_all_1000bp_read1_anti_treeview.png -o FoxA1/${BAM}_FoxA1_all_1000bp_merge_treeview.png
  java -jar $SCRIPTMANAGER figure-generation label-heatmap FoxA1/${BAM}_FoxA1_all_1000bp_merge_treeview.png -f 20 -l -500 -m 0 -r 500 -o FoxA1/${BAM}_FoxA1_all_1000bp_merge_treeview.svg
  rm FoxA1/${BAM}_FoxA1_*.cdt FoxA1/${BAM}_FoxA1_*.png
done

for file in $Input ; do
  BAM=$(basename "$file" ".bam")
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup $MOTIF/1000bp/FoxA1_all_1000bp.bed "$file" --cpu 4 -m -M FoxA1/${BAM}_FoxA1_all_1000bp_midpoint
  java -jar "$SCRIPTMANAGER" figure-generation heatmap -p .95 --black FoxA1/${BAM}_FoxA1_all_1000bp_midpoint_combined.cdt -o FoxA1/${BAM}_FoxA1_all_1000bp_midpoint_combined_treeview.png
  java -jar $SCRIPTMANAGER figure-generation label-heatmap FoxA1/${BAM}_FoxA1_all_1000bp_midpoint_combined_treeview.png -f 20 -l -500 -m 0 -r 500 -o FoxA1/${BAM}_FoxA1_all_1000bp_midpoint_combined_treeview.svg
  rm FoxA1/${BAM}_FoxA1_*.cdt FoxA1/${BAM}_FoxA1_*.png
done

for file in $K562FoxA2 $K562FoxA3 $K562IgG $Input ; do
  BAM=$(basename "$file" ".bam")
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup $MOTIF/1000bp/FoxA1_K562_Nuc_overlap_1000bp.bed "$file" --cpu 4 -1 -o FoxA1/${BAM}_FoxA1_K562_Nuc_overlap_1000bp_read1.out
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup $MOTIF/1000bp/FoxA1_K562_Nuc_away_1000bp.bed "$file" --cpu 4 -1 -o FoxA1/${BAM}_FoxA1_K562_Nuc_away_1000bp_read1.out
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup $MOTIF/1000bp/FoxA1_uniq_HepG2_Nuc_overlap_1000bp.bed "$file" --cpu 4 -1 -o FoxA1/${BAM}_FoxA1_uniq_HepG2_Nuc_overlap_1000bp_read1.out
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup $MOTIF/1000bp/FoxA1_uniq_HepG2_Nuc_away_1000bp.bed "$file" --cpu 4 -1 -o FoxA1/${BAM}_FoxA1_uniq_HepG2_Nuc_away_1000bp_read1.out
done

for file in $K562FoxA2 $K562FoxA3 $Input ; do
  BAM=$(basename "$file" ".bam")
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup $MOTIF/1000bp/FoxA1_K562_Nuc_overlap_1000bp.bed "$file" --cpu 4 -2 -o FoxA1/${BAM}_FoxA1_K562_Nuc_overlap_1000bp_read2.out
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup $MOTIF/1000bp/FoxA1_K562_Nuc_away_1000bp.bed "$file" --cpu 4 -2 -o FoxA1/${BAM}_FoxA1_K562_Nuc_away_1000bp_read2.out
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup $MOTIF/1000bp/FoxA1_uniq_HepG2_Nuc_overlap_1000bp.bed "$file" --cpu 4 -2 -o FoxA1/${BAM}_FoxA1_uniq_HepG2_Nuc_overlap_1000bp_read2.out
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup $MOTIF/1000bp/FoxA1_uniq_HepG2_Nuc_away_1000bp.bed "$file" --cpu 4 -2 -o FoxA1/${BAM}_FoxA1_uniq_HepG2_Nuc_away_1000bp_read2.out
done

for file in $K562FoxA2 $K562FoxA3 ; do
  BAM=$(basename "$file" ".bam")
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup $MOTIF/1000bp/FoxA1_K562_Nuc_overlap_1000bp.bed "$file" --cpu 4 -1 -n 100 -o FoxA1/${BAM}_FoxA1_K562_Nuc_overlap_1000bp_MIN100_read1.out
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup $MOTIF/1000bp/FoxA1_K562_Nuc_away_1000bp.bed "$file" --cpu 4 -1 -n 100 -o FoxA1/${BAM}_FoxA1_K562_Nuc_away_1000bp_MIN100_read1.out
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup $MOTIF/1000bp/FoxA1_uniq_HepG2_Nuc_overlap_1000bp.bed "$file" --cpu 4 -1 -n 100 -o FoxA1/${BAM}_FoxA1_uniq_HepG2_Nuc_overlap_1000bp_MIN100_read1.out
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup $MOTIF/1000bp/FoxA1_uniq_HepG2_Nuc_away_1000bp.bed "$file" --cpu 4 -1 -n 100 -o FoxA1/${BAM}_FoxA1_uniq_HepG2_Nuc_away_1000bp_MIN100_read1.out
done
WRK=/storage/group/bfp2/default/hxc585_HainingChen/
MOTIF=$WRK/Fox_NFIA_CTCF/data/RefPT-Motif
SCRIPTMANAGER=$WRK/2023_Chen_PIC3/bin/ScriptManager-v0.14.jar
GENOME=$WRK/2023_Chen_PIC3/data/hg19_files/hg19.fa
Nuc="$WRK/2024-Chen_Benzonase-ChIP-exo/03_RfMotif_sort_Nuc/K562_nuc_uHex_uTetra_1bp_sort_nostrand.bed"
K562FoxA2=$WRK/Fox_NFIA_CTCF/data/BAM/K562_FoxA2_BX_rep1_hg19.bam
K562FoxA3=$WRK/Fox_NFIA_CTCF/data/BAM/K562_FoxA3_BX_rep1_hg19.bam
for file in $K562FoxA2 $K562FoxA3  ; do
  BAM=$(basename "$file" ".bam")
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup $MOTIF/1000bp/FoxA1_K562_Nuc_overlap_1000bp.bed "$file" --cpu 4 -2 -n 100 -o FoxA1/${BAM}_FoxA1_K562_Nuc_overlap_1000bp_MIN100_read2.out
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup $MOTIF/1000bp/FoxA1_K562_Nuc_away_1000bp.bed "$file" --cpu 4 -2 -n 100 -o FoxA1/${BAM}_FoxA1_K562_Nuc_away_1000bp_MIN100_read2.out
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup $MOTIF/1000bp/FoxA1_uniq_HepG2_Nuc_overlap_1000bp.bed "$file" --cpu 4 -2 -n 100 -o FoxA1/${BAM}_FoxA1_uniq_HepG2_Nuc_overlap_1000bp_MIN100_read2.out
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup $MOTIF/1000bp/FoxA1_uniq_HepG2_Nuc_away_1000bp.bed "$file" --cpu 4 -2 -n 100 -o FoxA1/${BAM}_FoxA1_uniq_HepG2_Nuc_away_1000bp_MIN100_read2.out
done




#!/bin/bash

# get 
module load anaconda
source activate bioinfo
### CHANGE ME
WRK=/storage/group/bfp2/default/hxc585_HainingChen/
MOTIF=$WRK/Fox_NFIA_CTCF/data/RefPT-Motif
SCRIPTMANAGER=$WRK/2023_Chen_PIC3/bin/ScriptManager-v0.14.jar
GENOME=$WRK/2023_Chen_PIC3/data/hg19_files/hg19.fa
Nuc="$WRK/2024-Chen_Benzonase-ChIP-exo/03_RfMotif_sort_Nuc/K562_nuc_uHex_uTetra_1bp_sort_nostrand.bed"
# expand CTCF sites

wc -l $MOTIF/1bp/CTCF_Occupancy_1bp.bed
20818 /storage/group/bfp2/default/hxc585_HainingChen//Fox_NFIA_CTCF/data/RefPT-Motif/1bp/CTCF_Occupancy_1bp.bed
bedtools sort -i $MOTIF/1bp/scrambledCTCF_top10k_1bp.bed  | uniq | bedtools closest -a - -b $Nuc -d -D a | sort -k13,13n  > $MOTIF/1bp/scrambledCTCF_NucSort_1bp.bed


java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/1bp/CTCF_Occupancy_1bp.bed -o $MOTIF/1000bp/CTCF_Occupancy_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/1bp/scrambledCTCF_NucSort_1bp.bed -o $MOTIF/1000bp/scrambledCTCF_NucSort_1000bp.bed

exit

mkdir -p CTCF
CTCF1=$WRK/Fox_NFIA_CTCF/data/BAM/K562_CTCF_BX_rep1_hg19.bam
CTCF2=$WRK/Fox_NFIA_CTCF/data/BAM/K562_CTCF_BX_rep2_hg19.bam
SMC31=$WRK/Fox_NFIA_CTCF/data/BAM/K562_SMC3_BX_rep1_hg19.bam
SMC32=$WRK/Fox_NFIA_CTCF/data/BAM/K562_SMC3_BX_rep2_hg19.bam
RAD211=$WRK/Fox_NFIA_CTCF/data/BAM/K562_RAD21_BX_rep1_hg19.bam
RAD212=$WRK/Fox_NFIA_CTCF/data/BAM/K562_RAD21_BX_rep2_hg19.bam
K562IgG=$WRK/Fox_NFIA_CTCF/data/BAM/K562_IgG_BX_merge.bam
Input=$WRK/Fox_NFIA_CTCF/data/BAM/K562_-_BI_rep1_hg19.bam

for file in $CTCF1 $K562IgG ; do
  BAM=$(basename "$file" ".bam")
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup $MOTIF/1000bp/CTCF_Occupancy_1000bp.bed "$file" --cpu 4 -M CTCF/${BAM}_CTCF_Occupancy_1000bp_read1
  java -jar "$SCRIPTMANAGER" figure-generation heatmap -p .95 --blue CTCF/${BAM}_CTCF_Occupancy_1000bp_read1_sense.cdt -o CTCF/${BAM}_CTCF_Occupancy_1000bp_read1_sense_treeview.png
  java -jar "$SCRIPTMANAGER" figure-generation heatmap -p .95 --red CTCF/${BAM}_CTCF_Occupancy_1000bp_read1_anti.cdt -o CTCF/${BAM}_CTCF_Occupancy_1000bp_read1_anti_treeview.png
  java -jar "$SCRIPTMANAGER" figure-generation merge-heatmap CTCF/${BAM}_CTCF_Occupancy_1000bp_read1_sense_treeview.png CTCF/${BAM}_CTCF_Occupancy_1000bp_read1_anti_treeview.png -o CTCF/${BAM}_CTCF_Occupancy_1000bp_merge_treeview.png
  java -jar $SCRIPTMANAGER figure-generation label-heatmap CTCF/${BAM}_CTCF_Occupancy_1000bp_merge_treeview.png -f 20 -l -500 -m 0 -r 500 -o CTCF/${BAM}_CTCF_Occupancy_1000bp_merge_treeview.svg
  rm CTCF/${BAM}_CTCF_*.cdt CTCF/${BAM}_CTCF_*.png
done

for file in $Input ; do
  BAM=$(basename "$file" ".bam")
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup $MOTIF/1000bp/CTCF_Occupancy_1000bp.bed "$file" --cpu 4 -m -M CTCF/${BAM}_CTCF_Occupancy_1000bp_midpoint
  java -jar "$SCRIPTMANAGER" figure-generation heatmap -p .95 --black CTCF/${BAM}_CTCF_Occupancy_1000bp_midpoint_combined.cdt -o CTCF/${BAM}_CTCF_Occupancy_1000bp_midpoint_combined_treeview.png
  java -jar $SCRIPTMANAGER figure-generation label-heatmap CTCF/${BAM}_CTCF_Occupancy_1000bp_midpoint_combined_treeview.png -f 20 -l -500 -m 0 -r 500 -o CTCF/${BAM}_CTCF_Occupancy_1000bp_midpoint_combined_treeview.svg
  rm CTCF/${BAM}_CTCF_*.cdt CTCF/${BAM}_CTCF_*.png
done

for file in $CTCF1 $SMC31 $RAD211 $K562IgG ; do
  BAM=$(basename "$file" ".bam")
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup $MOTIF/1000bp/CTCF_Occupancy_1000bp.bed "$file" --cpu 4 -1 -o CTCF/${BAM}_CTCF_Occupancy_1000bp_read1.out
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup $MOTIF/1000bp/CTCF_Occupancy_1000bp.bed "$file" --cpu 4 -2 -o CTCF/${BAM}_CTCF_Occupancy_1000bp_read2.out
done

for file in $CTCF1 $SMC31 $RAD211 ; do
  BAM=$(basename "$file" ".bam")
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup $MOTIF/1000bp/CTCF_Occupancy_1000bp.bed "$file" --cpu 4 -1 -n 100 -o CTCF/${BAM}_CTCF_Occupancy_1000bp_MIN100_read1.out
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup $MOTIF/1000bp/CTCF_Occupancy_1000bp.bed "$file" --cpu 4 -2 -n 100 -o CTCF/${BAM}_CTCF_Occupancy_1000bp_MIN100_read2.out
done



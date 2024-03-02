#!/bin/bash

### CHANGE ME
WINDOW=100
WRK=/path/to/2023-Chen_Benzonase-ChIP-exo/02_Call_RefPT
###

# Dependencies
# - java
# - python

set -exo
module load anaconda3
source activate bioinfo

# Script shortcuts
SCRIPTMANAGER=$WRK/bin/ScriptManager-v0.14.jar

# Inputs and outputs
BLACKLIST=$WRK/data/hg19_files/hg19_exclude.bed
MOTIF=$WRK/data/RefPT-Motif

# find top 100 sites of WDR5_M1
head -n 100 FIMO/WDR5/WDR5_MOTIF1_sorted.bed | awk '{OFS="\t"} {print $1,$2,$3,$1"_"$2"_"$3,$6,$7}' >  $MOTIF/WDR5_Occupancy.bed
# expand 1bp and sort
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1 $MOTIF/WDR5_Occupancy.bed -o $MOTIF/1bp/WDR5_Occupancy_1bp.bed
bedtools sort -i $MOTIF/1bp/WDR5_Occupancy_1bp.bed > $MOTIF/1bp/WDR5_Occupancy_1bp_sort.bed
# find TSS close to WDR5_M1
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1 $WRK/data/BED/tss_all_k562.bed -o $WRK/data/BED/tss_all_k562_1bp.bed
bedtools sort -i $WRK/data/BED/tss_all_k562_1bp.bed | bedtools closest -a $MOTIF/1bp/WDR5_Occupancy_1bp_sort.bed -b - -d -D b -t first | awk '{print $7,$8,$9,$10,$11,$12,$1,$2,$3,$4,$5,$6,$13}' | sort -k13,13n > FIMO/WDR5/TSS_WDR5_sorted.bed
# Clean-up
rm $MOTIF/1bp/WDR5_Occupancy_1bp_sort.bed
# Find TSS that are close to WDR5 less than 500bp
awk '{if (($13 <= 500) && ($13 >= -500) && ($7 != ".")) print $0 > "../data/BED/TSS94_WDR5_sort.bed"}' FIMO/WDR5/TSS_WDR5_sorted.bed
rm FIMO/WDR5/TSS_WDR5_sorted.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/WDR5_Occupancy.bed -o $MOTIF/1000bp/WDR5_Occupancy_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 ../data/BED/TSS94_WDR5_sort.bed -o ../data/BED/TSS94_WDR5_sort_1000bp.bed

awk -v motif="$MOTIF" '{OFS="\t"} {print $7,$8,$9,$10,$11,$12,$1,$2,$3,$4,$5,$6,$13}' ../data/BED/TSS94_WDR5_sort.bed | \
awk -v motif="$MOTIF" '{if ($6 == $12) print $0 > motif "/WDR5_TSS_same.bed"; else print $0 > motif "/WDR5_TSS_oppo.bed" }'

awk '{if ($6 == $12) print $0 > "../data/BED/TSS_WDR5_same_1000bp.bed"; else print $0 > "../data/BED/TSS_WDR5_oppo_1000bp.bed" }' ../data/BED/TSS94_WDR5_sort_1000bp.bed

java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 32 $MOTIF/WDR5_TSS_same.bed -o $MOTIF/WDR5_TSS_same_32bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 32 $MOTIF/WDR5_TSS_oppo.bed -o $MOTIF/WDR5_TSS_oppo_32bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/WDR5_TSS_oppo.bed -o $MOTIF/1000bp/WDR5_TSS_oppo_1000bp.bed

java -jar $SCRIPTMANAGER peak-analysis peak-align-ref -o FIMO/WDR5/WDR5_TSS_same_32bp_TSS94_WDR5_sort_1000bp.cdt $MOTIF/WDR5_TSS_same_32bp.bed ../data/BED/TSS94_WDR5_sort_1000bp.bed

java -jar $SCRIPTMANAGER figure-generation heatmap -a 1 --blue FIMO/WDR5/WDR5_TSS_same_32bp_TSS94_WDR5_sort_1000bp.cdt -o FIMO/WDR5/WDR5_TSS_same_32bp_TSS94_WDR5_sort_1000bp.png

java -jar $SCRIPTMANAGER peak-analysis peak-align-ref -o FIMO/WDR5/WDR5_TSS_oppo_32bp_TSS94_WDR5_sort_1000bp.cdt $MOTIF/WDR5_TSS_oppo_32bp.bed ../data/BED/TSS94_WDR5_sort_1000bp.bed

java -jar $SCRIPTMANAGER figure-generation heatmap -a 1 --red FIMO/WDR5/WDR5_TSS_oppo_32bp_TSS94_WDR5_sort_1000bp.cdt -o FIMO/WDR5/WDR5_TSS_oppo_32bp_TSS94_WDR5_sort_1000bp.png

java -jar $SCRIPTMANAGER figure-generation merge-heatmap FIMO/WDR5/WDR5_TSS_same_32bp_TSS94_WDR5_sort_1000bp.png FIMO/WDR5/WDR5_TSS_oppo_32bp_TSS94_WDR5_sort_1000bp.png -o FIMO/WDR5/WDR5_TSS_merge.png

java -jar $SCRIPTMANAGER figure-generation label-heatmap FIMO/WDR5/WDR5_TSS_merge.png -f 20 -l -500 -m 0 -r 500 -o FIMO/WDR5/WDR5_TSS_merge.svg

# Clean-up

rm FIMO/WDR5/WDR5_TSS_oppo_32bp_TSS94_WDR5_sort_1000bp.cdt FIMO/WDR5/WDR5_TSS_same_32bp_TSS94_WDR5_sort_1000bp.cdt

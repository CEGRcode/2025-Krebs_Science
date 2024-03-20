#!/bin/bash

### CHANGE ME
WINDOW=100
WRK=/path/to/2024-Chen_Nature/02_Call_RefPT
WRK=/storage/home/owl5022/scratch/2024-Chen_Nature/02_Call_RefPT
###

# Dependencies
# - java
# - python
# - bedtools

set -exo
module load anaconda3
module load bedtools
source activate bx

# Script shortcuts
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.14.jar

# Inputs and outputs
BLACKLIST=$WRK/../data/hg19_files/hg19_exclude.bed
MOTIF=$WRK/../data/RefPT-Motif

# find top 100 sites of WDR5_M1
head -n 100 FIMO/WDR5/WDR5_MOTIF1_sorted.bed | awk '{OFS="\t"}{print $1,$2,$3,$1"_"$2"_"$3,$6,$7}' > $MOTIF/WDR5_Occupancy.bed

# Expand 1000bp and 1bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1 $MOTIF/WDR5_Occupancy.bed -o $MOTIF/1bp/WDR5_Occupancy_1bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/WDR5_Occupancy.bed -o $MOTIF/1000bp/WDR5_Occupancy_1000bp.bed

# Genomic sort for bedtools
bedtools sort -i $MOTIF/1bp/WDR5_Occupancy_1bp.bed > $MOTIF/1bp/WDR5_Occupancy_1bp_sort.bed
# find Annotated TSS or GRO-cap defined all TSS close to WDR5_M1 
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1 $WRK/data/RefPT-Other/tss_all_k562.bed -o $WRK/data/RefPT-Other/tss_all_k562_1bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1 $WRK/data/RefPT-Other/K562_CoPRO-expressed_Gene-refSeqTSS_1000bp.bed -o $WRK/data/RefPT-Other/K562_CoPRO-expressed_Gene-refSeqTSS_1bp.bed
bedtools sort -i $WRK/data/RefPT-Other/tss_all_k562_1bp.bed | bedtools closest -a $MOTIF/1bp/WDR5_Occupancy_1bp_sort.bed -b - -d -D b -t first | awk '{print $7,$8,$9,$10,$11,$12,$1,$2,$3,$4,$5,$6,$13}' | sort -k13,13n | awk '{if (($13 <= 500) && ($13 >= -500) && ($7 != ".")) print $0 > "../data/RefPT-Other/TSS94_WDR5_sort.bed"}'
bedtool sort -i $WRK/data/RefPT-Other/K562_CoPRO-expressed_Gene-refSeqTSS_1bp.bed | bedtools closest -a $MOTIF/1bp/WDR5_Occupancy_1bp_sort.bed -b - -d -D a -t first | sort -k13,13n | awk '{if (($13 <= 500) && ($13 >= -500) && ($2 != "-1")) print $0 > "../data/RefPT-Other/Annotated_TSS_WDR5_sort.bed"}'  

# Clean-up
rm $MOTIF/1bp/WDR5_Occupancy_1bp_sort.bed

# expand all WDR5 and related TSS
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 ../data/RefPT-Other/TSS94_WDR5_sort.bed -o ../data/RefPT-Other/TSS94_WDR5_sort_1000bp.bed

# find WDR5 orientation
awk -v motif="$MOTIF" '{OFS="\t"} {print $7,$8,$9,$10,$11,$12,$1,$2,$3,$4,$5,$6,$13}' ../data/RefPT-Other/TSS94_WDR5_sort.bed | \
awk -v motif="$MOTIF" '{if ($6 == $12) print $0 > motif "/WDR5_TSS_same.bed"; else print $0 > motif "/WDR5_TSS_oppo.bed" }'

# Selecting TSS with WDR5 same orientation
awk '{if ($6 == $12) print $0 > "../data/RefPT-Other/TSS_WDR5_same_1000bp.bed"; else print $0 > "../data/RefPT-Other/TSS_WDR5_oppo_1000bp.bed" }' ../data/RefPT-Other/TSS94_WDR5_sort_1000bp.bed

# Plot figure of WDR5 motif relative to TSS, orientation seperately
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 32 $MOTIF/WDR5_TSS_same.bed -o $MOTIF/WDR5_TSS_same_32bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 32 $MOTIF/WDR5_TSS_oppo.bed -o $MOTIF/WDR5_TSS_oppo_32bp.bed

java -jar $SCRIPTMANAGER peak-analysis peak-align-ref -o FIMO/WDR5/WDR5_TSS_same_32bp_TSS94_WDR5_sort_1000bp.cdt $MOTIF/WDR5_TSS_same_32bp.bed ../data/RefPT-Other/TSS94_WDR5_sort_1000bp.bed

java -jar $SCRIPTMANAGER figure-generation heatmap -a 1 --blue FIMO/WDR5/WDR5_TSS_same_32bp_TSS94_WDR5_sort_1000bp.cdt -o FIMO/WDR5/WDR5_TSS_same_32bp_TSS94_WDR5_sort_1000bp.png

java -jar $SCRIPTMANAGER peak-analysis peak-align-ref -o FIMO/WDR5/WDR5_TSS_oppo_32bp_TSS94_WDR5_sort_1000bp.cdt $MOTIF/WDR5_TSS_oppo_32bp.bed ../data/RefPT-Other/TSS94_WDR5_sort_1000bp.bed

java -jar $SCRIPTMANAGER figure-generation heatmap -a 1 --red FIMO/WDR5/WDR5_TSS_oppo_32bp_TSS94_WDR5_sort_1000bp.cdt -o FIMO/WDR5/WDR5_TSS_oppo_32bp_TSS94_WDR5_sort_1000bp.png

java -jar $SCRIPTMANAGER figure-generation merge-heatmap FIMO/WDR5/WDR5_TSS_same_32bp_TSS94_WDR5_sort_1000bp.png FIMO/WDR5/WDR5_TSS_oppo_32bp_TSS94_WDR5_sort_1000bp.png -o FIMO/WDR5/WDR5_TSS_merge.png

java -jar $SCRIPTMANAGER figure-generation label-heatmap FIMO/WDR5/WDR5_TSS_merge.png -f 20 -l -500 -m 0 -r 500 -o FIMO/WDR5/WDR5_TSS_merge.svg

# Clean-up

rm FIMO/WDR5/WDR5_TSS_oppo_32bp_TSS94_WDR5_sort_1000bp.cdt FIMO/WDR5/WDR5_TSS_same_32bp_TSS94_WDR5_sort_1000bp.cdt

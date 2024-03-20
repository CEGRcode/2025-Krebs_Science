#!/bin/bash

# Prep RefPT for Figure 4 and S3 (WDR and TSS)

### CHANGE ME
NSITES=100
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
OTHER=$WRK/../data/RefPT-Other

# =====Make WDR5_Occupancy RefPTs=====

# Report top 100 WDR5 sites for the "Occupancy" reference
head -n $NSITES FIMO/WDR5/WDR5_MOTIF1_sorted.bed | awk '{OFS="\t"}{print $1,$2,$3,$1"_"$2"_"$3,$5,$6}' > $MOTIF/WDR5_Occupancy.bed

# Expand 1000bp and 1bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1 $MOTIF/WDR5_Occupancy.bed -o $MOTIF/1bp/WDR5_Occupancy_1bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/WDR5_Occupancy.bed -o $MOTIF/1000bp/WDR5_Occupancy_1000bp.bed

# =====Make TSS_DistWDR5 RefPTs=====

# Expand 1bp window for precise TSS mark
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1 $OTHER/tss_all_k562.bed -o GROcap_1bp.bed # GROcap
#java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1 $OTHER/K562_CoPRO-expressed_Gene-refSeqTSS_1000bp.bed -o CoPRO_1bp.bed # CoPRO and annotated genes

# Genomic sort for bedtools closest
bedtools sort -i $MOTIF/1bp/WDR5_Occupancy_1bp.bed > WDR5.bed
bedtools sort -i GROcap_1bp.bed > GROcap.bed
#bedtools sort -i CoPRO_1bp.bed > CoPRO.bed # pick one of these...

# Get closest TSS annotation (within -500 to +500)
bedtools closest -a WDR5.bed -b GROcap.bed -d -D b -t first \
    | awk '{OFS="\t"}{if ($13 <= 500 && $13 >= -500) print $7,$8,$9,$10,$11,$12,$1,$2,$3,$4,$5,$6,$13 }' \
    | sort -nk13,13 \
    > GROcap_DistWDR5_Filter-closest.bed

cp GROcap_DistWDR5_Filter-closest.bed $OTHER/TSS_DistWDR5.bed

# Expand 1000bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $OTHER/TSS_DistWDR5.bed -o $OTHER/TSS_DistWDR5_1000bp.bed

# Make same-strand filtered TSS RefPTs
awk '{OFS="\t"}{if ($6 == $12) print $1,$2,$3,$4,$5,$6}' $OTHER/TSS_DistWDR5_Filter-SameStrand_1000bp.bed

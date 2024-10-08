#!/bin/bash

# get 
module load anaconda
source activate bioinfo
### CHANGE ME
WRK=/storage/group/bfp2/default/hxc585_HainingChen/Fox_NFIA_CTCF/
MOTIF=$WRK/data/RefPT-Motif/
FIMO=$WRK/02_Call_RefPT/FIMO
SCRIPTMANAGER=$WRK/bin/ScriptManager-v0.15.jar
GENOME=$WRK/data/hg38_files/hg38.fa
Nuc="$WRK/Motif_analysis/K562_Nucpeak_hg38liftover_sort.bed"
# expand CTCF sites

mkdir -p $WRK/Motif_analysis/CTCF

cd $WRK/Motif_analysis/CTCF

wc -l $FIMO/CTCF/CTCF_Occupancy_1bp.bed

awk '{OFS="\t"} {print $1,$2,$3,$4,$7,$6,"K562_CTCF"}' $FIMO/CTCF/CTCF_Occupancy_1bp.bed  > CTCF_Occuoancy.bed 

java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 CTCF_Occuoancy.bed -o $MOTIF/1000bp/CTCF_Occuoancy_1000bp.bed

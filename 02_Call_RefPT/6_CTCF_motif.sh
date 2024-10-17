#!/bin/bash

# get 
module load anaconda
source activate bioinfo
### CHANGE ME
WRK=/Path/to/Title
MOTIF=$WRK/data/RefPT-Motif/
FIMO=$WRK/02_Call_RefPT/FIMO
SCRIPTMANAGER=$WRK/bin/ScriptManager-v0.15.jar
GENOME=$WRK/data/hg38_files/hg38.fa
Nuc="$WRK/Motif_analysis/K562_Nucpeak_hg38liftover_sort.bed"
# expand CTCF sites

mkdir -p $WRK/Motif_analysis/CTCF

cd $WRK/Motif_analysis/CTCF

awk '{OFS="\t"} {print $1,$2,$3,$4,$7,$6,"K562_CTCF"}' $FIMO/CTCF/CTCF_Occupancy_1bp.bed  | awk '
{
    if ($1 !~/alt/ && $1 !~/random/ && $1 !~/Un/ ) {
        print $0 > "CTCF_Occuoancy_1bp.bed";
    } 
}'

wc -l CTCF_Occupancy_1bp.bed


java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 CTCF_Occuoancy.bed -o $MOTIF/1000bp/CTCF_Occuoancy_1000bp.bed

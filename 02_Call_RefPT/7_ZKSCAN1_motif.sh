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
Genome=$WRK/data/hg38_files/hg38.info.txt

ZKSCAN1=$WRK/Motif_analysis/ZKSCAN1

mkdir -p $ZKSCAN1FILE
cd $ZKSCAN1

# determine the nucleosome engagment direction on each NFIA sites

awk '
{
    if ($1 !~/alt/ && $1 !~/random/ && $1 !~/Un/ ) {
        print $0 > "ZKSCAN1_Occupancy_1bp.bed";
    } 
}' $FIMO/ZKSCAN1/ZKSCAN1_Occupancy_1bp.bed

awk '
{
    if ($6 == "+") {
        print $0 > "ZKSCAN1_Occupancy_1bp+.bed";
    } else {
        print $0 > "ZKSCAN1_Occupancy_1bp-.bed";
    }
}' ZKSCAN1_Occupancy_1bp.bed

# make ZKSCAN1 same orientation and 0 point as jordan's
# For the + strand
awk '{OFS="\t"; print $1, $2, $3, $4, $5, "-", $7}' ZKSCAN1_Occupancy_1bp+.bed | bedtools shift -i - -g $Genome -p 2 -m -2 > ZKSCAN1_Occupancy_1bp-_shift2.bed

# For the - strand
awk '{OFS="\t"; print $1, $2, $3, $4, $5, "+", $7}' ZKSCAN1_Occupancy_1bp-.bed | bedtools shift -i - -g $Genome -p 2 -m -2 > ZKSCAN1_Occupancy_1bp+_shift2.bed
cat ZKSCAN1_Occupancy_1bp-_shift2.bed ZKSCAN1_Occupancy_1bp+_shift2.bed | bedtools sort -i | uniq | sort -k7,7nr > ZKSCAN1_Occupancy_flip_shift2.bed

rm ZKSCAN1_Occupancy_1bp-_shift2.bed ZKSCAN1_Occupancy_1bp+_shift2.bed ZKSCAN1_Occupancy_1bp+.bed ZKSCAN1_Occupancy_1bp-.bed
wc -l ZKSCAN1_Occupancy_1bp.bed
# expand
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 ZKSCAN1_Occupancy_flip_shift2.bed  -o $MOTIF/1000bp/ZKSCAN1_Occupancy_flip_shift2_1000bp.bed




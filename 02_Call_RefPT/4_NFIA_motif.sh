#!/bin/bash

### CHANGE ME
WRK=/path/to/2024-Chen_Nature/02_Call_RefPT
WRK=/storage/home/owl5022/scratch/2024-Chen_Nature/02_Call_RefPT
###

set -exo
module load anaconda
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Inputs and outputs
MOTIF=../data/RefPT-Motif/
GENOME=../data/hg38_files/hg38.fa
GINFO=../data/hg38_files/hg38.info.txt

NUCLEOSOME=Motif_analysis/K562_Nucpeak_hg38liftover_sort.bed
NFIA1=../data/BAM/K562_NFIA_BX_rep1_hg38.bam

# Script shortcuts
SCRIPTMANAGER=../bin/ScriptManager-v0.15.jar

[ -d Motif_analysis/NFIA ] || mkdir Motif_analysis/NFIA
NFIA=Motif_analysis/NFIA

# TODO: Fix filepaths

# determine the nucleosome engagment direction on each NFIA sites
awk '
{
    if $1 !~/alt/ ) && ($1 !~/random/) && ($1 !~/Un/) {
        print $0 > "NFIA_Occupancy_1bp.bed";
    } 
} ' FIMO/NFIA/NFIA_Occupancy_1bp.bed

wc -l NFIA_Occupancy_1bp.bed
# 7965 NFIA_Occupancy_1bp.bed
bedtools shift -i NFIA_Occupancy_1bp.bed -g $GINFO -p 95 -m -95 > NFIA_Occupancy_down95bp.bed
bedtools shift -i NFIA_Occupancy_1bp.bed -g $GINFO -p -95 -m 95 > NFIA_Occupancy_up95bp.bed

mkdir -p $NFIA/SCORES

for file in NFIA_Occupancy_down95bp.bed NFIA_Occupancy_up95bp.bed ; do
  filename=$(basename $file ".bed")
  java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 150 $file -o ${filename}_150bp.bed
  java -jar $SCRIPTMANAGER read-analysis tag-pileup ${filename}_150bp.bed $NFIA1 -n 100 -1 --combined --cpu 4 -M SCORES/NFIA_${filename}_150bp_read1_MIN100
  rm ${filename}_150bp.bed
done

java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum SCORES/NFIA_NFIA_Occupancy_up95bp_150bp_read1_MIN100_combined.cdt -o SCORES/
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum SCORES/NFIA_NFIA_Occupancy_down95bp_150bp_read1_MIN100_combined.cdt -o SCORES/

cut -f 2 SCORES/NFIA_NFIA_Occupancy_down95bp_150bp_read1_MIN100_combined_SCORES.out \
    | paste SCORES/NFIA_NFIA_Occupancy_up95bp_150bp_read1_MIN100_combined_SCORES.out - \
    | tail -n +2 | cut -f 2-3 \
    | paste NFIA_Occupancy_1bp.bed - \
    awk '
{
    if ((($8 > $9 && $6 == "+") || ($8 > $9 && $6 == "-")) || ($8 == $9 && $6 == "+")) {
        print $0 > "NFIA_Nuc_upstream_temp.bed";
    } else {
        print $0 > "NFIA_Nuc_downstream_temp.bed";
    }
}
'

wc -l NFIA_Nuc_upstream_temp.bed
# 4093 NFIA_Nuc_upstream_temp.bed
wc -l NFIA_Nuc_downstream_temp.bed
# 3873 NFIA_Nuc_downstream_temp.bed

awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,"engageNuddown"}' NFIA_Nuc_downstream_temp.bed > NFIA_engagedNuc_downstream.bed
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,"engageNudup"}' NFIA_Nuc_upstream_temp.bed > NFIA_engagedNuc_upstream.bed
rm NFIA_Nuc_downstream_temp.bed NFIA_Nuc_upstream_temp.bed 
rm NFIA_Occupancy_*95bp.bed 

## convert NFIA all to engage Nuc downstream

awk '{ if ($6 == "+") print $0 > NFIA_Nuc_upstream_+.bed; else print $0 > NFIA_Nuc_upstream_-.bed }' NFIA_engagedNuc_upstream.bed
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,"+",$7,$9,$8,$10}' NFIA_Nuc_upstream_-.bed > NFIA_Nuc_upstream_-_convert.bed
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,"-",$7,$9,$8,$10}' NFIA_Nuc_upstream_+.bed > NFIA_Nuc_upstream_+_convert.bed

rm NFIA_Nuc_upstream_-.bed NFIA_Nuc_upstream_+.bed

cat NFIA_Nuc_upstream_-_convert.bed NFIA_Nuc_upstream_+_convert.bed \
    | bedtools shift -i - -p 1 -m -1 -g $GINFO \
    | cat - NFIA_engagedNuc_downstream.bed \
    | sort -k7,7nr > NFIA_downNuc.bed

rm NFIA_Nuc_upstream_-_convert.bed NFIA_Nuc_upstream_+_convert.bed


# Sort orientated NFIA by distance to closest nucleosome, seperate by up/down/overlap

bedtools sort -i NFIA_downNuc.bed | uniq | bedtools closest -a - -b $Nuc -d -D a -t first | sort -k17,17n > NFIA_engageNucdown_NucSort.bed 

awk '{
   if ($17 < -73 ) {
        print $0 >  "NFIA_NucSort-UPSTREAM.bed"
    } else if ($17 > 73 ) {
        print $0 >  "NFIA_NucSort-DOWNSTREAM.bed"
    } else {
      print $0 >  "NFIA_NucSort-OVERLAP.bed"
    }
}' NFIA_engageNucdown_NucSort.bed 


wc -l NFIA_NucSort-*.bed
# 1563 NFIA_NucSort-DOWNSTREAM.bed
# 5248 NFIA_NucSort-OVERLAP.bed
# 1155 NFIA_NucSort-UPSTREAM.bed

## random strand
cut -f 3-152 SCORES/NFIA_NFIA_Occupancy_up95bp_150bp_read1_MIN100_combined.cdt \
    | paste SCORES/NFIA_NFIA_Occupancy_down95bp_150bp_read1_MIN100_combined.cdt - \
    > SCORES/NFIA_NFIA_Occupancy_down95bp_150bp_up95bp_150bp_read1_combined.cdt

python ../shuffle_script.py SCORES/NFIA_NFIA_Occupancy_down95bp_150bp_up95bp_150bp_read1_combined.cdt SCORES/NFIA_NFIA_Occupancy_down95bp_150bp_up95bp_150bp_read1_combined_shuffled.cdt

cut -f 1-152 SCORES/NFIA_NFIA_Occupancy_down95bp_150bp_up95bp_150bp_read1_combined_shuffled.cdt > SCORES/NFIA_NFIA_Occupancy_down95bp_150bp_read1_MIN100_combined_shuffled.cdt
cut -f 1-2,153-303 SCORES/NFIA_NFIA_Occupancy_down95bp_150bp_up95bp_150bp_read1_combined_shuffled.cdt > SCORES/NFIA_NFIA_Occupancy_up95bp_150bp_read1_MIN100_combined_shuffled.cdt

java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum SCORES/NFIA_NFIA_Occupancy_down95bp_150bp_read1_MIN100_combined_shuffled.cdt -o $NFIA/SCORES
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum SCORES/NFIA_NFIA_Occupancy_up95bp_150bp_read1_MIN100_combined_shuffled.cdt -o $NFIA/SCORES

cut -f 2 SCORES/NFIA_NFIA_Occupancy_down95bp_150bp_read1_MIN100_combined_shuffled_SCORES.out \
    | paste SCORES/NFIA_NFIA_Occupancy_up95bp_150bp_read1_MIN100_combined_shuffled_SCORES.out - \
    | tail -n +2 | cut -f 2-3 \
    | paste NFIA_Occupancy_1bp.bed - \
    | awk '
{
    if ((($8 > $9 && $6 == "+") || ($8 > $9 && $6 == "-")) || ($8 == $9 && $6 == "+")) {
        print $0 > "NFIA_Nuc_upstream_shuffled_temp.bed";
    } else {
        print $0 > "NFIA_Nuc_downstream_shuffled_temp.bed";
    }

}'

wc -l NFIA_Nuc_upstream_shuffled_temp.bed
wc -l NFIA_Nuc_downstream_shuffled_temp.bed


awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,"engageNuddown"}' NFIA_Nuc_downstream_shuffled_temp.bed > NFIA_engagedNuc_shuffled_downstream.bed
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,"engageNudup"}' NFIA_Nuc_upstream_shuffled_temp.bed > NFIA_engagedNuc_shuffled_upstream.bed
rm NFIA_Nuc_downstream_shuffled_temp.bed NFIA_Nuc_upstream_shuffled_temp.bed 
## convert NFIA with randomized socre all to engage Nuc downstream

awk '{ if ($6 == "+") print $0 >  "NFIA_Nuc_upstream_shuffled_+.bed"; else print $0 >  "NFIA_Nuc_upstream_shuffled_-.bed" }' NFIA_engagedNuc_shuffled_upstream.bed
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,"+",$7,$9,$8,$10}' NFIA_Nuc_upstream_shuffled_-.bed >  NFIA_Nuc_upstream_shuffled_-_convert.bed
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,"-",$7,$9,$8,$10}'  NFIA_Nuc_upstream_shuffled_+.bed >  NFIA_Nuc_upstream_shuffled_+_convert.bed

rm NFIA_Nuc_upstream_shuffled_-.bed NFIA_Nuc_upstream_shuffled_+.bed

cat NFIA_Nuc_upstream_shuffled_-_convert.bed NFIA_Nuc_upstream_shuffled_+_convert.bed \
    | bedtools shift -i - -p 1 -m -1 -g $GINFO \
    | cat - NFIA_engagedNuc_shuffled_downstream.bed \
    | sort -k7,7nr > NFIA_shuffled_downNuc.bed

rm NFIA_Nuc_upstream_shuffled_-_convert.bed NFIA_Nuc_upstream_shuffled_+_convert.bed

bedtools shift -i NFIA_shuffled_downNuc.bed -g $GINFO -p 95 -m -95 > NFIA_shuffled_downNuc_down95.bed
bedtools shift -i NFIA_shuffled_downNuc.bed -g $GINFO -p -95 -m 95 > NFIA_shuffled_downNuc_up95.bed

for file in NFIA_shuffled_downNuc_down95.bed NFIA_shuffled_downNuc_up95.bed ; do
  filename=$(basename $file ".bed")
  java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 150 $file -o ${filename}_150bp.bed
  java -jar $SCRIPTMANAGER read-analysis tag-pileup ${filename}_150bp.bed $NFIA1 -n 100 -1 --combined --cpu 4 -M SCORES/NFIA_${filename}_150bp_read1_MIN100
  java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum SCORES/NFIA_${filename}_150bp_read1_MIN100_combined.cdt -o SCORES/
  rm ${filename}_150bp.bed
done

cut -f 2 SCORES/NFIA_NFIA_shuffled_downNuc_down95_150bp_read1_MIN100_combined_SCORES.out \
    | paste SCORES/NFIA_NFIA_shuffled_downNuc_up95_150bp_read1_MIN100_combined_SCORES.out - \
    | tail -n +2 | cut -f 2-3 \
    | paste NFIA_shuffled_downNuc.bed - \
    > NFIA_shuffled_downNuc_SCORES.bed
rm NFIA_shuffled_downNuc_down95.bed NFIA_shuffled_downNuc_up95.bed
## shift Orientated NFIA to downstream 125bp for better see 10bp pattern in heatmap
bedtools shift -i NFIA_downNuc.bed -p 125 -m -125 -g $GINFO > NFIA_downNuc_down125.bed
## take Q4 of BX NFIA sites

tail -n 1991 NFIA_downNuc.bed > NFIA_downNuc_4.bed

java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 NFIA_engageNucdown_NucSort.bed -o $MOTIF/1000bp/NFIA_engageNucdown_NucSort_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 500  NFIA_downNuc.bed -o $MOTIF/1000bp/NFIA_downNuc_500bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 250  NFIA_downNuc_down125.bed -o $MOTIF/1000bp/NFIA_downNuc_down125_250bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 NFIA_NucSort-UPSTREAM.bed -o $MOTIF/1000bp/NFIA_NucSort-UPSTREAM_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 NFIA_NucSort-DOWNSTREAM.bed -o $MOTIF/1000bp/NFIA_NucSort-DOWNSTREAM_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 NFIA_NucSort-OVERLAP.bed -o $MOTIF/1000bp/NFIA_NucSort-OVERLAP_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 NFIA_downNuc_4.bed -o $MOTIF/1000bp/NFIA_downNuc_4_1000bp.bed




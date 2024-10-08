#!/bin/bash

# get 
module load anaconda
source activate bioinfo
### CHANGE ME
WRK=/storage/group/bfp2/default/hxc585_HainingChen/Fox_NFIA_CTCF
MOTIF=$WRK/data/RefPT-Motif/
FIMO=$WRK/02_Call_RefPT/FIMO
SCRIPTMANAGER=$WRK/bin/ScriptManager-v0.15.jar
GENOME=$WRK/data/hg38_files/hg38.fa
Genome=$WRK/data//hg38_files/hg38.info.txt
Nuc="$WRK/Motif_analysis/K562_Nucpeak_hg38liftover_sort.bed"

mkdir -p $WRK/Motif_analysis/FOXA
cd $WRK/Motif_analysis/FOXA

# merge HepG2 FoxA1 sites

cat $FIMO/FOXA2/FOXA_HepG2_*_Occupancy_1bp.bed | \
bedtools sort -i - | \
uniq | \
awk '($1 !~ /alt|random|chrUn/){OFS="\t"; print $1, $2, $3, $1"_"$2"_"$3, $5, $6, "HepG2"}' | bedtools sort -i | uniq > FOXA_HepG2.bed

cat $FIMO/FOXA2/FOXA_K562_*_Occupancy_1bp.bed | \
bedtools sort -i - | \
uniq | \
awk '($1 !~ /alt|random|chrUn/){OFS="\t"; print $1, $2, $3, $1"_"$2"_"$3, $5, $6, "K562"}' | bedtools sort -i | uniq > FOXA_K562.bed
# test overlap
bedtools intersect -v -a FOXA_HepG2.bed -b FOXA_K562.bed | awk '{OFS="\t"}{FS="\t"}{print $1,$2,$3,$4,$5,$6,"HepG2_uniq"}'  > FOXA_uniq_HepG2.bed
bedtools intersect -u -a FOXA_HepG2.bed -b FOXA_K562.bed | awk '{OFS="\t"}{FS="\t"}{print $1,$2,$3,$4,$5,$6,"K562_HepG2_overlap"}'  > FOXA_K562_HepG2_overlap.bed
bedtools intersect -v -a FOXA_K562.bed  -b FOXA_HepG2.bed | awk '{OFS="\t"}{FS="\t"}{print $1,$2,$3,$4,$5,$6,"K562_uniq"}'  > FOXA_uniq_K562.bed
wc -l FOXA_K562.bed
wc -l FOXA_HepG2.bed
wc -l FOXA_uniq_HepG2.bed
wc -l FOXA_K562_HepG2_overlap.bed
wc -l  FOXA_uniq_K562.bed
2642 FOXA_uniq_HepG2.bed
614 FOXA_uniq_K562.bed
1406 FOXA_K562.bed
3434 FOXA_HepG2.bed
792 FOXA_K562_HepG2_overlap.bed

rm FOXA_K562.bed
cat FOXA_uniq_K562.bed FOXA_K562_HepG2_overlap.bed | bedtools sort -i | uniq >  FOXA_K562.bed
## test if FoxA1 interact with Nucleosome downstream

HepG2FoxA1=$WRK/data/BAM/HepG2_FOXA1_BX_rep1_hg38.bam
K562FoxA1=$WRK/data/BAM/K562_FOXA1_BX_rep1_hg38.bam

bedtools shift -i FOXA_K562.bed -g $Genome -p 50 -m -50 > FOXA_K562_down50.bed
bedtools shift -i FOXA_K562.bed -g $Genome -p 150 -m -150 > FOXA_K562_down150.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 FOXA_K562_down50.bed -o FOXA_K562_down50_100bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 FOXA_K562_down150.bed -o FOXA_K562_down150_100bp.bed

mkdir -p SCORES

java -jar "$SCRIPTMANAGER" read-analysis tag-pileup FOXA_K562_down50_100bp.bed $K562FoxA1 -2 --cpu 4 -M SCORES/K562FoxA1_FOXA_K562_down50_100bp_read2
java -jar "$SCRIPTMANAGER" read-analysis tag-pileup FOXA_K562_down150_100bp.bed $K562FoxA1 -2 --cpu 4 -M  SCORES/K562FoxA1_FOXA_K562_down150_100bp_read2

java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum SCORES/K562FoxA1_FOXA_K562_down50_100bp_read2_anti.cdt -o SCORES/
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum SCORES/K562FoxA1_FOXA_K562_down150_100bp_read2_anti.cdt -o SCORES/
tail -n +2 SCORES/K562FoxA1_FOXA_K562_down50_100bp_read2_anti_SCORES.out | cut -f 2 | paste FOXA_K562.bed - > FOXA_K562_read2_down100.bed
tail -n +2 SCORES/K562FoxA1_FOXA_K562_down150_100bp_read2_anti_SCORES.out | cut -f 2 | paste FOXA_K562_read2_down100.bed -  | sort -k9,9nr > FOXA_K562_NucengageSort.bed
rm FOXA_K562_down50.bed FOXA_K562_down150.bed FOXA_K562_down50_100bp.bed FOXA_K562_down150_100bp.bed
  

bedtools shift -i FOXA_uniq_HepG2.bed -g $Genome -p 50 -m -50 > FOXA_uniq_HepG2_down50.bed
bedtools shift -i FOXA_uniq_HepG2.bed -g $Genome -p 150 -m -150 > FOXA_uniq_HepG2_down150.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 FOXA_uniq_HepG2_down50.bed -o FOXA_uniq_HepG2_down50_100bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 FOXA_uniq_HepG2_down150.bed -o FOXA_uniq_HepG2_down150_100bp.bed
java -jar "$SCRIPTMANAGER" read-analysis tag-pileup FOXA_uniq_HepG2_down50_100bp.bed $HepG2FoxA1 -2 --cpu 4 -M SCORES/HepG2FoxA1_FOXA_uniq_HepG2_down50_100bp_read2
java -jar "$SCRIPTMANAGER" read-analysis tag-pileup FOXA_uniq_HepG2_down150_100bp.bed $HepG2FoxA1 -2 --cpu 4 -M SCORES/HepG2FoxA1_FOXA_uniq_HepG2_down150_100bp_read2
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum SCORES/HepG2FoxA1_FOXA_uniq_HepG2_down50_100bp_read2_anti.cdt -o SCORES/
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum SCORES/HepG2FoxA1_FOXA_uniq_HepG2_down150_100bp_read2_anti.cdt -o SCORES/
tail -n +2 SCORES/HepG2FoxA1_FOXA_uniq_HepG2_down50_100bp_read2_anti_SCORES.out | cut -f 2 | paste FOXA_uniq_HepG2.bed - > FOXA_uniq_HepG2_read2_down100.bed
tail -n +2 SCORES/HepG2FoxA1_FOXA_uniq_HepG2_down150_100bp_read2_anti_SCORES.out | cut -f 2 | paste FOXA_uniq_HepG2_read2_down100.bed - | sort -k9,9nr > FOXA_uniq_HepG2_NucengageSort.bed
rm FOXA_uniq_HepG2_down50.bed FOXA_uniq_HepG2_down150.bed FOXA_uniq_HepG2_down50_100bp.bed FOXA_uniq_HepG2_down150_100bp.bed FOXA_uniq_HepG2_read2_down100.bed


bedtools sort -i FOXA_K562_NucengageSort.bed | uniq | bedtools closest -a - -b $Nuc -d -D a | sort -k16,16n > FOXA_K562_NucSort.bed
bedtools sort -i FOXA_uniq_HepG2_NucengageSort.bed | uniq | bedtools closest -a - -b $Nuc -d -D a | sort -k16,16n > FOXA_uniq_HepG2_NucSort.bed


awk '{
    if ($16 >= -73 && $16 <= 73) {
        print $0 >> "FOXA_K562_NucSort-OVERLAP.bed"
    } else {
        print $0 >> "FOXA_K562_NucSort-NFR.bed"
    }
}' FOXA_K562_NucSort.bed


awk '{
    if ($16 >= -73 && $16 <= 73) {
        print $0 >> "FOXA_uniq_HepG2_NucSort-OVERLAP.bed"
    } else {
        print $0 >> "FOXA_uniq_HepG2_NucSort-NFR.bed"
    }
}' FOXA_uniq_HepG2_NucSort.bed

wc -l *NucSort-*.bed
571 FOXA_K562_NucSort-NFR.bed
 837 FOXA_K562_NucSort-OVERLAP.bed
1164 FOXA_uniq_HepG2_NucSort-NFR.bed
1480 FOXA_uniq_HepG2_NucSort-OVERLAP.bed

cat FOXA_K562_NucSort.bed FOXA_uniq_HepG2_NucSort.bed > FOXA_all.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 FOXA_all.bed -o $MOTIF/1000bp/FOXA_all_1000bp.bed

java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 500 FOXA_K562_NucengageSort.bed -o $MOTIF/1000bp/FOXA_K562_NucengageSort_500bp.bed

head -n 500 FOXA_K562_NucengageSort.bed > FOXA_K562_Nucengage.bed
tail -n +501 FOXA_K562_NucengageSort.bed > FOXA_K562_NoNuc.bed

head -n 900 FOXA_uniq_HepG2_NucengageSort.bed > FOXA_HepG2_Nucengage.bed
tail -n +901 FOXA_uniq_HepG2_NucengageSort.bed > FOXA_HepG2_NoNuc.bed


java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 FOXA_K562_NucSort.bed -o $MOTIF/1000bp/FOXA_K562_NucSort_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 FOXA_uniq_HepG2_NucSort.bed -o $MOTIF/1000bp/FOXA_uniq_HepG2_NucSort_1000bp.bed

java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 FOXA_K562_NucSort-OVERLAP.bed -o $MOTIF/1000bp/FOXA_K562_NucSort-OVERLAP_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 FOXA_K562_NucSort-NFR.bed -o $MOTIF/1000bp/FOXA_K562_NucSort-NFR_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 FOXA_uniq_HepG2_NucSort-OVERLAP.bed -o $MOTIF/1000bp/FOXA_uniq_HepG2_NucSort-OVERLAP_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 FOXA_uniq_HepG2_NucSort-NFR.bed -o $MOTIF/1000bp/FOXA_uniq_HepG2_NucSort-NFR_1000bp.bed

java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 FOXA_K562_Nucengage.bed -o $MOTIF/1000bp/FOXA_K562_Nucengage_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 FOXA_K562_NoNuc.bed -o $MOTIF/1000bp/FOXA_K562_NoNuc_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 FOXA_HepG2_Nucengage.bed -o $MOTIF/1000bp/FOXA_HepG2_Nucengage_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 FOXA_HepG2_NoNuc.bed -o $MOTIF/1000bp/FOXA_HepG2_NoNuc_1000bp.bed




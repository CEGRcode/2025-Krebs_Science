#!/bin/bash

# Build FOXA motif-related reference points

# data/RefPT-Motif
#   |--FOXA_LABEL-K562_SORT-NucleosomeEngagement.bed                                (see 500bp)
#   |--FOXA_LABEL-K562_SORT-NucleosomeEngagement_GROUP-Engaged.bed                  (see 500bp)
#   |--FOXA_LABEL-K562_SORT-NucleosomeEngagement_GROUP-LessEngaged.bed              (see 500bp)
#   |--FOXA_LABEL-uHepG2_SORT-NucleosomeEngagement.bed                              (see 500bp)
#   |--FOXA_LABEL-uHepG2_SORT-NucleosomeEngagement_GROUP-Engaged.bed                (see 500bp)
#   |--FOXA_LABEL-uHepG2_SORT-NucleosomeEngagement_GROUP-LessEngaged.bed            (see 500bp)
#   |--FOXA_LABEL-K562_SORT-ClosestDyad.bed                                         (see 1000bp)
#   |--FOXA_LABEL-K562_SORT-ClosestDayd_GROUP-NoOverlap.bed                         (see 1000bp)
#   |--FOXA_LABEL-K562_SORT-ClosestDyad_GROUP-NFR.bed                               (see 1000bp)
#   |--FOXA_LABEL-uHepG2_SORT-ClosestDyad.bed                                       (see 1000bp)
#   |--FOXA_SORT-ClosestDyad_STACK-K562-uHepG2.bed                                  (see 1000bp)
#   |--500bp
#      |--FOXA_LABEL-K562_SORT-NucleosomeEngagement_500bp.bed                       (FOXA_K562_NucengageSort_500bp.bed)
#      |--FOXA_LABEL-K562_SORT-NucleosomeEngagement_GROUP-Engaged_500bp.bed         (FOXA_K562_Nucengage_1000bp.bed)
#      |--FOXA_LABEL-K562_SORT-NucleosomeEngagement_GROUP-LessEngaged_500bp.bed     (FOXA_K562_NoNuc_1000bp.bed)
#      |--FOXA_LABEL-uHepG2_SORT-NucleosomeEngagement_500bp.bed                     (FOXA_uniq_HepG2_NucengageSort_500bp.bed)
#      |--FOXA_LABEL-uHepG2_SORT-NucleosomeEngagement_GROUP-Engaged_500bp.bed       (FOXA_HepG2_Nucengage_1000bp.bed)
#      |--FOXA_LABEL-uHepG2_SORT-NucleosomeEngagement_GROUP-LessEngaged_500bp.bed   (FOXA_HepG2_NoNuc_1000bp.bed)
#   |--1000bp
#      |--FOXA_LABEL-K562_SORT-ClosestDyad.bed                                      (FOXA_K562_NucSort.bed)
#      |--FOXA_LABEL-K562_SORT-ClosestDayd_GROUP-NoOverlap_1000bp.bed               (FOXA_K562_NucSort-OVERLAP_1000bp.bed)
#      |--FOXA_LABEL-K562_SORT-ClosestDyad_GROUP-NFR_1000bp.bed                     (FOXA_K562_NucSort-NFR_1000bp.bed)
#      |--FOXA_SORT-ClosestDyad_STACK-K562-uHepG2_1000bp.bed                        (FOXA_all_1000bp.bed)


### CHANGE ME
WRK=/path/to/2024-Krebs_Science/04_Call_Motifs
WRK=/storage/home/owl5022/scratch/2024-Krebs_Science/04_Call_Motifs
###

# Dependencies
# - bedtools
# - java

set -exo
module load bedtools
module load anaconda
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Inputs and outputs
MOTIF=../data/RefPT-Motif
GENOME=../data/hg38_files/hg38.fa
GINFO=../data/hg38_files/hg38.chrom.sizes.txt
NUCLEOSOME=../data/RefPT-Krebs/BNase-Nucleosomes.bed

[ -d $MOTIF ] || mkdir $MOTIF
[ -d $MOTIF/500bp ] || mkdir $MOTIF/500bp
[ -d $MOTIF/1000bp ] || mkdir $MOTIF/1000bp

# Script shortcuts
SCRIPTMANAGER=../bin/ScriptManager-v0.15.jar

TEMP=temp-5_FoxA_motif
[ -d $TEMP ] || mkdir $TEMP

## =====Intersect and subtract K562 vs HepG2 bound sites=====

HEPG2_FOXA1=temp-3_Filter_and_Sort_by_occupancy/FOXA2_FOXA1-HepG2_M1_100bp_7-Occupancy_BOUND_1bp.bed
HEPG2_FOXA2=temp-3_Filter_and_Sort_by_occupancy/FOXA2_FOXA2-HepG2_M1_100bp_7-Occupancy_BOUND_1bp.bed
K562_FOXA1=temp-3_Filter_and_Sort_by_occupancy/FOXA2_FOXA1-K562_M1_100bp_7-Occupancy_BOUND_1bp.bed
K562_FOXA2=temp-3_Filter_and_Sort_by_occupancy/FOXA2_FOXA2-K562_M1_100bp_7-Occupancy_BOUND_1bp.bed

# Merge and unique HepG2 sites with genomic sort
cat $HEPG2_FOXA1 $HEPG2_FOXA2 \
    | sort -uk4,4 \
    | awk 'BEGIN{FS="\t";OFS="\t"}{print $0,"HepG2"}' \
    | bedtools sort -i \
    > $TEMP/FOXA_HepG2.bed

# Merge and unique K562 sites with genomic sort
cat $K562_FOXA1 $K562_FOXA2 \
    | sort -uk4,4 \
    | awk 'BEGIN{FS="\t";OFS="\t"}{print $0,"K562"}' \
    | bedtools sort -i \
    > $TEMP/FOXA_K562.bed

# QC: Stat each group
wc -l $TEMP/FOXA_K562.bed
wc -l $TEMP/FOXA_HepG2.bed

# Venn Intersect
bedtools intersect -v -a $TEMP/FOXA_HepG2.bed -b $TEMP/FOXA_K562.bed  | awk '{OFS="\t"}{FS="\t"}{print $1,$2,$3,$4,$5,$6,"uHepG2"}' > $TEMP/FOXA_uHepG2.bed
bedtools intersect -u -a $TEMP/FOXA_HepG2.bed -b $TEMP/FOXA_K562.bed  | awk '{OFS="\t"}{FS="\t"}{print $1,$2,$3,$4,$5,$6,"K562-HepG2"}' > $TEMP/FOXA_K562-HepG2.bed
bedtools intersect -v -a $TEMP/FOXA_K562.bed  -b $TEMP/FOXA_HepG2.bed | awk '{OFS="\t"}{FS="\t"}{print $1,$2,$3,$4,$5,$6,"uK562"}' > $TEMP/FOXA_uK562.bed

# QC: Stat each group
wc -l $TEMP/FOXA_uHepG2.bed
wc -l $TEMP/FOXA_K562-HepG2.bed
wc -l $TEMP/FOXA_uK562.bed

# Merge with Venn labels (uK562 + K562-HepG2)
cat $TEMP/FOXA_uK562.bed $TEMP/FOXA_K562-HepG2.bed | bedtools sort -i | uniq > $TEMP/FOXA_LABEL-K562.bed
mv $TEMP/FOXA_uHepG2.bed $TEMP/FOXA_LABEL-uHepG2.bed


## =====Sort FOXA sites by downstream nucleosome engagement=====
# for each of K562 and uHepG2 groups, count 5' end read 2 tags at 150bp downstream

K562_FOXA1_BAM=../data/BAM/K562_FOXA1_BX_rep1_hg38.bam
HEPG2_FOXA1_BAM=../data/BAM/HepG2_FOXA1_BX_rep1_hg38.bam

# Shift 150bp downstream
bedtools shift -i $TEMP/FOXA_LABEL-K562.bed -g $GINFO -p 150 -m -150 > $TEMP/FOXA_LABEL-K562_Shift-d150.bed
bedtools shift -i $TEMP/FOXA_LABEL-uHepG2.bed -g $GINFO -p 150 -m -150 > $TEMP/FOXA_LABEL-uHepG2_Shift-d150.bed

# Expand 100bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $TEMP/FOXA_LABEL-K562_Shift-d150.bed -o $TEMP/FOXA_LABEL-K562_Shift-d150_100bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $TEMP/FOXA_LABEL-uHepG2_Shift-d150.bed -o $TEMP/FOXA_LABEL-uHepG2_Shift-d150_100bp.bed

# Tag Pileup
java -jar $SCRIPTMANAGER read-analysis tag-pileup $TEMP/FOXA_LABEL-K562_Shift-d150_100bp.bed $K562_FOXA1_BAM -2 --cpu 4 -M $TEMP/FOXA_LABEL-K562_Shift-d150_100bp_read2
java -jar $SCRIPTMANAGER read-analysis tag-pileup $TEMP/FOXA_LABEL-uHepG2_Shift-d150_100bp.bed $HEPG2_FOXA1_BAM -2 --cpu 4 -M $TEMP/FOXA_LABEL-uHepG2_Shift-d150_100bp_read2

# Sum antisense scores
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum $TEMP/FOXA_LABEL-K562_Shift-d150_100bp_read2_anti.cdt -o $TEMP/FOXA_LABEL-K562_Shift-d150_100bp_read2_anti.out
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum $TEMP/FOXA_LABEL-uHepG2_Shift-d150_100bp_read2_anti.cdt -o $TEMP/FOXA_LABEL-uHepG2_Shift-d150_100bp_read2_anti.out

# Append scores to initial BED file
tail -n +2 $TEMP/FOXA_LABEL-K562_Shift-d150_100bp_read2_anti.out | cut -f 2 | paste $TEMP/FOXA_LABEL-K562.bed - | sort -k8,8nr > $MOTIF/FOXA_LABEL-K562_SORT-NucleosomeEngagement.bed
tail -n +2 $TEMP/FOXA_LABEL-uHepG2_Shift-d150_100bp_read2_anti.out | cut -f 2 | paste $TEMP/FOXA_LABEL-uHepG2.bed - | sort -k8,8nr > $MOTIF/FOXA_LABEL-uHepG2_SORT-NucleosomeEngagement.bed

# Expand 500bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 500 $MOTIF/FOXA_LABEL-K562_SORT-NucleosomeEngagement.bed -o $MOTIF/500bp/FOXA_LABEL-K562_SORT-NucleosomeEngagement_500bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 500 $MOTIF/FOXA_LABEL-uHepG2_SORT-NucleosomeEngagement.bed -o $MOTIF/500bp/FOXA_LABEL-uHepG2_SORT-NucleosomeEngagement_500bp.bed

# Seperate by nucleosome engagement level (Engaged/LessEngaged, Top/Bottom)
head -n 500 $MOTIF/FOXA_LABEL-K562_SORT-NucleosomeEngagement.bed > $MOTIF/FOXA_LABEL-K562_SORT-NucleosomeEngagement_GROUP-Engaged.bed
tail -n +501 $MOTIF/FOXA_LABEL-K562_SORT-NucleosomeEngagement.bed > $MOTIF/FOXA_LABEL-K562_SORT-NucleosomeEngagement_GROUP-LessEngaged.bed
head -n 900 $MOTIF/FOXA_LABEL-uHepG2_SORT-NucleosomeEngagement.bed > $MOTIF/FOXA_LABEL-uHepG2_SORT-NucleosomeEngagement_GROUP-Engaged.bed
tail -n +901 $MOTIF/FOXA_LABEL-uHepG2_SORT-NucleosomeEngagement.bed > $MOTIF/FOXA_LABEL-uHepG2_SORT-NucleosomeEngagement_GROUP-LessEngaged.bed

# Expand 500bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 500 $MOTIF/FOXA_LABEL-K562_SORT-NucleosomeEngagement_GROUP-Engaged.bed       -o $MOTIF/500bp/FOXA_LABEL-K562_SORT-NucleosomeEngagement_GROUP-Engaged_500bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 500 $MOTIF/FOXA_LABEL-K562_SORT-NucleosomeEngagement_GROUP-LessEngaged.bed   -o $MOTIF/500bp/FOXA_LABEL-K562_SORT-NucleosomeEngagement_GROUP-LessEngaged_500bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 500 $MOTIF/FOXA_LABEL-uHepG2_SORT-NucleosomeEngagement_GROUP-Engaged.bed     -o $MOTIF/500bp/FOXA_LABEL-uHepG2_SORT-NucleosomeEngagement_GROUP-Engaged_500bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 500 $MOTIF/FOXA_LABEL-uHepG2_SORT-NucleosomeEngagement_GROUP-LessEngaged.bed -o $MOTIF/500bp/FOXA_LABEL-uHepG2_SORT-NucleosomeEngagement_GROUP-LessEngaged_500bp.bed


## =====Sort FOXA sites by distance to closest nucleosome and separate into NFR/nucleosome groups=====

# Sort BED files for closest operation
bedtools sort -i $TEMP/FOXA_LABEL-K562.bed > $TEMP/FOXA_LABEL-K562_SORT-Genomic.bed
bedtools sort -i $TEMP/FOXA_LABEL-uHepG2.bed > $TEMP/FOXA_LABEL-uHepG2_SORT-Genomic.bed
bedtools sort -i $NUCLEOSOME > $TEMP/Nucleosomes_SORT-Genomic.bed

# Call closest nucleosome
bedtools closest -d -D a -t all -a $TEMP/FOXA_LABEL-K562_SORT-Genomic.bed -b $TEMP/Nucleosomes_SORT-Genomic.bed > $TEMP/FOXA_LABEL-K562_SORT-DistClosestDyad_Redundant.tsv
bedtools closest -d -D a -t all -a $TEMP/FOXA_LABEL-uHepG2_SORT-Genomic.bed -b $TEMP/Nucleosomes_SORT-Genomic.bed > $TEMP/FOXA_LABEL-uHepG2_SORT-DistClosestDyad_Redundant.tsv

# Random select from ties and sort by distance w/respect to NFIA
shuf $TEMP/FOXA_LABEL-K562_SORT-DistClosestDyad_Redundant.tsv | sort -uk4,4 | sort -nk14,14 > $TEMP/FOXA_LABEL-K562_SORT-ClosestDyad.tsv
shuf $TEMP/FOXA_LABEL-uHepG2_SORT-DistClosestDyad_Redundant.tsv | sort -uk4,4 | sort -nk14,14 > $TEMP/FOXA_LABEL-uHepG2_SORT-ClosestDyad.tsv

# Group by distance to closest nucleosome bounded by -73 and +73 
awk -v DIR="$TEMP" 'BEGIN{OFS="\t";FS="\t"}{
    if ($14 >= -73 && $14 <= 73) {
      print $0 > DIR"/FOXA_LABEL-K562_SORT-ClosestDyad_GROUP-NoOverlap.tsv"
    } else {
      print $0 > DIR"/FOXA_LABEL-K562_SORT-ClosestDyad_GROUP-NFR.tsv"
    }
}' $TEMP/FOXA_LABEL-K562_SORT-ClosestDyad.tsv


# Group by Nucleosome-overlapping or NFR
awk -v DIR="$TEMP" 'BEGIN{OFS="\t";FS="\t"}{
    if ($14 >= -73 && $14 <= 73) {
        print $0 >DIR"/FOXA_LABEL-uHepG2_SORT-ClosestDyad_GROUP-NoOverlap.tsv"
    } else {
        print $0 > DIR"/FOXA_LABEL-uHepG2_SORT-ClosestDyad_GROUP-NFR.tsv"
    }
}' $TEMP/FOXA_LABEL-uHepG2_SORT-ClosestDyad.tsv

# Slice to BED6
cut -f1-6 $TEMP/FOXA_LABEL-K562_SORT-ClosestDyad.tsv                   > $MOTIF/FOXA_LABEL-K562_SORT-ClosestDyad.bed
cut -f1-6 $TEMP/FOXA_LABEL-K562_SORT-ClosestDyad_GROUP-NoOverlap.tsv   > $MOTIF/FOXA_LABEL-K562_SORT-ClosestDyad_GROUP-NoOverlap.bed
cut -f1-6 $TEMP/FOXA_LABEL-K562_SORT-ClosestDyad_GROUP-NFR.tsv         > $MOTIF/FOXA_LABEL-K562_SORT-ClosestDyad_GROUP-NFR.bed
cut -f1-6 $TEMP/FOXA_LABEL-uHepG2_SORT-ClosestDyad.tsv                 > $MOTIF/FOXA_LABEL-uHepG2_SORT-ClosestDyad.bed
cut -f1-6 $TEMP/FOXA_LABEL-uHepG2_SORT-ClosestDyad_GROUP-NoOverlap.tsv > $MOTIF/FOXA_LABEL-uHepG2_SORT-ClosestDyad_GROUP-NoOverlap.bed
cut -f1-6 $TEMP/FOXA_LABEL-uHepG2_SORT-ClosestDyad_GROUP-NFR.tsv       > $MOTIF/FOXA_LABEL-uHepG2_SORT-ClosestDyad_GROUP-NFR.bed

# Concat for stacked BED
cat $MOTIF/FOXA_LABEL-K562_SORT-ClosestDyad.bed $MOTIF/FOXA_LABEL-uHepG2_SORT-ClosestDyad.bed > $MOTIF/FOXA_SORT-ClosestDyad_STACK-K562-uHepG2.bed

# QC: Stat NucSort
wc -l $MOTIF/*SORT-ClosestDyad*.bed

# Expand 1000bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/FOXA_SORT-ClosestDyad_STACK-K562-uHepG2.bed -o $MOTIF/1000bp/FOXA_SORT-ClosestDyad_STACK-K562-uHepG2_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/FOXA_LABEL-K562_SORT-ClosestDyad.bed -o $MOTIF/1000bp/FOXA_LABEL-K562_SORT-ClosestDyad_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/FOXA_LABEL-uHepG2_SORT-ClosestDyad.bed -o $MOTIF/1000bp/FOXA_LABEL-uHepG2_SORT-ClosestDyad_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/FOXA_LABEL-K562_SORT-ClosestDyad_GROUP-NFR.bed -o $MOTIF/1000bp/FOXA_LABEL-K562_SORT-ClosestDyad_GROUP-NFR_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/FOXA_LABEL-K562_SORT-ClosestDyad_GROUP-NoOverlap.bed -o $MOTIF/1000bp/FOXA_LABEL-K562_SORT-ClosestDyad_GROUP-NoOverlap_1000bp.bed

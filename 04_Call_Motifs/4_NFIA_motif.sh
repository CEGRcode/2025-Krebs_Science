#!/bin/bash

# Determine motif orientation for the palindromic NFIA motif and create RefPT files by occupancy and dist to closest dyad

# data/RefPT-Motif
#   |--NFIA_SORT-DistClosestDyad.bed                                        (see 1000bp)
#   |--NFIA_SORT-DistClosestDyad_GROUP-Downstream.bed                       (see 1000bp)
#   |--NFIA_SORT-DistClosestDyad_GROUP-Overlap.bed                          (see 1000bp)
#   |--NFIA_SORT-DistClosestDyad_GROUP-Upstream.bed                         (see 1000bp)
#   |--NFIA_SORT-Occupancy.bed                                              (see 500bp)
#   |--NFIA-u95_SORT-Occupancy.bed
#   |--NFIA-d95_SORT-Occupancy.bed
#   |--NFIA_REORIENT-Random_SORT-Occupancy.bed                              (NFIA_randomly_orientated.bed)
#   |--NFIA-u95_REORIENT-Random_SORT-Occupancy.bed
#   |--NFIA-d95_REORIENT-Random_SORT-Occupancy.bed
#   |--150bp
#      |--NFIA-u95_SORT-Occupancy_150bp.bed
#      |--NFIA-d95_SORT-Occupancy_150bp.bed
#      |--NFIA-u95_REORIENT-Random_SORT-Occupancy_150bp.bed
#      |--NFIA-d95_REORIENT-Random_SORT-Occupancy_150bp.bed
#   |--250bp
#      |--NFIA-d250bp_SORT-Occupancy_250bp.bed                              (NFIA_downNuc_down125_250bp.bed)
#   |--500bp
#      |--NFIA_SORT-Occupancy_500bp.bed                                     (NFIA_downNuc_500bp.bed)
#   |--1000bp
#      |--NFIA_SORT-DistClosestDyad_1000bp.bed                              (NFIA_engageNucdown_NucSort_1000bp.bed)
#      |--NFIA_SORT-DistClosestDyad_GROUP-Downstream_1000bp.bed             (NFIA_NucSort-DOWNSTREAM_1000bp)
#      |--NFIA_SORT-DistClosestDyad_GROUP-Overlap_1000bp.bed                (NFIA_NucSort-OVERLAP_1000bp)
#      |--NFIA_SORT-DistClosestDyad_GROUP-Upstream_1000bp.bed               (NFIA_NucSort-UPSTREAM_1000bp)
#      |--NFIA_REORIENT-DownstreamDyad_SORT-Occupancy_GROUP-Q4_1000bp.bed   (NFIA_downNuc_4_1000bp.bed)

### CHANGE ME
WRK=/path/to/2024-Krebs_Science/04_Call_Motifs
WRK=/scratch/owl5022/2024-Krebs_Science/04_Call_Motifs
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
BAMFILE=../data/BAM/K562_NFIA_BX_rep1_hg38.bam
NUCLEOSOME=../data/RefPT-Krebs/BNase-Nucleosomes.bed
BOUND_NFIA=temp-3_Filter_and_Sort_by_occupancy/NFIA_NFIA-K562_M1_100bp_7-Occupancy_BOUND_1bp.bed

[ -d $MOTIF ] || mkdir $MOTIF
[ -d $MOTIF/150bp ] || mkdir $MOTIF/150bp
[ -d $MOTIF/250bp ] || mkdir $MOTIF/250bp
[ -d $MOTIF/500bp ] || mkdir $MOTIF/500bp
[ -d $MOTIF/1000bp ] || mkdir $MOTIF/1000bp

# Script shortcuts
SCRIPTMANAGER=../bin/ScriptManager-v0.15.jar
SHUFFLE=../bin/shuffle_script.py

TEMP=temp-4_NFIA
[ -d $TEMP ] || mkdir $TEMP

## =====Set orientation for palindromic motif======
# by occupancy of exo stop sites (upstream vs downstream)

# Create RefPT shifting up and down 95bp
bedtools shift -i $BOUND_NFIA -g $GINFO -p 95 -m -95 | cut -f1-6 > $TEMP/NFIA_down95bp.bed
bedtools shift -i $BOUND_NFIA -g $GINFO -p -95 -m 95 | cut -f1-6 > $TEMP/NFIA_up95bp.bed

# Extract tag counts for up and down shift groups...
for GROUP in "down95bp" "up95bp";
do
    BASE=$TEMP/NFIA_$GROUP
    # Expand 150bp
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 150 $BASE.bed -o ${BASE}_150bp.bed
    # Tag Pileup
    java -jar $SCRIPTMANAGER read-analysis tag-pileup ${BASE}_150bp.bed $BAMFILE -n 100 -1 --combined --cpu 4 -M ${BASE}_150bp_5read1_MIN100
    # Aggregate
    java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum ${BASE}_150bp_5read1_MIN100_combined.cdt -o ${BASE}_150bp_5read1_MIN100_combined.out
done

# Merge up and down tag count scores
paste $TEMP/NFIA_up95bp_150bp_5read1_MIN100_combined.out $TEMP/NFIA_down95bp_150bp_5read1_MIN100_combined.out \
    | tail -n +2 | cut -f 2,4 > $TEMP/NFIA_UpDownTagCounts.out

# Append to original BED file ($BOUND_NFIA)
paste $BOUND_NFIA $TEMP/NFIA_UpDownTagCounts.out > $TEMP/NFIA_7-Occupancy_8-up_9-down.bed

# Compare up/down scores and split by up > down or not
awk '{if ($8 >= $9) print $0,"engageNucDown" }' $TEMP/NFIA_7-Occupancy_8-up_9-down.bed > $TEMP/NFIA_EngagedDownstream.bed
awk '{if ($8 <  $9) print $0,"engageNucUp"}' $TEMP/NFIA_7-Occupancy_8-up_9-down.bed > $TEMP/NFIA_EngagedUpstream.bed

# QC: Stat line counts
wc -l $TEMP/NFIA_EngagedDownstream.bed
# 4093 temp-4_NFIA/NFIA_EngagedDownstream.bed
wc -l $TEMP/NFIA_EngagedUpstream.bed
# 3873 temp-4_NFIA/NFIA_EngagedUpstream.bed

# Flip strands of upstream (puts nucleosome downstream)
awk 'BEGIN{OFS="\t";FS="\t"}{
    if ($6 == "-") {
        $6="+";
    } else {
        $6="-";
    }
    print $1,$2,$3,$4,$5,$6,$7,$9,$8,$10
}' $TEMP/NFIA_EngagedUpstream.bed > $TEMP/NFIA_EngagedUpstream_FlipStrand.bed

# Shift 1bp downstream
bedtools shift -p 1 -m -1 -g $GINFO -i $TEMP/NFIA_EngagedUpstream_FlipStrand.bed > $TEMP/NFIA_EngagedUpstream_FlipStrand_shift1bp.bed

# Merge reoriented up with down slices
cat $TEMP/NFIA_EngagedUpstream_FlipStrand_shift1bp.bed $TEMP/NFIA_EngagedDownstream.bed > $TEMP/NFIA-Reoriented_7-Occupancy.bed


## =====Create occupancy RefPTs======

# Re-sort by original occupancy (col 7) and slice to BED6
sort -k7,7nr $TEMP/NFIA-Reoriented_7-Occupancy.bed | cut -f1-6 > $MOTIF/NFIA_SORT-Occupancy.bed

# Expand 500bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 500 $MOTIF/NFIA_SORT-Occupancy.bed -o $MOTIF/500bp/NFIA_SORT-Occupancy_500bp.bed

# Slop upstream window
bedtools slop -l -250 -r 0 -s -i $MOTIF/500bp/NFIA_SORT-Occupancy_500bp.bed -g $GINFO > $MOTIF/250bp/NFIA-d250bp_SORT-Occupancy_250bp.bed


## =====Shift RefPTs up/down for violins=====

# Create RefPT shifting up and down 95bp
bedtools shift -i $MOTIF/NFIA_SORT-Occupancy.bed -g $GINFO -p 95 -m -95 > $MOTIF/NFIA-u95_SORT-Occupancy.bed
bedtools shift -i $MOTIF/NFIA_SORT-Occupancy.bed -g $GINFO -p -95 -m 95 > $MOTIF/NFIA-d95_SORT-Occupancy.bed

# Expand 150bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 150 $MOTIF/NFIA-u95_SORT-Occupancy.bed -o $MOTIF/150bp/NFIA-u95_SORT-Occupancy_150bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 150 $MOTIF/NFIA-d95_SORT-Occupancy.bed -o $MOTIF/150bp/NFIA-d95_SORT-Occupancy_150bp.bed


## =====Slice bottom quartile=====

# # Take bottom quartile (Q4) of NFIA_downNuc
# tail -n 1991 $MOTIF/NFIA_REORIENT-DownstreamDyad_SORT-Occupancy.bed > $MOTIF/NFIA_REORIENT-DownstreamDyad_SORT-Occupancy_GROUP-Q4.bed

# # Expand 1000bp
# java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/NFIA_downNuc_GROUP-Q4.bed -o $MOTIF/1000bp/NFIA_downNuc_GROUP-Q4_1000bp.bed


# =====Sort by closest nucleosome and group by distance=====

# Sort BED files for closest operation
bedtools sort -i $MOTIF/NFIA_SORT-Occupancy.bed > $TEMP/NFIA_SORT-Genomic.bed
bedtools sort -i $NUCLEOSOME > $TEMP/Nucleosomes_SORT-Genomic.bed

# Call closest nucleosome
bedtools closest -d -D a -t all -a $TEMP/NFIA_SORT-Genomic.bed -b $TEMP/Nucleosomes_SORT-Genomic.bed > $TEMP/NFIA_SORT-DistClosestDyad_Redundant.tsv

# Random select from ties and sort by distance w/respect to NFIA
shuf $TEMP/NFIA_SORT-DistClosestDyad_Redundant.tsv | sort -uk4,4 | sort -nk13,13 > $TEMP/NFIA_SORT-DistClosestDyad.tsv

# Group by distance to closest nucleosome bounded by -73 and +73 (3 groups)
awk -v DIR="$MOTIF" 'BEGIN{OFS="\t";FS="\t"}{
   if ($13 < -73 ) {
      print $0 > DIR"/NFIA_SORT-DistClosestDyad_GROUP-Upstream.bed"
    } else if ($13 > 73 ) {
      print $0 > DIR"/NFIA_SORT-DistClosestDyad_GROUP-Downstream.bed"
    } else {
      print $0 > DIR"/NFIA_SORT-DistClosestDyad_GROUP-Overlap.bed"
    }
}' $TEMP/NFIA_SORT-DistClosestDyad.tsv

# Slice tsv into BED6
cut -f1-6 $TEMP/NFIA_SORT-DistClosestDyad.tsv > $MOTIF/NFIA_SORT-DistClosestDyad.bed

# Expand 1000bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/NFIA_SORT-DistClosestDyad.bed -o $MOTIF/1000bp/NFIA_SORT-DistClosestDyad_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/NFIA_SORT-DistClosestDyad_GROUP-Upstream.bed -o $MOTIF/1000bp/NFIA_SORT-DistClosestDyad_GROUP-Upstream_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/NFIA_SORT-DistClosestDyad_GROUP-Downstream.bed -o $MOTIF/1000bp/NFIA_SORT-DistClosestDyad_GROUP-Downstream_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/NFIA_SORT-DistClosestDyad_GROUP-Overlap.bed -o $MOTIF/1000bp/NFIA_SORT-DistClosestDyad_GROUP-Overlap_1000bp.bed

# QC: Stat line counts
wc -l $MOTIF/NFIA_SORT-DistClosestDyad_GROUP-*.bed
# 1563 ../data/RefPT-Motif/NFIA_SORT-DistClosestDyad_GROUP-Downstream.bed
# 5248 ../data/RefPT-Motif/NFIA_SORT-DistClosestDyad_GROUP-Overlap.bed
# 1155 ../data/RefPT-Motif/NFIA_SORT-DistClosestDyad_GROUP-Upstream.bed

## =====Randomly reorient nucleosomes=====
# by shuffling tag pileup matrix values and re-calling orientation by exo occupancy

# Concatenate up95/down95 pileup CDTs to create an extended matrix
cut -f 3-152 $TEMP/NFIA_up95bp_150bp_5read1_MIN100_combined.cdt \
    | paste $TEMP/NFIA_down95bp_150bp_5read1_MIN100_combined.cdt - \
    > $TEMP/NFIA_CONCAT-updown.cdt

# Shuffle the CDT values from up and down on a per-motif basis
python $SHUFFLE $TEMP/NFIA_CONCAT-updown.cdt $TEMP/NFIA_shuffled.cdt

# Sum the rows for each
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum <(cut -f 1-152 $TEMP/NFIA_shuffled.cdt) -o $TEMP/NFIA_shuffled-down.out
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum <(cut -f 1-2,153-302 $TEMP/NFIA_shuffled.cdt) -o $TEMP/NFIA_shuffled-up.out

# Merge up and down tag count scores
paste $TEMP/NFIA_shuffled-up.out $TEMP/NFIA_shuffled-down.out \
    | tail -n +2 | cut -f2,4 > $TEMP/NFIA_shuffled-UpDownTagCounts.out

# Append to original BED file ($BOUND_NFIA)
paste $BOUND_NFIA $TEMP/NFIA_shuffled-UpDownTagCounts.out > $TEMP/NFIA_shuffled_7-Occupancy_8-upshuffle_9-downshuffle.bed

# Compare up/down scores and split by up > down or not
awk '{if ($8 >= $9) print $0,"engageNucDown" }' $TEMP/NFIA_shuffled_7-Occupancy_8-upshuffle_9-downshuffle.bed > $TEMP/NFIA_shuffled_EngagedDownstream.bed
awk '{if ($8 <  $9) print $0,"engageNucUp"}' $TEMP/NFIA_shuffled_7-Occupancy_8-upshuffle_9-downshuffle.bed > $TEMP/NFIA_shuffled_EngagedUpstream.bed

# QC: Stat line counts
wc -l $TEMP/NFIA_shuffled_EngagedUpstream.bed
wc -l $TEMP/NFIA_shuffled_EngagedDownstream.bed

# Flip strands of upstream (puts nucleosome downstream)
awk 'BEGIN{OFS="\t";FS="\t"}{
    if ($6 == "-") {
        $6="+";
    } else {
        $6="-";
    }
    print $1,$2,$3,$4,$5,$6,$7,$9,$8,$10
}' $TEMP/NFIA_shuffled_EngagedUpstream.bed > $TEMP/NFIA_shuffled_EngagedUpstream_FlipStrand.bed

# Shift 1bp downstream
bedtools shift -p 1 -m -1 -g $GINFO -i $TEMP/NFIA_shuffled_EngagedUpstream_FlipStrand.bed > $TEMP/NFIA_shuffled_EngagedUpstream_FlipStrand_shift1bp.bed

# Merge reoriented up with down slices
cat $TEMP/NFIA_shuffled_EngagedUpstream_FlipStrand_shift1bp.bed $TEMP/NFIA_shuffled_EngagedDownstream.bed | sort -k7,7nr > $MOTIF/NFIA_REORIENT-Random_SORT-Occupancy.bed

## =====Shift randomly reoriented RefPTs up/down for violins=====

# Create RefPT shifting up and down 95bp
bedtools shift -i $MOTIF/NFIA_REORIENT-Random_SORT-Occupancy.bed -g $GINFO -p 95 -m -95 | cut -f1-6 > $MOTIF/NFIA-u95_REORIENT-Random_SORT-Occupancy.bed
bedtools shift -i $MOTIF/NFIA_REORIENT-Random_SORT-Occupancy.bed -g $GINFO -p -95 -m 95 | cut -f1-6 > $MOTIF/NFIA-d95_REORIENT-Random_SORT-Occupancy.bed

# Expand 150bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 150 $MOTIF/NFIA-u95_REORIENT-Random_SORT-Occupancy.bed -o $MOTIF/150bp/NFIA-u95_REORIENT-Random_SORT-Occupancy_150bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 150 $MOTIF/NFIA-d95_REORIENT-Random_SORT-Occupancy.bed -o $MOTIF/150bp/NFIA-d95_REORIENT-Random_SORT-Occupancy_150bp.bed

# Violin plot commands should be moved to Z_Figures
# for file in $TEMP/NFIA_shuffled_downNuc_down95.bed $TEMP/NFIA_shuffled_downNuc_up95.bed ; do
#   filename=$(basename $file ".bed")
#   java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 150 $file -o ${filename}_150bp.bed
#   java -jar $SCRIPTMANAGER read-analysis tag-pileup ${filename}_150bp.bed $NFIA1 -n 100 -1 --combined --cpu 4 -M SCORES/NFIA_${filename}_150bp_read1_MIN100
#   java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum SCORES/NFIA_${filename}_150bp_read1_MIN100_combined.cdt -o SCORES/
#   rm ${filename}_150bp.bed
# done

# cut -f 2 SCORES/NFIA_NFIA_shuffled_downNuc_down95_150bp_read1_MIN100_combined_SCORES.out \
#     | paste SCORES/NFIA_NFIA_shuffled_downNuc_up95_150bp_read1_MIN100_combined_SCORES.out - \
#     | tail -n +2 | cut -f 2-3 \
#     | paste NFIA_shuffled_downNuc.bed - \
#     > NFIA_shuffled_downNuc_SCORES.bed
# rm NFIA_shuffled_downNuc_down95.bed NFIA_shuffled_downNuc_up95.bed
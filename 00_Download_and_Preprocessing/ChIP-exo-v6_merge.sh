#!/bin/bash

# Script to hardcode the merging/renaming of PEGR BAM & MEME files into a standard file naming system

### CHANGE ME
WRK=/storage/group/bfp2/default/hxc585_HainingChen/2025_Chen_TF-Nuc/00_Download_and_Preprocessing
WRK=/storage/home/owl5022/scratch/2024-Krebs_Science/00_Download_and_Preprocessing
###

# Dependencies
# - java
# - samtools

set -exo
module load samtools

# Inputs and outputs
BAMDIR=$WRK/../data/BAM

# Script shortcuts
PICARD=$WRK/../bin/picard.jar

[ -d $BAMDIR ] || mkdir $BAMDIR

cd ../data/sample-BAM

# Merged IgG (K562-BX) (HepG2-BX)
java -jar $PICARD MergeSamFiles -O K562_IgG_BX_merge_hg38.bam \
        -I 33925_IgG_i5006_K562_-_IMDM_-_BX.bam \
        -I 33963_IgG_i5006_K562_-_IMDM_-_BX.bam \
        -I 34031_IgG_i5006_K562_-_IMDM_-_BX.bam \
        -I 34055_IgG_i5006_K562_-_IMDM_-_BX.bam \
        -I 34175_IgG_i5006_K562_-_IMDM_-_BX.bam \
        -I 34469_IgG_i5006_K562_-_IMDM_-_BX.bam \
        -I 36713_IgG_i5006_K562_-_IMDM_-_BX.bam \
        -I 36900_IgG_i5006_K562_-_IMDM_-_BX.bam \
        -I 37449_IgG_i5006_K562_-_IMDM_-_BX.bam \

        # -I 34668_IgG_i5006_K562_-_IMDM_-_BX.bam \
        # -I 38288_IgG_i5006_K562_-_IMDM_-_BX.bam
# HAINING QUESTION: What are these samples doing?


java -jar $PICARD MergeSamFiles -O HepG2_IgG_BX_merge_hg38.bam \
        -I 38481_IgG_i5006_HepG2_-_-_-_BX.bam \
        -I 38482_IgG_i5006_HepG2_-_-_-_BX.bam \
        -I 38483_IgG_i5006_HepG2_-_-_-_BX.bam

# Use Picard to merge resequenced technical replicates, otherwise rename BAM using cp

# Input
cp 28326_Input_-_K562_-_-_50Unuclease10min-0cycSonic-newQuenchingbuffer_BI.bam K562_-_BI_rep1_hg38.bam

# CTCF
cp 33924_CTCF_07-729_K562_-_IMDM_-_BX.bam K562_CTCF_BX_rep1_hg38.bam
cp 34174_CTCF_07-729_K562_-_IMDM_-_BX.bam K562_CTCF_BX_rep2_hg38.bam

# non-crosslinked CTCF 
java -jar $PICARD MergeSamFiles -O K562_CTCF_nonXLBX_rep1_hg38.bam \
        -I 29009_CTCF_07-729_K562_-_-_2stepNonXLlysis-50Unuclease-10min-0cycSonic-quenchONLY_BX.bed \
        -I 29010_CTCF_07-729_K562_-_-_2stepNonXLlysis-50Unuclease-10min-0cycSonic-quenchONLY_BX.bed

# RAD21
java -jar $PICARD MergeSamFiles -O K562_RAD21_BX_rep1_hg38.bam \
        -I 32070_RAD21_HPA020044_K562_-_-_-_BX.bam \
        -I 32119_RAD21_HPA020044_K562_-_-_-_BX.bam

cp 34621_RAD21_HPA020044_K562_-_-_-_BX.bam K562_RAD21_BX_rep2_hg38.bam
# SMC3
cp 34622_Smc3_ab9263_K562_-_-_-_BX.bam K562_SMC3_BX_rep1_hg38.bam
cp 38476_Smc3_ab9263_K562_-_IMDM_-_BX.bam K562_SMC3_BX_rep2_hg38.bam


# NFIA
java -jar $PICARD MergeSamFiles -O K562_NFIA_BX_rep1_hg38.bam \
        -I 32067_NFIA_HPA008884_K562_-_-_-_BX.bed \
        -I 32116_NFIA_HPA008884_K562_-_-_-_BX.bed


cp 34598_NFIA_HPA008884_K562_-_-_-BX.bed K562_NFIA_BX_rep2_hg38.bam

# FOXA1 (HepG2)
cp 32354_FoxA1_ab23738_HepG2_-_-_-_BX.bam HepG2_FOXA1_BX_rep1_hg38.bam
cp 33219_FOXA1_ab23738_HepG2_-_-_-_BX.bam HepG2_FOXA1_BX_rep2_hg38.bam

# FOXA1 (K562)
cp 38477_FoxA1_ab23738_K562_-_IMDM_-_BX.bam K562_FOXA1_BX_rep1_hg38.bam
cp 38478_FoxA1_ab23738_K562_-_IMDM_-_BX.bam K562_FOXA1_BX_rep2_hg38.bam

# FOXA2 (HepG2)
cp 32467_FOXA2_ab256493_HepG2_-_-_-_BX.bam HepG2_FOXA2_BX_rep2_hg38.bam
cp 38480_FOXA2_ab256493_HepG2_-_-_-_BX.bam HepG2_FOXA2_BX_rep1_hg38.bam

# HNF4A (HepG2)
cp 34158_HNF4A_HPA004712_HepG2_-_-_-_BX.bam HepG2_HNF4A_BX_rep1_hg38.bam

# FOXA2 (K562)
cp 32659_FOXA2_ab256493_K562_-_-_-_BX.bam K562_FOXA2_BX_rep1_hg38.bam
cp 38479_FOXA2_ab256493_K562_-_IMDM_-_BX.bam K562_FOXA2_BX_rep2_hg38.bam

# ZKSCAN1 K562
cp  34048_ZKSCAN1_HPA006672_K562_-_IMDM_-_BX.bam K562_ZKSCAN1_BX_rep1_hg38.bam

cd $WRK/../data
mv sample-BAM/K562_*.bam $BAMDIR
mv sample-BAM/HepG2_*.bam $BAMDIR

# Index set of BAM files
for FILE in $BAMDIR/*.bam;
do
  [ -f $FILE.bai ] || samtools index $FILE
done

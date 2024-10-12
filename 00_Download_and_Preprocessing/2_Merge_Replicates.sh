#!/bin/bash

# Script to hardcode the merging/renaming of PEGR BAM & MEME files into a standard file naming system

### CHANGE ME
WRK=/Path/to/Title/00_Download_and_Preprocessing
###

module load samtools

# Script shortcuts
PICARD=$WRK/../bin/picard.jar

cd $WRK/../data/sample-BAM
[ -d ../BAM ] || mkdir ../BAM


# Merged IgG (K562-BX) (HepG2-BX)
java -jar $PICARD MergeSamFiles -I 33925_IgG_i5006_K562_-_IMDM_-_BX.bam \
                                -I 34031_IgG_i5006_K562_-_IMDM_-_BX.bam \
                                -I 34055_IgG_i5006_K562_-_IMDM_-_BX.bam \
                                -I 34175_IgG_i5006_K562_-_IMDM_-_BX.bam \
                                -I 33963_IgG_i5006_K562_-_IMDM_-_BX.bam \
                                -I 36713_IgG_i5006_K562_-_IMDM_-_BX.bam \
                                -I 37449_IgG_i5006_K562_-_IMDM_-_BX.bam \
                                -I 34668_IgG_i5006_K562_-_-_-_BX.bam \
                                -I 34469_IgG_i5006_K562_-_IMDM_-_BX.bam \
                                -I 36900_IgG_i5006_K562_-_IMDM_-_BX.bam \
                                -O K562_IgG_BX_merge_hg38.bam

java -jar $PICARD MergeSamFiles -I 38481_IgG_i5006_HepG2_-_-_-_BX.bam \
                                -I 38482_IgG_i5006_HepG2_-_-_-_BX.bam \
                                -I 38483_IgG_i5006_HepG2_-_-_-_BX.bam \
                                -O HepG2_IgG_BX_merge_hg38.bam

# Use Picard to merge resequenced technical replicates, otherwise rename BAM using cp

# Input
cp 28326_Input_-_K562_-_-_50Unuclease10min-0cycSonic-newQuenchingbuffer_BI.bam K562_-_BI_rep1_hg38.bam

# CTCF
cp 33924_CTCF_07-729_K562_-_IMDM_-_BX.bam K562_CTCF_BX_rep1_hg38.bam
cp 34174_CTCF_07-729_K562_-_IMDM_-_BX.bam K562_CTCF_BX_rep2_hg38.bam
# RAD21
java -jar $PICARD MergeSamFiles -I 32070_RAD21_HPA020044_K562_-_-_-_BX.bam \
                                -I 32119_RAD21_HPA020044_K562_-_-_-_BX.bam \
                                -O K562_RAD21_BX_rep1_hg38.bam
cp 34621_RAD21_HPA020044_K562_-_-_-_BX.bam K562_RAD21_BX_rep2_hg19.bam
# SMC3
cp 34622_Smc3_ab9263_K562_-_-_-_BX.bam K562_SMC3_BX_rep1_hg38.bam
cp 38476_Smc3_ab9263_K562_-_IMDM_-_BX.bam K562_SMC3_BX_rep2_hg38.bam


# NFIA
java -jar $PICARD MergeSamFiles -I 32067_NFIA_HPA008884_K562_-_-_-_BX_FilteredBAM.bed \
                                -I 32116_NFIA_HPA008884_K562_-_-_-_BX_FilteredBAM.bed \
                                -O K562_NFIA_BX_rep1_hg38.bam
cp 34598_NFIA_HPA008884_K562_-_-_-BX_FilteredBAM.bed K562_NFIA_BX_rep2_hg38.bam

# FoxA1 HepG2, K562
cp 32354_FoxA1_ab23738_HepG2_-_-_-_BX.bam HepG2_FOXA1_BX_rep1_hg38.bam
cp 33219_FOXA1_ab23738_HepG2_-_-_-_BX.bam HepG2_FOXA1_BX_rep2_hg38.bam
cp 38477_FoxA1_ab23738_K562_-_IMDM_-_BX.bam K562_FOXA1_BX_rep1_hg38.bam
cp 38478_FoxA1_ab23738_K562_-_IMDM_-_BX.bam K562_FOXA1_BX_rep2_hg38.bam

# FoxA2 HepG2 K562
cp  32467_FOXA2_ab256493_HepG2_-_-_-_BX.bam  HepG2_FOXA2_BX_rep2_hg38.bam
cp  38480_FOXA2_ab256493_HepG2_-_-_-_BX.bam  HepG2_FOXA2_BX_rep1_hg38.bam
cp  32659_FOXA2_ab256493_K562_-_-_-_BX.bam K562_FOXA2_BX_rep1_hg38.bam
cp  38479_FOXA2_ab256493_K562_-_IMDM_-_BX.bam K562_FOXA2_BX_rep2_hg38.bam

cd $WRK/../data
mv sample-BAM/K562_*.bam BAM/
mv sample-BAM/HepG2_*.bam BAM/

# Index set of BAM files
for FILE in BAM/*.bam;
do
  [ -f $FILE.bai ] || samtools index $FILE
done

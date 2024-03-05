#!/bin/bash

# Script to hardcode the merging/renaming of PEGR BAM & MEME files into a standard file naming system

### CHANGE ME
WRK=/path/to/2023-Chen_Benzonase-ChIP-exo/00_Download_and_Preprocessing
###

module load samtools

# Script shortcuts
PICARD=$WRK/../bin/picard.jar

cd $WRK/../data/sample-BAM
[ -d ../BAM ] || mkdir ../BAM

# Merged IgG (K562-BX)
java -jar $PICARD MergeSamFiles -I 29104_IgG_i5006_K562_-_-_-_BX.bam \
                                -I 29246_IgG_i5006_K562_-_-_-_BX.bam \
                                -I 32091_IgG_i5006_K562_-_RPMI_-_BX.bam \
                                -I 32110_IgG_i5006_K562_-_-_-_BX.bam \
                                -I 32636_IgG_i5006_K562_-_-_Lysis-50Unuclease-10cycSonic-washesWithoutPVA-Exo-3splintLigation-shortExo-HCelution-16h-NoAmpure-HCday3-splintOligoswith3primeddC-40ulPCR-standardBuffer-62degreeAnnealing_BX.bam \
                                -I 32704_IgG_i5006_K562_-_IMDM_-_BX.bam \
                                -O K562_IgG_BX_merge_hg19.bam

# Use Picard to merge resequenced technical replicates, otherwise rename BAM using cp

# Input
# cp 28326_Input_-_K562_-_-_50Unuclease10min-0cycSonic-newQuenchingbuffer_BI.bam K562_-_BI_rep1_hg19.bam
# CTCF
cp 30478_CTCF_07-729_K562_-_-_Lysis-50Unuclease-4cycSonic-NoPapainDigestion-bothAmpureSteps_BX.bam K562_CTCF_BX_rep1_hg19.bam
cp 33924_CTCF_07-729_K562_-_IMDM_-_BX.bam K562_CTCF_BX_rep2_hg19.bam
# EP300
java -jar $PICARD MergeSamFiles -I 32064_EP300_HPA003128_K562_-_-_-_BX.bam \
                                -I 32113_EP300_HPA003128_K562_-_-_-_BX.bam \
                                -O K562_EP300_BX_rep1_hg19.bam
cp 32705_p300_A300-358A_K562_-_IMDM_-_BX.bam K562_EP300_BX_rep2_hg19.bam
#GABPA
java -jar $PICARD MergeSamFiles -I 32065_GABPA_HPA003258_K562_-_-_-_BX.bam \
                                -I 32114_GABPA_HPA003258_K562_-_-_-_BX.bam \
                                -O K562_GABPA_BX_rep1_hg19.bam
cp 34607_GABPA_HPA003258_K562_-_-_-_BX.bam K562_GABPA_BX_rep2_hg19.bam
# GTF2B
java -jar $PICARD MergeSamFiles -I 32066_GTF2B_HPA061626_K562_-_-_-_BX.bam \
                                -I 32115_GTF2B_HPA061626_K562_-_-_-_BX.bam \
                                -O K562_GTF2B_BX_rep1_hg19.bam
cp 34617_GTF2B_HPA061626_K562_-_-_-_BX.bam K562_GTF2B_BX_rep2_hg19.bam



# cp K562_H3K27ac_BX_rep1_hg19	28810
# cp K562_H3K4me3_BX_rep1_hg19	28809


# NFIA
java -jar $PICARD MergeSamFiles -I 32067_NFIA_HPA008884_K562_-_-_-_BX.bam \
                                -I 32116_NFIA_HPA008884_K562_-_-_-_BX.bam \
                                -O K562_NFIA_BX_rep1_hg19.bam
cp 34598_NFIA_HPA008884_K562_-_-_-_BX.bam K562_NFIA_BX_rep2_hg19.bam
# Pol2
java -jar $PICARD MergeSamFiles -I 31645_PolII_ab76123_K562_-_RPMI_-_BX.bam \
                                -I 32389_PolII_ab76123_K562_-_-_-_BX.bam \
                                -O K562_Pol2_BX_rep1_hg19.bam
cp 31646_PolII_ab76123_K562_-_RPMI_-_BX.bam K562_Pol2_BX_rep2_hg19.bam
# RAD21
java -jar $PICARD MergeSamFiles -I 32070_RAD21_HPA020044_K562_-_-_-_BX.bam \
                                -I 32119_RAD21_HPA020044_K562_-_-_-_BX.bam \
                                -O K562_RAD21_BX_rep1_hg19.bam
cp 34621_RAD21_HPA020044_K562_-_-_-_BX.bam K562_RAD21_BX_rep2_hg19.bam
# RBBP5
cp 34601_RBBP5_A300-109A_K562_-_-_-_BX.bam K562_RBBP5_BX_rep1_hg19.bam
cp 34600_RBBP5_A300-109A_K562_-_-_-_BX.bam K562_RBBP5_BX_rep2_hg19.bam
# SMC3
cp 34622_Smc3_ab9263_K562_-_-_-_BX.bam K562_SMC3_BX_rep1_hg19.bam
cp 29518_Smc3_ab9263_K562_-_-_benzonase_XO.bam K562_SMC3_XO_rep1_hg19.bam
# SP1
java -jar $PICARD MergeSamFiles -I 32071_Sp1_HPA001853_K562_-_-_-_BX.bam \
                                -I 32120_Sp1_HPA001853_K562_-_-_-_BX.bam \
                                -O K562_SP1_BX_rep1_hg19.bam
cp 34604_Sp1_HPA001853_K562_-_-_-_BX.bam K562_SP1_BX_rep2_hg19.bam
# TBP
java -jar $PICARD MergeSamFiles -I 32662_TBP_hTBPserum_K562_-_-_-_BX.bam \
                                -I 32755_TBP_hTBPserum_K562_-_-_-_BX.bam \
                                -O K562_TBP_BX_rep1_hg19.bam
cp 32661_TBP_hTBPpurified_K562_-_-_-_BX.bam K562_TBP_BX_rep2_hg19.bam
# WDR5
cp 33918_WDR5_HPA047182_K562_-_IMDM_-_BX.bam K562_WDR5_BX_rep1_hg19.bam
cp 34606_WDR5_HPA047182_K562_-_-_-_BX.bam K562_WDR5_BX_rep2_hg19.bam

# Move renamed files into data/BAM
cd $WRK/../data
mv sample-BAM/K562_*.bam BAM/

# Index set of BAM files
for FILE in BAM/*.bam;
do
  [ -f $FILE.bai ] || samtools index $FILE
done

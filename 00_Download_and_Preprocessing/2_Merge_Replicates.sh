#!/bin/bash

# Script to hardcode the merging/renaming of PEGR BAM & MEME files into a standard file naming system

### CHANGE ME
WRK=/path/to/2023-Chen_Benzonase-ChIP-exo/00_Download_and_Preprocessing
###

# Script shortcuts
PICARD=$WRK/../bin/picard.jar

cd $WRK/../data/sample-BAM
[ -d ../BAM ] || mkdir ../BAM

# Resequencing replicates merged
java -jar $PICARD MergeSamFiles -I 31645_PolII_ab76123_K562_-_RPMI_-_BX_FilteredBAM.bam \
                                -I 32389_PolII_ab76123_K562_-_-_-_BX_FilteredBAM.bam \
                                -O K562_Pol2_BX_rep1_hg19.bam

java -jar $PICARD MergeSamFiles -I 32067_NFIA_HPA008884_K562_-_-_-_BX_FilteredBAM.bam \
                                -I 32116_NFIA_HPA008884_K562_-_-_-_BX_FilteredBAM.bam \
                                -O K562_NFIA_BX_rep1_hg19.bam

java -jar $PICARD MergeSamFiles -I 32071_Sp1_HPA001853_K562_-_-_-_BX_FilteredBAM.bam \
                                -I 32120_Sp1_HPA001853_K562_-_-_-_BX_FilteredBAM.bam \
                                -O K562_SP1_BX_rep1_hg19.bam

java -jar $PICARD MergeSamFiles -I 32064_EP300_HPA003128_K562_-_-_-_BX_FilteredBAM.bam \
                                -I 32113_EP300_HPA003128_K562_-_-_-_BX_FilteredBAM.bam \
                                -O K562_EP300_BX_rep1_hg19.bam

java -jar $PICARD MergeSamFiles -I 32114_GABPA_HPA003258_K562_-_-_-_BX_FilteredBAM.bam \
                                -I 32065_GABPA_HPA003258_K562_-_-_-_BX_FilteredBAM.bam \
                                -O K562_GABPA_BX_rep1_hg19.bam

# Merged IgG (K562-BX)
java -jar $PICARD MergeSamFiles -I 29104_IgG_i5006_K562_-_-_-_BX_FilteredBAM.bam \
                                -I 29246_IgG_i5006_K562_-_-_-_BX_FilteredBAM.bam \
                                -I 32091_IgG_i5006_K562_-_RPMI_-_BX_FilteredBAM.bam \
                                -I 32636_IgG_i5006_K562_-_-_Lysis-50Unuclease-10cycSonic-washesWithoutPVA-Exo-3splintLigation-shortExo-HCelution-16h-NoAmpure-HCday3-splintOligoswith3primeddC-40ulPCR-standardBuffer-62degreeAnnealing_BX_FilteredBAM.bam \
                                -I 32704_IgG_i5006_K562_-_IMDM_-_BX_FilteredBAM.bam \
                                -I 32110_IgG_i5006_K562_-_-_-_BX_FilteredBAM.bam \
                                -O K562_IgG_BX_merge.bam
# Merged IgG (K562-XO)
java -jar $PICARD MergeSamFiles -I 24140_IgG_i5006_K562_-_-_-_XO_FilteredBAM.bam \
                                -I 24262_IgG_i5006_K562_-_-_-_XO_FilteredBAM.bam \
                                -I 24430_IgG_i5006_K562_-_-_-_XO_FilteredBAM.bam \
                                -I 24643_IgG_i5006_K562_-_-_-_XO_FilteredBAM.bam \
                                -I 24965_IgG_i5006_K562_-_-_-_XO_FilteredBAM.bam \
                                -I 24969_IgG_i5006_K562_-_-_-_XO_FilteredBAM.bam \
                                -O K562_IgG_XO_merge.bam

# Use Picard to merge resequenced technical replicates, otherwise rename BAM using cp

# Input
cp 28326_Input_-_K562_-_-_50Unuclease10min-0cycSonic-newQuenchingbuffer_BI_FilteredBAM.bam K562_-_BI_rep1_hg19.bam
# H3K4me3
cp 28809_H3K4me3_ab8580_K562_-_-_50Unuclease10min-0cycSonic_BX_FilteredBAM.bam K562_H3K4me3_BX_rep1_hg19.bam
# H3K27ac
cp 28810_H3K27ac_ab4729_K562_-_-_50Unuclease10min-0cycSonic_BX_FilteredBAM.bam K562_H3K27ac_BX_rep1_hg19.bam
# PolII
cp 31646_PolII_ab76123_K562_-_RPMI_-_BX_FilteredBAM.bam K562_Pol2_BX_rep2_hg19.bam
# WDR5
cp 33918_WDR5_HPA047182_K562_-_IMDM_-_BX_FilteredBAM.bam K562_WDR5_BX_rep1_hg19.bam
cp 34606_WDR5_HPA047182_K562_-_-_-_BX_FilteredBAM.bam K562_WDR5_BX_rep2_hg19.bam
# CTCF
cp 24644_CTCF_07-729_K562_-_-_-_XO_FilteredBAM.bam K562_CTCF_XO_rep1_hg19.bam
cp 30478_CTCF_07-729_K562_-_-_Lysis-50Unuclease-4cycSonic-NoPapainDigestion-bothAmpureSteps_BX_FilteredBAM.bam K562_CTCF_BX_rep1_hg19.bam
cp 33924_CTCF_07-729_K562_-_IMDM_-_BX_FilteredBAM.bam K562_CTCF_BX_rep2_hg19.bam
# GABPA
cp 34607_GABPA_HPA003258_K562_-_-_-_BX_FilteredBAM.bam K562_GABPA_BX_rep2_hg19.bam
# RBBP5
cp 25084_RBBP5_A300-109A_K562_-_IMDM_-_XO_FilteredBAM.bam K562_RBBP5_XO_hg19.bam
cp 34600_RBBP5_A300-109A_K562_-_-_-_BX_FilteredBAM.bam K562_RBBP5_BX_rep2_hg19.bam
cp 34601_RBBP5_A300-109A_K562_-_-_-_BX_FilteredBAM.bam K562_RBBP5_BX_rep1_hg19.bam
# NFIA
cp 34598_NFIA_HPA008884_K562_-_-_-_BX_FilteredBAM.bam K562_NFIA_BX_rep2_hg19.bam
# SP1
cp 34604_Sp1_HPA001853_K562_-_-_-_BX_FilteredBAM.bam K562_SP1_BX_rep2_hg19.bam
# EP300
cp 32705_p300_A300-358A_K562_-_IMDM_-_BX_FilteredBAM.bam K562_EP300_BX_rep2_hg19.bam 
# E4F1
cp 25023_E4F1_A300-832A_K562_-_IMDM_-_XO_FilteredBAM.bam K562_E4F1_XO_rep1_hg19.bam
cp 34595_E4F1_A300-832A_K562_-_-_-_BX_FilteredBAM.bam K562_E4F1_BX_rep1_hg19.bam
cp 32241_E4F1_A300-832A_K562_-_-_-_BX_FilteredBAM.bam K562_E4F1_BX_rep2_hg19.bam
# MEIS2
cp 22624_MEIS2_MEIS2-1A11_K562_-_IMDM_-_XO_FilteredBAM.bam K562_MEIS2_XO_rep1_hg19.bam
cp 32096_MEIS2_1A11_K562_-_-_-_BX_FilteredBAM.bam K562_MEIS2_BX_rep1_hg19.bam
cp 31973_MEIS2_1A11_K562_-_-_-_BX_FilteredBAM.bam K562_MEIS2_BX_rep2_hg19.bam
# NRF1
cp 29445_NRF1_PCRP-NRF1-3H1-s_K562_-_-_-_BX_FilteredBAM.bam K562_NRF1_BX_rep2_hg19.bam
cp 31970_NRF1_3H1-s_K562_-_-_-_BX_FilteredBAM.bam K562_NRF1_BX_rep1_hg19.bam
cp 29445_NRF1_3H1-s_K562_-_-_-_BX_FilteredBAM.bam K562_NRF1_BX_rep2_hg19.bam
# SRF
cp 25104_SRF_A303-172A_K562_-_IMDM_-_XO_FilteredBAM.bam K562_SRF_XO_rep1_hg19.bam
cp 32849_SRF_A303-172A_K562_-_-_-_BX_FilteredBAM.bam K562_SRF_BX_rep1_hg19.bam
cp 34605_SRF_A303-172A_K562_-_-_-_BX_FilteredBAM.bam K562_SRF_BX_rep2_hg19.bam
# USF1
cp 24645_USF1_170427_K562_-_-_-_XO_FilteredBAM.bam K562_USF1_XO_rep1_hg19.bam
cp 29261_USF1_170427_K562_-_-_-_BX_FilteredBAM.bam K562_USF1_BX_rep1_hg19.bam
cp 32099_USF1_1B8_K562_-_-_-_BX_FilteredBAM.bam K562_USF1_BX_rep2_hg19.bam
# USF2
cp 28259_USF2_1A11_K562_-_IMDM_-_XO_FilteredBAM.bam K562_USF2_XO_rep1_hg19.bam
cp 32201_USF2_1A11_K562_-_-_-_BX_FilteredBAM.bam K562_USF2_BX_rep1_hg19.bam
cp 31972_USF2_1A11_K562_-_-_-_BX_FilteredBAM.bam K562_USF2_BX_rep2_hg19.bam
# YY1
cp 24263_YY1_1B2_K562_-_-_-_XO_FilteredBAM.bam K562_YY1_XO_rep1_hg19.bam
cp 32097_YY1_1B2_K562_-_-_-_BX_FilteredBAM.bam K562_YY1_BX_rep1_hg19.bam
cp 32101_YY1_1B2_K562_-_-_-_BX_FilteredBAM.bam K562_YY1_BX_rep2_hg19.bam

# Move renamed files into data/BAM
cd $WRK/../data
mv sample-BAM/K562_*.bam BAM/

# Index set of BAM files
for FILE in BAM/*.bam;
do
  [ -f $FILE.bai ] || samtools index $FILE
done

#!/bin/bash

# Script to hardcode the merging/renaming of PEGR BAM & MEME files into a standard file naming system

### CHANGE ME
WRK=/storage/group/bfp2/default/hxc585_HainingChen/Fox_NFIA_CTCF/00_Download_and_Preprocessing
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

# cp K562_H3K27ac_BX_rep1_hg19	28810
# cp K562_H3K4me3_BX_rep1_hg19	28809

# NFIA
java -jar $PICARD MergeSamFiles -I 32067_NFIA_HPA008884_K562_-_-_-_BX.bam \
                                -I 32116_NFIA_HPA008884_K562_-_-_-_BX.bam \
                                -O K562_NFIA_BX_rep1_hg19.bam
cp 34598_NFIA_HPA008884_K562_-_-_-_BX.bam K562_NFIA_BX_rep2_hg19.bam
# RAD21
java -jar $PICARD MergeSamFiles -I 32070_RAD21_HPA020044_K562_-_-_-_BX.bam \
                                -I 32119_RAD21_HPA020044_K562_-_-_-_BX.bam \
                                -O K562_RAD21_BX_rep1_hg19.bam
cp 34621_RAD21_HPA020044_K562_-_-_-_BX.bam K562_RAD21_BX_rep2_hg19.bam
# SMC3
cp 34622_Smc3_ab9263_K562_-_-_-_BX.bam K562_SMC3_BX_rep1_hg19.bam
cp 29518_Smc3_ab9263_K562_-_-_benzonase_XO.bam K562_SMC3_XO_rep1_hg19.bam

# FoxA1 HepG2, will add rep2 later
cp 32354_FoxA1_ab23738_HepG2_-_-_-_BX_FilteredBAM.bam HepG2_FoxA1_BX_rep1_hg19.bam
# FoxA1 K562 
# WAITING FOR DATA
# FoxA2 HepG2, need abcam antibody to be Rep1 and Rep2
cp  35655_FoxA2_HPA066846_HepG2_-_-_-_BX_FilteredBAM.bam HepG2_FoxA1_BX_rep3_hg19.bam

# FoxA2 K562, need data from rep2
cp  32659_FoxA2_ab256493_K562_-_-_-_BX_FilteredBAM.bam K562_FoxA1_BX_rep1_hg19.bam

# FoxA3 HepG2, need data rep2
java -jar $PICARD MergeSamFiles -I 36847_FOXA3_HPA054034_HepG2_-_IMDM_-_BX_FilteredBAM.bam \
                                -I 37123_FOXA3_HPA054034_HepG2_-_IMDM_-_BX_FilteredBAM.bam \
                                -I 37807_FOXA3_HPA054034_HepG2_-_IMDM_-_BX_FilteredBAM.bam \
                                -O K562_FOXA3_BX_rep1_hg19.bam
# Move renamed files into data/BAM
cd $WRK/../data
mv sample-BAM/K562_*.bam BAM/

# Index set of BAM files
for FILE in BAM/*.bam;
do
  [ -f $FILE.bai ] || samtools index $FILE
done

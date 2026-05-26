#!/bin/bash

# Script to hardcode the merging/renaming of PEGR BAM & MEME files into a standard file naming system

# Dependencies
# - java
# - samtools

set -exo
module load anaconda
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Inputs and outputs
BAMDIR=../data/BAM

# Script shortcuts
PICARD=../bin/picard.jar

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


java -jar $PICARD MergeSamFiles -O HepG2_IgG_BX_merge_hg38.bam \
        -I 38481_IgG_i5006_HepG2_-_-_-_BX.bam \
        -I 38482_IgG_i5006_HepG2_-_-_-_BX.bam \
        -I 38483_IgG_i5006_HepG2_-_-_-_BX.bam

# Use Picard to merge resequenced technical replicates, otherwise rename BAM using cp

# CTCF
cp 33924_CTCF_07-729_K562_-_IMDM_-_BX.bam K562_CTCF_BX_rep1_hg38.bam
cp 34174_CTCF_07-729_K562_-_IMDM_-_BX.bam K562_CTCF_BX_rep2_hg38.bam
cp 41780_CTCF_07-729_SEM_-_RPMI_-_BX.bam SEM_CTCF_BX_WT_hg38.bam
cp 41781_CTCF_07-729_SEM_-_RPMI_CTCF depletion_BX.bam SEM_CTCF_BX_CTCTKO_hg38.bam

# CTCF (non-crosslinked)
cp 29009_CTCF_07-729_K562_-_-_2stepNonXLlysis-50Unuclease-10min-0cycSonic-quenchONLY_BX.bam K562_CTCF_nonXLBX_rep1_hg38.bam
cp 29010_CTCF_07-729_K562_-_-_2stepNonXLlysis-50Unuclease-10min-0cycSonic-quenchONLY_BX.bam K562_CTCF_nonXLBX_rep2_hg38.bam

# CTCF (KO in SEM cell)
cp 41780_CTCF_07-729_SEM_-_RPMI_-_BX_hg38.bam SEM-CTCFAID2-DMSO_CTCF_BX_rep1_hg38.bam
cp 41781_CTCF_07-729_SEM_-_RPMI_CTCFdepletion_BX_hg38.bam SEM-CTCFAID2-1µM5PhIAA_CTCF_BX_rep1_hg38.bam
cp 42078_CTCF_07-729_SEM_-_RPMI_-_BX_hg38.bam SEM-CTCFAID2-DMSO_CTCF_BX_rep2_hg38.bam
cp 42079_CTCF_07-729_SEM_-_RPMI_CTCFdepletion_BX_hg38.bam SEM-CTCFAID2-1µM5PhIAA_CTCF_BX_rep2_hg38.bam

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
cp 32354_FoxA1_ab23738_HepG2_-_-_-_BX_hg38.bam HepG2_FOXA1_BX_rep1_hg38.bam
cp 33219_FOXA1_ab23738_HepG2_-_-_-_BX_hg38.bam HepG2_FOXA1_BX_rep2_hg38.bam

# FOXA1 (K562)
cp 38477_FoxA1_ab23738_K562_-_IMDM_-_BX_hg38.bam K562_FOXA1_BX_rep1_hg38.bam
cp 38478_FoxA1_ab23738_K562_-_IMDM_-_BX_hg38.bam K562_FOXA1_BX_rep2_hg38.bam

# FOXA2 (HepG2)
cp 32467_FOXA2_ab256493_HepG2_-_-_-_BX_hg38.bam HepG2_FOXA2_BX_rep2_hg38.bam
cp 38480_FOXA2_ab256493_HepG2_-_-_-_BX_hg38.bam HepG2_FOXA2_BX_rep1_hg38.bam

# HNF4A (HepG2)
cp 34158_HNF4A_HPA004712_HepG2_-_-_-_BX.bam HepG2_HNF4A_BX_rep1_hg38.bam

# FOXA2 (K562)
cp 32659_FOXA2_ab256493_K562_-_-_-_BX.bam K562_FOXA2_BX_rep1_hg38.bam
cp 38479_FOXA2_ab256493_K562_-_IMDM_-_BX.bam K562_FOXA2_BX_rep2_hg38.bam

# Native K562 and HepG2 combination assay
cp 41762_FoxA1_ab23738_K562_-_IMDM_-_BX_hg38.bam   KH_FOXA1_BX_rep2_hg38.bam
cp 41965_FoxA1_ab23738_K562_-_IMDM_-_BX_hg38.bam   KH_FOXA1_BX_rep1_hg38.bam
cp 41978_IgG_i5006_K562_-_IMDM_-_BX_hg38.bam       KH_IgG_BX_rep1_hg38.bam

cp 42304_IgG_i5006_K562_-_IMDM_-_BX_hg38.bam       SWCK562_IgG_BX_rep1_hg38.bam
cp 41764_FoxA1_ab23738_K562_-_IMDM_-_BX_hg38.bam   SWCK562_FOXA1_BX_rep1_hg38.bam
cp 42581_FoxA1_ab23738_K562_-_IMDM_-_BX_hg38.bam   SWCK562_FOXA1_BX_rep2_hg38.bam

cp 39314_FoxA1_ab23738_K562_-_IMDM_-_BX_hg38.bam   nativeK562_FoxA1_BX_rep2_hg38.bam
cp 39487_FoxA1_ab23738_K562_-_IMDM_-_BX_hg38.bam   nativeK562_FoxA1_BX_rep1_hg38.bam
cp 38819_IgG_i5006_K562_-_IMDM_-_BX_hg38.bam       nativeK562_IgG_BX_rep1_hg38.bam

cp 39490_FoxA1_ab23738_HepG2_-_IMDM_-_BX_hg38.bam  nativeHepG2_FoxA1_BX_rep1_hg38.bam
cp 39899_FoxA1_ab23738_HepG2_-_IMDM_-_BX_hg38.bam  nativeHepG2_FoxA1_BX_rep2_hg38.bam
cp 39492_IgG_i5006_HepG2_-_IMDM_-_BX_hg38.bam      nativeHepG2_IgG_BX_rep1_hg38.bam

cp 41974_IgG_i5006_K562_-_IMDM_-_BX_hg38.bam       KK_IgG_BX_rep1_hg38.bam
cp 41760_FoxA1_ab23738_K562_-_IMDM_-_BX_hg38.bam   KK_FOXA1_BX_rep2_hg38.bam
cp 41971_FoxA1_ab23738_K562_-_IMDM_-_BX_hg38.bam   KK_FOXA1_BX_rep1_hg38.bam

cp 41963_FoxA1_ab23738_K562_-_IMDM_-_BX_hg38.bam   WCEK562_FOXA1_BX_rep1_hg38.bam

## naked DNA 
cp 40726_Input_-_K562_-_IMDM_-_BI_hg38_hg38.bam    NakedDNA_-_0.04U-10min-BI_hg38.bam
cp 40727_Input_-_K562_-_IMDM_-_BI_hg38_hg38.bam    NakedDNA_-_0.125U-10min-BI_hg38.bam
cp 40728_Input_-_K562_-_IMDM_-_BI_hg38_hg38.bam    NakedDNA_-_0.375U-10min-BI_hg38.bam

cd ../data
mv sample-BAM/K562_*.bam $BAMDIR
mv sample-BAM/HepG2_*.bam $BAMDIR
mv sample-BAM/NakedDNA_*.bam $BAMDIR
mv sample-BAM/KK_*.bam $BAMDIR
mv sample-BAM/WCEK562*.bam $BAMDIR
mv sample-BAM/nativeHepG2_I*.bam $BAMDIR

# Index set of BAM files
for FILE in $BAMDIR/*.bam;
do
  [ -f $FILE.bai ] || samtools index $FILE
done

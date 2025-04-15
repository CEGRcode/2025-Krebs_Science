#!/bin/bash

# Make PE insert size histogram of BNase-seq 2 replicates at each of three digestion levels (50U-3min, 50U-10min, 50U-30min)

### CHANGE ME
WRK=/storage/group/bfp2/default/hxc585_HainingChen/2025_Chen_TF-Nuc/Z_Figures
#WRK=/ocean/projects/see180003p/owlang/2024-Krebs_Science/Z_Figures
#WRK=/storage/home/owl5022/scratch/2024-Krebs_Science/Z_Figures
###

# Dependencies
# - java
# - pandas
# - python
# - seaborn

set -exo
module load anaconda3
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Inputs and outputs
BAMDIR=$WRK/../data/BAM

# Script shortcuts
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar
HISTOGRAM=$WRK/../bin/make_fragment_histograms1.py

[ -d S1/a ] || mkdir -p S1/a
#enzyme cut nucleotide analysis
cp $WRK/../0X_Bulk_Processing/Library/UpstreamKmerAnalysis/DiTally/DNase-seq_ENCFF518XTC_rep1_hg38_SUBSAMPLE_NT-l50r100-R1.svg  S1/a
cp $WRK/../0X_Bulk_Processing/Library/UpstreamKmerAnalysis/DiTally/MNase-seq_21U_rep1_hg38_SUBSAMPLE_NT-l50r100-R1.svg S1/a
cp $WRK/../0X_Bulk_Processing/Library/UpstreamKmerAnalysis/DiTally/BNase-seq_50U-10min_merge_hg38_SUBSAMPLE_NT-l50r100-R1.svg S1/a
cp $WRK/../0X_Bulk_Processing/Library/UpstreamKmerAnalysis/DiTallyMPE-seq_20min_rep2_mm10_SUBSAMPLE_NT-l50r100-R1.svg S1/a

[ -d S1/b ] || mkdir -p S1/b

# Run Paired-end Statistics
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/BNase-seq_50U-3min_1_hg38.bam -o S1/b/BNase-seq_50U-3min_1_hg38
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/BNase-seq_50U-3min_2_hg38.bam -o S1/b/BNase-seq_50U-3min_2_hg38
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/BNase-seq_50U-10min_1_hg38.bam -o S1/b/BNase-seq_50U-10min_1_hg38
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/BNase-seq_50U-10min_2_hg38.bam -o S1/b/BNase-seq_50U-10min_2_hg38
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/BNase-seq_50U-30min_1_hg38.bam -o S1/b/BNase-seq_50U-30min_1_hg38
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/BNase-seq_50U-30min_2_hg38.bam -o S1/b/BNase-seq_50U-30min_2_hg38

# Generate insert size frequency histograms
python $HISTOGRAM --ymax 100000 --reflines 30 40 50 60 70 80 90 100 110 120 130 140 150 160  -i S1/b/BNase-seq_50U-3min_1_hg38_InsertHistogram.out -o S1/b/BNase-seq_50U-3min_1_hg38_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax 100000 --reflines 30 40 50 60 70 80 90 100 110 120 130 140 150 160  -i S1/b/BNase-seq_50U-3min_2_hg38_InsertHistogram.out -o S1/b/BNase-seq_50U-3min_2_hg38_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax 400000 --reflines 30 40 50 60 70 80 90 100 110 120 130 140 150 160  -i S1/b/BNase-seq_50U-10min_1_hg38_InsertHistogram.out -o S1/b/BNase-seq_50U-10min_1_hg38_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax 5100000 --reflines 30 40 50 60 70 80 90 100 110 120 130 140 150 160  -i S1/b/BNase-seq_50U-10min_2_hg38_InsertHistogram.out -o S1/b/BNase-seq_50U-10min_2_hg38_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax 140000 --reflines 30 40 50 60 70 80 90 100 110 120 130 140 150 160  -i S1/b/BNase-seq_50U-30min_1_hg38_InsertHistogram.out -o S1/b/BNase-seq_50U-30min_1_hg38_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax 140000 --reflines 30 40 50 60 70 80 90 100 110 120 130 140 150 160  -i S1/b/BNase-seq_50U-30min_2_hg38_InsertHistogram.out -o S1/b/BNase-seq_50U-30min_2_hg38_InsertHistogram_HIST.svg

[ -d S1/c ] || mkdir -p S1/c

# Run Paired-end Statistics
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/NakedDNA_BNase-seq_0.04U_1_hg38.bam -o S1/c/NakedDNA_BNase-seq_0.04U_1_hg38
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/NakedDNA_BNase-seq_0.125U_2_hg38.bam -o S1/c/NakedDNA_BNase-seq_0.125U_2_hg38
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/NakedDNA_BNase-seq_0.375U_3_hg38.bam -o S/c/NakedDNA_BNase-seq_0.375U_3_hg38
# Generate insert size frequency histograms
python $HISTOGRAM --ymax 200000  -i S1/c/NakedDNA_BNase-seq_0.04U_1_hg38_InsertHistogram.out -o S1/c/BNase-seq_50U-3min_1_hg38_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax 60000   -i S1/c/NakedDNA_BNase-seq_0.125U_2_hg38_InsertHistogram.out -o S1/c/BNase-seq_50U-3min_2_hg38_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax 15000  -i S1/c//NakedDNA_BNase-seq_0.375U_3_hg38_InsertHistogram.out -o S1/c/BNase-seq_50U-10min_1_hg38_InsertHistogram_HIST.svg

[ -d S1/d ] || mkdir -p S1/d

# Run Paired-end Statistics
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 440 $BAMDIR/DNase-FLASH_SRR801880.bam -o S1/d/DNase-FLASH_SRR801880
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 440 $BAMDIR/DNase-FLASH_SRR801881.bam -o S1/d/DNase-FLASH_SRR801881
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 440 $BAMDIR/MPE-seq_10min_merge_rep1_mm10.bam -o S1/d/MPE-seq_10min_merge_rep1_mm10
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 440 $BAMDIR/MPE-seq_20min_merge_rep1_mm10.bam -o S1/d/MPE-seq_20min_merge_rep1_mm10
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 440 $BAMDIR/MPE-seq_30min_merge_rep1_mm10.bam -o S1/d/MPE-seq_30min_merge_rep1_mm10
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 440 $BAMDIR/MPE-seq_10min_rep2_mm10.bam -o S1/d/MPE-seq_10min_rep2_mm10
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 440 $BAMDIR/MPE-seq_20min_rep2_mm10.bam -o S1/d/MPE-seq_20min_rep2_mm10
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 440 $BAMDIR/MPE-seq_30min_rep2_mm10.bam -o S1/d/MPE-seq_30min_rep2_mm10

# Generate insert size frequency histograms
python $HISTOGRAM --ymax 1800000 --reflines 30 40 50 60 70 80 90 100 110 120 130  -i S1/d/DNase-FLASH_SRR801880_InsertHistogram.out -o S1/d/DNase-FLASH_SRR801880_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax 1400000 --reflines 30 40 50 60 70 80 90 100 110 120 130  -i S1/d/DNase-FLASH_SRR801881_InsertHistogram.out -o S1/d/DNase-FLASH_SRR801881_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax 500000 --reflines 80 90 100 110 120 130  -i S1/d/MPE-seq_10min_merge_rep1_mm10_InsertHistogram.out -o S1/d/MPE-seq_10min_merge_rep1_mm10_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax 170000 --reflines 80 90 100 110 120 130  -i S1/d/MPE-seq_10min_rep2_mm10_InsertHistogram.out -o S1/d/MPE-seq_10min_rep2_mm10_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax 500000 --reflines 80 90 100 110 120 130  -i S1/d/MPE-seq_20min_merge_rep1_mm10_InsertHistogram.out -o S1/d/MPE-seq_20min_merge_rep1_mm10_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax 200000 --reflines 60 70 80 90 100 110 120 130 140 150  -i S1/d/MPE-seq_20min_rep2_mm10_InsertHistogram.out -o S1/d/MPE-seq_20min_rep2_mm10_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax 400000 --reflines 80 90 100 110 120 130  -i S1/d/MPE-seq_30min_merge_rep1_mm10_InsertHistogram.out -o S1/d/MPE-seq_30min_merge_rep1_mm10_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax 120000 --reflines 90 100 110 120 130 140 150 160  -i S1/d/MPE-seq_30min_rep2_mm10_InsertHistogram.out -o S1/d/MPE-seq_30min_rep2_mm10_InsertHistogram_HIST.svg
#!/bin/bash

# Dependencies
# - java
# - pandas
# - python
# - seaborn

set -exo
source activate bx

# Inputs and outputs
BAMDIR=../data/BAM

# Script shortcuts
SCRIPTMANAGER=../bin/ScriptManager-v0.15.jar
HISTOGRAM=../bin/make_fragment_histograms.py

[ -d F1 ] || mkdir F1

# ===============================================================================================================================

# Figure 1A-B are graphics

# ===============================================================================================================================

[ -d F1/C ] || mkdir -p F1/C

# Run Paired-end Statistics
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/BNase-seq_50U-10min_merge_hg38.bam -o F1/C/BNase-seq_50U-10min_merge_hg38
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/DNase-seq_ENCFF518XTC_rep1_hg38.bam -o F1/C/DNase-seq_ENCFF518XTC_rep1_hg38
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/MNase-seq_21U_rep1_hg38.bam -o F1/C/MNase-seq_21U_rep1_hg38
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/MNase-seq_304U_rep1_hg38.bam -o F1/C/MNase-seq_304U_rep1_hg38

# Generate insert size frequency histograms
python $HISTOGRAM --ymax 16000000 --reflines 30 40 50 60 70 80 90 100 110 120 130 140 150 160 -i F1/C/BNase-seq_50U-10min_merge_hg38_InsertHistogram.out -o F1/C/BNase-seq_50U-10min_merge_hg38_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax  1200000 --reflines 70 80 90 100 110 120 130 140 -i F1/C/DNase-seq_ENCFF518XTC_rep1_hg38_InsertHistogram.out -o F1/C/DNase-seq_ENCFF518XTC_rep1_hg38_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax  1400000 --reflines 106 113 128 134 149 -i F1/C/MNase-seq_21U_rep1_hg38_InsertHistogram.out -o F1/C/MNase-seq_21U_rep1_hg38_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax  1000000 --reflines 106 127 145 -i F1/C/MNase-seq_304U_rep1_hg38_InsertHistogram.out -o F1/C/MNase-seq_304U_rep1_hg38_InsertHistogram_HIST.svg
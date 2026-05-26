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
HGENOME=../data/hg38_files/hg38.fa
MGENOME=../data/mm10_files/mm10.fa

# Index genome if doesnt exist
[ -f $HGENOME.fai ] || samtools faidx $HGENOME
[ -f $MGENOME.fai ] || samtools faidx $MGENOME

# Script shortcuts
SCRIPTMANAGER=../bin/ScriptManager-v0.15.jar
HISTOGRAM=../bin/make_fragment_histograms.py
UPDOWN_KMER=../bin/updownstream_di-nt_tally.py
KMER2NT=../bin/dint_to_nt_positional_count_matrix.py
STACKNT=../bin/make_stack_barchart_TSV.py
SPLIT_KMER=../bin/split_kmer_to_composite.py

[ -d S1 ] || mkdir S1

# ===============================================================================================================================

[ -d S1/A ] || mkdir S1/A

## Perform sequence content analysis around cut sites
TEMP=temp-S1A
[ -d $TEMP ] || mkdir $TEMP

# # Get BAM info
# samtools flagstat $BAMDIR/BNase-seq_50U-10min_merge_hg38.bam  > $TEMP/BNase-seq_50U-10min_merge_hg38.bam.flagstat
# samtools flagstat $BAMDIR/DNase-seq_ENCFF518XTC_rep1_hg38.bam > $TEMP/DNase-seq_ENCFF518XTC_rep1_hg38.bam.flagstat
# samtools flagstat $BAMDIR/MNase-seq_21U_rep1_hg38.bam         > $TEMP/MNase-seq_21U_rep1_hg38.bam.flagstat
# samtools flagstat $BAMDIR/MPE-seq_20min_rep2_mm10.bam         > $TEMP/MPE-seq_20min_rep2_mm10.bam.flagstat

# Subsample BAM - hardcoded this for different BAM files (Seed=2)
samtools view -b -s 2.0001 ${BAMDIR}/BNase-seq_50U-10min_merge_hg38.bam  > $TEMP/BNase-seq_50U-10min_merge_hg38_DOWNSAMPLE.bam   # ~1.3B --> ~125K
samtools view -b -s 2.001  ${BAMDIR}/DNase-seq_ENCFF518XTC_rep1_hg38.bam > $TEMP/DNase-seq_ENCFF518XTC_rep1_hg38_DOWNSAMPLE.bam  # ~100M --> ~100K
samtools view -b -s 2.001  ${BAMDIR}/MNase-seq_21U_rep1_hg38.bam         > $TEMP/MNase-seq_21U_rep1_hg38_DOWNSAMPLE.bam          # ~170M --> ~170K
samtools view -b -s 2.002  ${BAMDIR}/MPE-seq_20min_rep2_mm10.bam         > $TEMP/MPE-seq_20min_rep2_mm10_DOWNSAMPLE.bam          #  ~55M --> ~110K

# Index
samtools index $TEMP/BNase-seq_50U-10min_merge_hg38_DOWNSAMPLE.bam
samtools index $TEMP/DNase-seq_ENCFF518XTC_rep1_hg38_DOWNSAMPLE.bam
samtools index $TEMP/MNase-seq_21U_rep1_hg38_DOWNSAMPLE.bam
samtools index $TEMP/MPE-seq_20min_rep2_mm10_DOWNSAMPLE.bam

# # Get subsample BAM info
# samtools flagstat $TEMP/BNase-seq_50U-10min_merge_hg38_DOWNSAMPLE.bam  > $TEMP/BNase-seq_50U-10min_merge_hg38_DOWNSAMPLE.bam.flagstat
# samtools flagstat $TEMP/DNase-seq_ENCFF518XTC_rep1_hg38_DOWNSAMPLE.bam > $TEMP/DNase-seq_ENCFF518XTC_rep1_hg38_DOWNSAMPLE.bam.flagstat
# samtools flagstat $TEMP/MNase-seq_21U_rep1_hg38_DOWNSAMPLE.bam         > $TEMP/MNase-seq_21U_rep1_hg38_DOWNSAMPLE.bam.flagstat
# samtools flagstat $TEMP/MPE-seq_20min_rep2_mm10_DOWNSAMPLE.bam         > $TEMP/MPE-seq_20min_rep2_mm10_DOWNSAMPLE.bam.flagstat

# Perform upstream kmer analysis around 5' cut sites (di-nucleotides, -50 to +100)
python $UPDOWN_KMER -l 50 -r 100 --read1 -p -g $HGENOME -i $TEMP/BNase-seq_50U-10min_merge_hg38_DOWNSAMPLE.bam  -o S1/A/BNase-seq_50U-10min_merge_hg38_DINT-l50r100-R1.tsv
python $UPDOWN_KMER -l 50 -r 100 --read1 -p -g $HGENOME -i $TEMP/DNase-seq_ENCFF518XTC_rep1_hg38_DOWNSAMPLE.bam -o S1/A/DNase-seq_ENCFF518XTC_rep1_hg38_DINT-l50r100-R1.tsv
python $UPDOWN_KMER -l 50 -r 100 --read1 -p -g $HGENOME -i $TEMP/MNase-seq_21U_rep1_hg38_DOWNSAMPLE.bam         -o S1/A/MNase-seq_21U_rep1_hg38_DINT-l50r100-R1.tsv
python $UPDOWN_KMER -l 50 -r 100 --read1 -p -g $MGENOME -i $TEMP/MPE-seq_20min_rep2_mm10_DOWNSAMPLE.bam         -o S1/A/MPE-seq_20min_rep2_mm10_DINT-l50r100-R1.tsv

# Re-tally for single nucleotide counts
python $KMER2NT -i S1/A/BNase-seq_50U-10min_merge_hg38_DINT-l50r100-R1.tsv  -o S1/A/BNase-seq_50U-10min_merge_hg38_NT-l50r100-R1.tsv
python $KMER2NT -i S1/A/DNase-seq_ENCFF518XTC_rep1_hg38_DINT-l50r100-R1.tsv -o S1/A/DNase-seq_ENCFF518XTC_rep1_hg38_NT-l50r100-R1.tsv
python $KMER2NT -i S1/A/MNase-seq_21U_rep1_hg38_DINT-l50r100-R1.tsv         -o S1/A/MNase-seq_21U_rep1_hg38_NT-l50r100-R1.tsv
python $KMER2NT -i S1/A/MPE-seq_20min_rep2_mm10_DINT-l50r100-R1.tsv         -o S1/A/MPE-seq_20min_rep2_mm10_NT-l50r100-R1.tsv

# Generate Figure: stack single nucleotides (enforce same-frequency)
python $STACKNT --entropy -i <(cut -f1,40-71 S1/A/BNase-seq_50U-10min_merge_hg38_NT-l50r100-R1.tsv)  --title BNase-seq_50U-10min_merge_hg38  -o S1/A/BNase-seq_50U-10min_merge_hg38_NT-l50r100-R1.svg
python $STACKNT --entropy -i <(cut -f1,40-71 S1/A/DNase-seq_ENCFF518XTC_rep1_hg38_NT-l50r100-R1.tsv) --title DNase-seq_ENCFF518XTC_rep1_hg38 -o S1/A/DNase-seq_ENCFF518XTC_rep1_hg38_NT-l50r100-R1.svg
python $STACKNT --entropy -i <(cut -f1,40-71 S1/A/MNase-seq_21U_rep1_hg38_NT-l50r100-R1.tsv)         --title MNase-seq_21U_rep1_hg38         -o S1/A/MNase-seq_21U_rep1_hg38_NT-l50r100-R1.svg
python $STACKNT --entropy -i <(cut -f1,40-71 S1/A/MPE-seq_20min_rep2_mm10_NT-l50r100-R1.tsv)         --title MPE-seq_20min_rep2_mm10         -o S1/A/MPE-seq_20min_rep2_mm10_NT-l50r100-R1.svg

# ===============================================================================================================================

[ -d S1/B ] || mkdir -p S1/B

# Run Paired-end Statistics
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/BNase-seq_50U-3min_1_hg38.bam -o S1/B/BNase-seq_50U-3min_1_hg38
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/BNase-seq_50U-3min_2_hg38.bam -o S1/B/BNase-seq_50U-3min_2_hg38
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/BNase-seq_50U-10min_1_hg38.bam -o S1/B/BNase-seq_50U-10min_1_hg38
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/BNase-seq_50U-10min_2_hg38.bam -o S1/B/BNase-seq_50U-10min_2_hg38
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/BNase-seq_50U-30min_1_hg38.bam -o S1/B/BNase-seq_50U-30min_1_hg38
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/BNase-seq_50U-30min_2_hg38.bam -o S1/B/BNase-seq_50U-30min_2_hg38

# Generate insert size frequency histograms
python $HISTOGRAM --ymax 100000 --reflines 30 40 50 60 70 80 90 100 110 120 130 140 150 160  -i S1/B/BNase-seq_50U-3min_1_hg38_InsertHistogram.out -o S1/B/BNase-seq_50U-3min_1_hg38_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax 100000 --reflines 30 40 50 60 70 80 90 100 110 120 130 140 150 160  -i S1/B/BNase-seq_50U-3min_2_hg38_InsertHistogram.out -o S1/B/BNase-seq_50U-3min_2_hg38_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax 400000 --reflines 30 40 50 60 70 80 90 100 110 120 130 140 150 160  -i S1/B/BNase-seq_50U-10min_1_hg38_InsertHistogram.out -o S1/B/BNase-seq_50U-10min_1_hg38_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax 5100000 --reflines 30 40 50 60 70 80 90 100 110 120 130 140 150 160  -i S1/B/BNase-seq_50U-10min_2_hg38_InsertHistogram.out -o S1/B/BNase-seq_50U-10min_2_hg38_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax 140000 --reflines 30 40 50 60 70 80 90 100 110 120 130 140 150 160  -i S1/B/BNase-seq_50U-30min_1_hg38_InsertHistogram.out -o S1/B/BNase-seq_50U-30min_1_hg38_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax 140000 --reflines 30 40 50 60 70 80 90 100 110 120 130 140 150 160  -i S1/B/BNase-seq_50U-30min_2_hg38_InsertHistogram.out -o S1/B/BNase-seq_50U-30min_2_hg38_InsertHistogram_HIST.svg

# ===============================================================================================================================

[ -d S1/C ] || mkdir -p S1/C

# Run Paired-end Statistics
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/NakedDNA_BNase-seq_0.04U-10min-BI_hg38.bam  -o S1/C/NakedDNA_BNase-seq_0.04U-10min-BI_hg38
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/NakedDNA_BNase-seq_0.125U-10min-BI_hg38.bam -o S1/C/NakedDNA_BNase-seq_0.125U-10min-BI_hg38
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/NakedDNA_BNase-seq_0.375U-10min-BI_hg38.bam -o S1/C/NakedDNA_BNase-seq_0.375U-10min-BI_hg38
# Generate insert size frequency histograms
python $HISTOGRAM --ymax 200000 -i S1/C/NakedDNA_BNase-seq_0.04U-10min-BI_hg38_InsertHistogram.out  -o S1/C/NakedDNA_BNase-seq_0.04U-10min-BI_hg38_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax 60000  -i S1/C/NakedDNA_BNase-seq_0.125U-10min-BI_hg38_InsertHistogram.out  -o S1/C/NakedDNA_BNase-seq_0.125U-10min-BI_hg38_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax 15000  -i S1/C//NakedDNA_BNase-seq_0.375U-10min-BI_hg38_InsertHistogram.out -o S1/C/NakedDNA_BNase-seq_0.375U-10min-BI_hg38_InsertHistogram_HIST.svg

# ===============================================================================================================================

[ -d S1/D ] || mkdir -p S1/D

# Run Paired-end Statistics
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 440 $BAMDIR/DNase-FLASH_SRR801880.bam -o S1/D/DNase-FLASH_SRR801880
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 440 $BAMDIR/DNase-FLASH_SRR801881.bam -o S1/D/DNase-FLASH_SRR801881
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 440 $BAMDIR/MPE-seq_10min_merge_rep1_mm10.bam -o S1/D/MPE-seq_10min_merge_rep1_mm10
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 440 $BAMDIR/MPE-seq_20min_merge_rep1_mm10.bam -o S1/D/MPE-seq_20min_merge_rep1_mm10
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 440 $BAMDIR/MPE-seq_30min_merge_rep1_mm10.bam -o S1/D/MPE-seq_30min_merge_rep1_mm10
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 440 $BAMDIR/MPE-seq_10min_rep2_mm10.bam -o S1/D/MPE-seq_10min_rep2_mm10
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 440 $BAMDIR/MPE-seq_20min_rep2_mm10.bam -o S1/D/MPE-seq_20min_rep2_mm10
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 440 $BAMDIR/MPE-seq_30min_rep2_mm10.bam -o S1/D/MPE-seq_30min_rep2_mm10

# Generate insert size frequency histograms
python $HISTOGRAM --ymax 1800000 --xmax-440 --reflines 30 40 50 60 70 80 90 100 110 120 130  -i S1/D/DNase-FLASH_SRR801880_InsertHistogram.out -o S1/D/DNase-FLASH_SRR801880_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax 1400000 --xmax-440 --reflines 30 40 50 60 70 80 90 100 110 120 130  -i S1/D/DNase-FLASH_SRR801881_InsertHistogram.out -o S1/D/DNase-FLASH_SRR801881_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax 500000  --xmax-440 --reflines 80 90 100 110 120 130  -i S1/D/MPE-seq_10min_merge_rep1_mm10_InsertHistogram.out -o S1/D/MPE-seq_10min_merge_rep1_mm10_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax 170000  --xmax-440 --reflines 80 90 100 110 120 130  -i S1/D/MPE-seq_10min_rep2_mm10_InsertHistogram.out -o S1/D/MPE-seq_10min_rep2_mm10_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax 500000  --xmax-440 --reflines 80 90 100 110 120 130  -i S1/D/MPE-seq_20min_merge_rep1_mm10_InsertHistogram.out -o S1/D/MPE-seq_20min_merge_rep1_mm10_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax 200000  --xmax-440 --reflines 60 70 80 90 100 110 120 130 140 150  -i S1/D/MPE-seq_20min_rep2_mm10_InsertHistogram.out -o S1/D/MPE-seq_20min_rep2_mm10_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax 400000  --xmax-440 --reflines 80 90 100 110 120 130  -i S1/D/MPE-seq_30min_merge_rep1_mm10_InsertHistogram.out -o S1/D/MPE-seq_30min_merge_rep1_mm10_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax 120000  --xmax-440 --reflines 90 100 110 120 130 140 150 160  -i S1/D/MPE-seq_30min_rep2_mm10_InsertHistogram.out -o S1/D/MPE-seq_30min_rep2_mm10_InsertHistogram_HIST.svg
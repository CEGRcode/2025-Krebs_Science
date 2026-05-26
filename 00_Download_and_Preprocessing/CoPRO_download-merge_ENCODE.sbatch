#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=48gb
#SBATCH -t 5:00:00
#SBATCH -A open
#SBATCH -o logs/CoPRO_download-merge_ENCODE.log.out
#SBATCH -e logs/CoPRO_download-merge_ENCODE.log.err

# Download capped CoPRO data from ENCODE and merge reps

# data/BAM
#   |--CoPRO_Capped_merge_hg38.bam

# Dependencies
# - java
# - samtools
# - wget

set -exo
source activate bx

# Inputs and outputs
BAMDIR=../data/BAM

# Script shortcuts
PICARD=../bin/picard.jar

[ -d ENCODE ] || mkdir ENCODE

# Download BAM files
wget -c -O ENCODE/ENCFF663UAN.bam https://www.encodeproject.org/files/ENCFF663UAN/@@download/ENCFF663UAN.bam
wget -c -O ENCODE/ENCFF877SQU.bam https://www.encodeproject.org/files/ENCFF877SQU/@@download/ENCFF877SQU.bam

# Merge BAM files
java -jar $PICARD MergeSamFiles -O $BAMDIR/CoPRO_Capped_merge_hg38.bam \
	-I ENCODE/ENCFF663UAN.bam \
	-I ENCODE/ENCFF877SQU.bam

# Index
samtools index $BAMDIR/CoPRO_Capped_merge_hg38.bam

# Clean-up
rm -r ENCODE
#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=48gb
#SBATCH -t 5:00:00
#SBATCH -A open
#SBATCH -o logs/ATAC-seq_download_ENCODE.log.out
#SBATCH -e logs/ATAC-seq_download_ENCODE.log.err

# Download and rename the ENCODE ATAC-seq BAM file

# data/BAM
#   |--ATAC-seq_ENCFF077FBI_hg38.bam

# Dependencies
# - samtools
# - wget

set -exo
source activate bx

# Inputs and outputs
BAMDIR=../data/BAM

# Download BAM files
wget -c -O $BAMDIR/ATAC-seq_ENCFF077FBI_hg38.bam https://www.encodeproject.org/files/ENCFF077FBI/@@download/ENCFF077FBI.bam

# Index
samtools index $BAMDIR/ATAC-seq_ENCFF077FBI_hg38.bam

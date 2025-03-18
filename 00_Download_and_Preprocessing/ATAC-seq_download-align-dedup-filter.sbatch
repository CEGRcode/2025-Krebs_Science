#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=48gb
#SBATCH -t 5:00:00
#SBATCH -A open
#SBATCH -o logs/MNase-seq-ENCODE_download-align-filter-merge.log.out
#SBATCH -e logs/MNase-seq-ENCODE_download-align-filter-merge.log.err

# Download deep MNase sequencing (single-end, ABI SOLiD) and merge from ENCODE
# No need to MarkDuplicates b/c ENCODE says these are bio replicates (and
# pipeline already filters)

# data/BAM
#   |--MNase-seq_ENCODE_merge_hg38.bam

### CHANGE ME
WRK=WRK=/storage/group/bfp2/default/hxc585_HainingChen/2025_Chen_TF-Nuc/00_Download_and_Preprocessing
WRK=/storage/home/owl5022/scratch/2024-Krebs_Science/00_Download_and_Preprocessing
THREADS=8
###

# Dependencies
# - bowtie v1.2.3
# - java
# - samtools
# - wget

set -exo
module load samtools
module load anaconda
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/colorspace
cd $WRK

# Inputs and outputs
BAMDIR=$WRK/../data/BAM
GENOME=$WRK/../data/hg38_files/hg38.fa.colorspace

# Script shortcuts
PICARD=$WRK/../bin/picard.jar

[ -d ENCODE ] || mkdir ENCODE

# Download BAM files
wget -c -O $BAMDIR/ATAC-seq_ENCFF077FBI.bam https://www.encodeproject.org/files/ENCFF077FBI/@@download/ENCFF077FBI.bam

# Index
samtools index $BAMDIR/MNase-seq_ENCODE_merge_hg38.bam

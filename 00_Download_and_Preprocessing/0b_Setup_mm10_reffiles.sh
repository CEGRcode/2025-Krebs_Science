#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=24gb
#SBATCH -t 6:00:00
#SBATCH -A open
#SBATCH -o logs/0b_Setup_mm10_reffiles.log.out
#SBATCH -e logs/0b_Setup_mm10_reffiles.log.err

# Download mouse (mm10) ref genome and build alignment indexes
# - primarily just used for MPE-seq data processing

set -exo
source activate bx

GENOME=../data/mm10_files/mm10.fa

[ -d ../data/mm10_files ] || mkdir -p ../data/mm10_files

# Download genome FASTA
wget -O $GENOME.gz https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
gzip -d $GENOME.gz

# Create genome indexes
samtools faidx $GENOME
bwa index $GENOME

# Download RefSeq annotations
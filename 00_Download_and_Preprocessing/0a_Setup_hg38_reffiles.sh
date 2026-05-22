#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=24gb
#SBATCH -t 6:00:00
#SBATCH -A open
#SBATCH -o logs/0a_Setup_hg38_reffiles.log.out
#SBATCH -e logs/0a_Setup_hg38_reffiles.log.err

# Download human (hg38) ref genome with annotations and build alignment indexes

set -exo
source activate bx

GENOME=../data/hg38_files/hg38.fa

[ -d ../data/hg38_files ] || mkdir -p ../data/hg38_files

# Download genome FASTA
wget -O $GENOME.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gzip -d $GENOME.gz

# Download hg38 Blacklist from ENCODE
wget -O ../data/hg38_files/ENCFF356LFX_hg38_exclude.bed.gz https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz

# Reslice into a BED6 format
gzip -dc ../data/hg38_files/ENCFF356LFX_hg38_exclude.bed.gz | awk 'BEGIN{OFS="\t";FS="\t"}{print $0,".","0","."}' > ../data/hg38_files/ENCFF356LFX_hg38_exclude.bed

# Create genome indexes
samtools faidx $GENOME
bowtie2-build $GENOME $GENOME
bowtie-build -C $GENOME $GENOME.colorspace
bwa index $GENOME

# Download RefSeq annotations
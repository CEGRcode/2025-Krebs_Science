#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=24gb
#SBATCH -t 6:00:00
#SBATCH -A open
#SBATCH -o logs/0_Setup_hg38_reffiles.log.out
#SBATCH -e logs/0_Setup_hg38_reffiles.log.err

module load anaconda
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/align

GENOME=../data/mm10_files/mm108.fa

[ -d ../data/mm10_files ] || mkdir ../data/mm10_files

# Download genome FASTA
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm100/bigZips/mm10.fa.gz
mv hg38.fa.gz $GENOME.gz
gzip -d $GENOME.gz

# Create genome indexes
samtools faidx $GENOME
bowtie2-build $GENOME $GENOME
bowtie-build -C $GENOME $GENOME.colorspace
bwa index $GENOME

# Download RefSeq annotations
#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=48gb
#SBATCH -t 24:00:00
#SBATCH -A open
#SBATCH -o logs/DNase-FLASH_download-align-dedup-filter.log.out
#SBATCH -e logs/DNase-FLASH_download-align-dedup-filter.log.err

# Download DNase FLASH sequencing data (paired-end, Illumina) then
# align, dedup, and filter the reads.

# data/BAM
#   |--DNase-FLASH_SRR801881_hg38.bam
#   |--DNase-FLASH_SRR801880_hg38.bam

THREADS=4

# Dependencies
# - bwa
# - java
# - samtools
# - SRA toolkit (fasterq-dump)

set -exo
source activate bx

# Inputs and outputs
BAMDIR=../data/BAM
GENOME=../data/hg38_files/hg38.fa

# Script shortcuts
PICARD=../bin/picard.jar

[ -d FQ ] || mkdir FQ
[ -d BAM ] || mkdir BAM

# Download raw FASTQ file
for SRR in "SRR801881" "SRR801880";
do
	fasterq-dump -O FQ --split-files $SRR &
done
wait

# Align each sample
for SRR in "SRR801881" "SRR801880";
do
	bwa mem -t $THREADS -v 1 $GENOME FQ/${SRR}_1.fastq FQ/${SRR}_2.fastq | samtools sort -@ $THREADS > BAM/${SRR}_sorted.bam
done

# Mark and remove duplicates for technical replicates pairs and single (non-technical sequencing) replicates
for SRR in "SRR801881" "SRR801880";
do
	java -jar $PICARD MarkDuplicates -I BAM/${SRR}_sorted.bam -O BAM/${SRR}_dedup.bam -M BAM/${SRR}_dedup.metrics \
		--REMOVE_DUPLICATES 'false' --ASSUME_SORTED 'true' --DUPLICATE_SCORING_STRATEGY 'SUM_OF_BASE_QUALITIES' \
		--READ_NAME_REGEX '[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*.' --OPTICAL_DUPLICATE_PIXEL_DISTANCE '100'
		--VALIDATION_STRINGENCY 'LENIENT' --TAGGING_POLICY All --QUIET true --VERBOSITY ERROR &
done
wait

# Filter resulting BAM to match PEGR
for SRR in "SRR801881" "SRR801880";
do
	samtools view -h -b -f 0x1 -F 0x404 BAM/${SRR}_dedup.bam -o BAM/${SRR}_hg38.bam &
done
wait

# Index
for SRR in "SRR801881" "SRR801880";
do
	samtools index BAM/DNase-FLASH_${SRR}_hg38.bam 
done

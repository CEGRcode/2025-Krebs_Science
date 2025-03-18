#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=48gb
#SBATCH -t 24:00:00
#SBATCH -A open
#SBATCH -o logs/MNase-ChIP_download-align-dedup-filter-merge.log.out
#SBATCH -e logs/MNase-ChIP_download-align-dedup-filter-merge.log.err

# Download MNase ChIP sequencing data (paired-end, Illumina) that targets H3K4me3
# and align, dedup, filter, and merge the reads.

# data/BAM
#   |--MNase-ChIP_H3K4me3_merge_hg38.bam

### CHANGE ME
WRK=storage/group/bfp2/default/hxc585_HainingChen/2025_Chen_TF-Nuc/00_Download_and_Preprocessing
WRK=/storage/home/owl5022/scratch/2024-Krebs_Science/00_Download_and_Preprocessing
THREADS=4
###

# Dependencies
# - bwa
# - java
# - samtools
# - SRA toolkit (fasterq-dump)

set -exo
module load anaconda3
module load samtools
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/align

# Inputs and outputs
BAMDIR=$WRK/../data/BAM
GENOME=$WRK/../data/hg38_files/hg38.fa

# Script shortcuts
PICARD=$WRK/../bin/picard.jar

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
	samtools view -h -b -f 0x1 -F 0x404 BAM/${SRR}_dedup.bam -o BAM/$SRR.bam &
done
wait

# Index
for SRR in "SRR801881" "SRR801880";
do
	samtools index BAM/DNase-FLASH_${SRR}.bam 
done

#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 01:00:00
#SBATCH -A open
#SBATCH -o logs/1_align_data.log.out-%a
#SBATCH -e logs/1_align_data.log.err-%a
#SBATCH --array 1-NSAMPLES

# Script to align samples according to our Galaxy Workflow

### CHANGE ME
WRK=/path/to/2024-Chen_Nature/00_Download_and_Preprocessing
WRK=/storage/home/owl5022/scratch/2024-Chen_Nature/00_Download_and_Preprocessing
###

# Dependencies
# - BWA
# - CWpair
# - Genetrack
# - java
# - MEME
# - samtools

set -exo
module load bwa
module load samtools
module load anaconda
source activate bx

# Inputs and outputs
FASTQ1=`ls ../data/FASTQ/*_R1.fastq.gz | head -n $SLURM_ARRAY_TASK_ID | tail -1`
SAMPLE=`basename $FASTQ1 "R1.fastq.gz"`
FASTQ2=../data/FASTQ/$SAMPLE\_R2.fastq.gz
BAM=../data/sample-BAM/$SAMPLE
MEME=../data/sample-MEME/$SAMPLE
GENOME=../data/hg38_files/hg38.fa

# Script shortcuts
PICARD=../bin/picard.jar
SCRIPTMANAGER=../bin/ScriptManager-v0.15.jar

[ -d logs ] || mkdir logs

# Initial alignment
bwa mem -t 4 $GENOME $FASTQ1 $FASTQ2 | samtools view -Sb - | samtools sort -o $BAM.bam -
# Mark duplicates
java -jar -Xmx8G $PICARD MarkDuplicates \
    REMOVE_DUPLICATES='false' ASSUME_SORTED='true' \
    DUPLICATE_SCORING_STRATEGY='SUM_OF_BASE_QUALITIES' \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE='100' \
    VALIDATION_STRINGENCY='LENIENT' \
    TAGGING_POLICY=All QUIET=true VERBOSITY=ERROR \
    INPUT=$BAM.bam OUTPUT=$BAM\_MARKDUP.bam METRICS_FILE=$BAM\_duplication.out
# Remove duplicates
samtools view -f 0x1 -F 0x404 $BAM\_MARKDUP.bam > $BAM.bam
# Index results
samtools index $BAM.bam

## ADD MEME MOTIF CALLING HERE
# generate $MEME.meme.txt

# Clean-up
rm $BAM\_MARKDUP.bam

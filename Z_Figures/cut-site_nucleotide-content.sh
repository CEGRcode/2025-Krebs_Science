#!/bin/bash

# Perform sequence content analysis around cut sites

### CHANGE ME
WRK=/path/to/2024-Krebs_Science/Z_Figures
WRK=/storage/group/bfp2/default/owl5022-OliviaLang/2024-Krebs_Science/Z_Figures
###

# Dependencies
# - java
# - opencv
# - python

set -exo
module load anaconda
module load samtools
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

UPDOWN_KMER=../bin/updownstream_di-nt_tally.py
KMER2NT=../bin/dint_to_nt_positional_count_matrix.py
STACKNT=../bin/make_stack_barchart_TSV.py

BAMDIR=../data/BAM

# Index genome if doesnt exist
GENOME=../data/hg38_files/hg38.fa
[ -f $GENOME.fai ] || samtools faidx $GENOME

# Create directory for output
TEMP=UpstreamKmerAnalysis
[ -d $TEMP ] || mkdir $TEMP
[ -d $TEMP/BAM ] || mkdir $TEMP/BAM
[ -d $TEMP/DiTally ] || mkdir $TEMP/DiTally

# Get BAM info
# samtools flagstat $BAMDIR/MNase-seq_21U_rep1_hg38.bam         > $TEMP/BAM/MNase-seq_21U_rep1_hg38.bam.flagstat
# samtools flagstat $BAMDIR/BNase-seq_50U-10min_merge_hg38.bam  > $TEMP/BAM/BNase-seq_50U-10min_merge_hg38.bam.flagstat
# samtools flagstat $BAMDIR/DNase-seq_ENCFF425WDA_rep1_hg38.bam > $TEMP/BAM/DNase-seq_ENCFF425WDA_rep1_hg38.bam.flagstat
# samtools flagstat $BAMDIR/DNase-seq_ENCFF518XTC_rep1_hg38.bam > $TEMP/BAM/DNase-seq_ENCFF518XTC_rep1_hg38.bam.flagstat
# samtools flagstat $BAMDIR/MNase-seq_5U_rep1_hg38.bam          > MNase-seq_5U_rep1_hg38.bam.flagstat
# samtools flagstat $BAMDIR/MNase-seq_21U_rep1_hg38.bam         > MNase-seq_21U_rep1_hg38.bam.flagstat
# samtools flagstat $BAMDIR/MNase-seq_79U_rep1_hg38.bam         > MNase-seq_79U_rep1_hg38.bam.flagstat
# samtools flagstat $BAMDIR/MNase-seq_304U_rep1_hg38.bam        > $TEMP/BAM/MNase-seq_304U_rep1_hg38.bam.flagstat

# Subsample BAM - hardcoded this for different BAM files (Seed=2)
# samtools view -b -s 2.0001 $BAMDIR/MNase-seq_21U_rep1_hg38.bam        > $TEMP/BAM/MNase-seq_21U_rep1_hg38_SUBSAMPLE.bam # 170M --> 17K reads reads
# samtools view -b -s 2.002 $BAMDIR/DNase-seq_ENCFF425WDA_rep1_hg38.bam    > $TEMP/BAM/DNase-seq_ENCFF425WDA_rep1_hg38_SUBSAMPLE.bam  # 18M --> 35K reads
# samtools view -b -s 2.0002 $BAMDIR/DNase-seq_ENCFF518XTC_rep1_hg38.bam    > $TEMP/BAM/DNase-seq_ENCFF518XTC_rep1_hg38_SUBSAMPLE.bam  # 98M --> 19K reads
# samtools view -b --threaded 4 -s 2.00001 $BAMDIR/BNase-seq_50U-10min_merge_hg38.bam > $TEMP/BAM/BNase-seq_50U-10min_merge_hg38_SUBSAMPLE.bam # 1.3B --> 12K reads
# samtools view -b -s 2.001 $BAMDIR/BNase-seq_50U-10min_1_hg38.bam > $TEMP/BAM/BNase-seq_50U-10min_1_hg38_SUBSAMPLE.bam # 22M --> 22K reads
# samtools view -b -s 2.002 $BAMDIR/BNase-seq_50U-30min_merge_hg38.bam > $TEMP/BAM/BNase-seq_50U-30min_merge_hg38_SUBSAMPLE.bam # 13M --> 26K reads

# Add your BAM files here
for BAMFILE in "${BAMDIR}/BNase-seq_50U-10min_1_hg38.bam" "${BAMDIR}/MNase-seq_21U_rep1_hg38.bam" "${BAMDIR}/BNase-seq_50U-30min_merge_hg38.bam" "${BAMDIR}/DNase-seq_ENCFF425WDA_rep1_hg38.bam" "${BAMDIR}/DNase-seq_ENCFF518XTC_rep1_hg38.bam";
do
    BAM=`basename $BAMFILE ".bam"`

    # Subsample BAM to 1%
    # samtools view -b -s 0.0001 $BAMFILE > $TEMP/BAM/${BAM}_SUBSAMPLE.bam

    # Get sample BAM info
    samtools flagstat $TEMP/BAM/${BAM}_SUBSAMPLE.bam > $TEMP/BAM/${BAM}_SUBSAMPLE.bam.flagstat

    # ===Read 1===

    # Perform upstream kmer analysis around 5' cut sites (di, -50 to +100)
    python $UPDOWN_KMER -l 50 -r 100 --read1 -p \
         -i $TEMP/BAM/${BAM}_SUBSAMPLE.bam -g $GENOME \
         -o $TEMP/DiTally/${BAM}_SUBSAMPLE_DINT-l50r100-R1.tsv

    # Re-tally for single nucleotide counts
    python $KMER2NT -i $TEMP/DiTally/${BAM}_SUBSAMPLE_DINT-l50r100-R1.tsv -o $TEMP/DiTally/${BAM}_SUBSAMPLE_NT-l50r100-R1.tsv

    # Generate Figure: stack single nucleotides (enforce same-frequency)
    python $STACKNT --entropy -i <(cut -f1,40-70 $TEMP/DiTally/${BAM}_SUBSAMPLE_NT-l50r100-R1.tsv) --title ${BAM}_SUBSAMPLE -o $TEMP/DiTally/${BAM}_SUBSAMPLE_NT-l50r100-R1.svg

    # ===Read 2===

    # Perform upstream kmer analysis around 5' cut sites (di, -50 to +100)
    python $UPDOWN_KMER -l 50 -r 100 --read2 -p \
         -i $TEMP/BAM/${BAM}_SUBSAMPLE.bam -g $GENOME \
         -o $TEMP/DiTally/${BAM}_SUBSAMPLE_DINT-l50r100-R2.tsv

    # Re-tally for single nucleotide counts
    python $KMER2NT -i $TEMP/DiTally/${BAM}_SUBSAMPLE_DINT-l50r100-R2.tsv -o $TEMP/DiTally/${BAM}_SUBSAMPLE_NT-l50r100-R2.tsv

    # Generate Figure: stack single nucleotides (enforce same-frequency)
    python $STACKNT --entropy -i <(cut -f1,40-70 $TEMP/DiTally/${BAM}_SUBSAMPLE_NT-l50r100-R2.tsv) --title ${BAM}_SUBSAMPLE -o $TEMP/DiTally/${BAM}_SUBSAMPLE_NT-l50r100-R2.svg



done
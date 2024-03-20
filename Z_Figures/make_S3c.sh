#!/bin/bash

# Get set of associated gene products for functional analysis (S3b)

### CHANGE ME
WRK=/path/to/2024-Chen_Nature/Z_Figures
WRK=/storage/home/owl5022/scratch/2024-Chen_Nature/Z_Figures
###

# Dependencies
# - bedtools
# - java
# - samtools

set -exo
module load bedtools
module load samtools

# Fill in placeholder constants with your directories
MOTIF=$WRK/../data/RefPT-Motif
OTHER=$WRK/../data/RefPT-Other
ANNOTATIONS=../data/RefPT-Other/hg19.refGene.gtf 

# Setup ScriptManager for job array
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.14.jar

# Set up output directories
[ -d logs ] || mkdir logs
[ -d S3/b ] || mkdir -p S3/b

# Reformat to just keep CDS?
#awk '{OFS="\t"}{FS="\t"}{if ($3=="5UTR" || $3=="CDS" ||  $3=="exon") print}' $ANNOTATIONS > reformatAnnotations.gtf

# Convert GTF/GFF to BED format
#java -jar $SCRIPTMANAGER coordinate-manipulation gff-to-bed reformatAnnotations.gtf -o hg19_refGene.bed
java -jar $SCRIPTMANAGER coordinate-manipulation gff-to-bed $ANNOTATIONS -o $OTHER/hg19_refGene.bed

# Sort genomically
#bedtools sort -i hg19_refGene.bed > RefSeq.bed
bedtools sort -i $OTHER/hg19_refGene.bed > RefSeq.bed
bedtools sort -i $OTHER/TSS_DistWDR5.bed | awk '{OFS="\t"}{FS="\t"}{print $1,$2,$3,$1"_"$2"_"$3"_"$6,$5,$6}' > TSS.bed

# Get closest downstream RefSeq annotation gene id (no more than 500bp downstream)
bedtools closest -a TSS.bed -b RefSeq.bed -iu -d -D a -t first \
	| awk '{OFS="\t"}{FS="\t"}{if($13>-500 && $13<500) print}' \
	> Genes.bed
paste <(awk -F"\"" '{print $2}' Genes.bed) Genes.bed | sort -k1,1 > S3/b/Genes.txt

# Clean-up
rm Genes.bed TSS.bed RefSeq.bed
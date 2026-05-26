#!/bin/bash

# Download annotations and scores for hg38 genome build for later processing

set -exo
source activate bx

# Inputs and outputs
CHRSZ=../data/hg38_files/hg38.chrom.sizes
CDIR=../data/Conservation-SNP
REFPT=../data/RefPT-Other

# Script shortcuts
BBTOBED=../bin/bigBedToBed
BGTOBW=../bin/bedGraphToBigWig

[ -d $CDIR ] || mkdir $CDIR

# =====Conservation processing=====

# Download
wget -c -O $CDIR/hg38.phyloP30way.bw https://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP30way/hg38.phyloP30way.bw

# =====SNP processing=====

# Download BigBed
wget -c http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp153.bb

# Decompress to standard BED format
$BBTOBED dbSnp153.bb dbSnp153.bed

# Filter to keep snv (remove ins,del,delins,mnv) and remove non-standard chr
awk '{OFS="\t"}{FS="\t"}{ if ($14=="snv") print; }' dbSnp153.bed \
	| grep -v '_fix' | grep -v '_random' | grep -v 'chrUn_' | grep -v '_alt' \
	> dbSnp153_snv.bed

# Rebuild each filtered BED file as a BedGraph file
#   - strand agnostic count of chr-start-stop
#   - reformat counts into score column of bedgraph
cut -f1-3 dbSnp153_snv.bed | sort | uniq -c \
	| awk '{OFS="\t"}{print $2,$3,$4,$1}' > dbSnp153_snv.bedGraph

# Convert BedGraph to BigWig (for BigWig-style pileup)
$BGTOBW dbSnp153_snv.bedGraph $CHRSZ $CDIR/dbSnp153_snv.bw

# Clean-up
rm dbSnp153.bb dbSnp153.bed dbSnp153_snv.bed dbSnp153_snv.bedGraph


# =====CpG Islands=====

# Download ../data/RefPT-Other/CpGIslands.tsv.gz
#  1. Navigate to https://genome.ucsc.edu/cgi-bin/hgTables
#  2. Select dataset
#    - Clade: "Mammal"
#    - Genome: "Human"
#    - Assembly: "Dec. 2013 (GRCh38/hg38)"
#    - Group: "Regulation"
#    - Track: "CpG Islands"
#    - Table: "cpgIslandExt"
#  3. Define region of interest
#    - Region: "Genome"
#  4. Click "Get output"
#  5. Save file to ../data/RefPT-Other/CpGIslands.tsv

# Download BigWig for more dynamic signal
# wget -c -O ../data/RefPT-Other/gc5Base.bw https://hgdownload.soe.ucsc.edu/gbdb/hg38/bbi/gc5BaseBw/gc5Base.bw

# Reformat columns into BED6 format
gzip -dc $REFPT/CpGIslands.tsv.gz | sed '1d' | awk 'BEGIN{OFS="\t";FS="\t"}{print $2,$3,$4,$5,$6,"."}' > $REFPT/CpG_Islands.bed



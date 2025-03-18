#!/bin/bash

### CHANGE ME
WRK=/storage/group/bfp2/default/hxc585_HainingChen/2025_Chen_TF-Nuc/00_Download_and_Preprocessing
WRK=/storage/home/owl5022/scratch/2024-Krebs_Science/00_Download_and_Preprocessing
###

# Inputs and outputs
CDIR=$WRK/../data/Conservation-SNP
CHRSZ=$WRK/../data/hg38_files/hg38.chrom.sizes

# Script shortcuts
BBTOBED=$WRK/../bin/bigBedToBed
BGTOBW=$WRK/../bin/bedGraphToBigWig

[ -d $CDIR ] || mkdir $CDIR
cd $CDIR

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
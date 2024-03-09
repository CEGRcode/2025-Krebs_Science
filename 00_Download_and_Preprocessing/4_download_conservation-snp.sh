#!/bin/bash

BWTOBG=../bin/bigWigToBedGraph
BBTOBED=../bin/bigBedToBed
BGTOBW=../bin/bedGraphToBigWig
SPLITALT=../bin/split_snps_alt.py
SPLITTYPE=../bin/split_snps_type.py

CDIR=../data/Conservation-SNP
CHRSZ=../data/hg19_files/hg19.chrom.sizes

# Conservation processing - Download
wget -c https://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals.phyloP46way.bw
wget -c https://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates.phyloP46way.bw
wget -c https://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate.phyloP46way.bw

mv *phylocP46way.bw $CDIR/

# $BWTOBG placentalMammals.phyloP46way.bw placentalMammals.phyloP46way.bedgraph
# $BWTOBG primates.phyloP46way.bw primates.phyloP46way.bedgraph
# $BWTOBG vertebrate.phyloP46way.bw vertebrate.phyloP46way.bedgraph

# SNP processing - download
wget -c http://hgdownload.soe.ucsc.edu/gbdb/hg19/snp/dbSnp153.bb
$BBTOBED dbSnp153.bb dbSnp153.bed

# SNP processsing - split by SNV or INDELS
awk '{OFS="\t"}{FS="\t"}{
 if ($14=="snv") {
   print $0 > "dbSnp153_snv.bed";
 } else if ($14=="delins") {
   print $0 > "dbSnp153_delins.bed";
 } else if ($14=="del") {
   print $0 > "dbSnp153_del.bed";
 } else if ($14=="ins") {
   print $0 > "dbSnp153_ins.bed";
 } else if ($14=="mnv") {
   print $0 > "dbSnp153_mnv.bed";
 }
}' dbSnp153.bed

# SNP processsing - split by SNV or INDELS (same as above but one-by-one version. Good if you are only interested in one of the types and limited on storage)
#awk '{OFS="\t"}{FS="\t"}{if($14=="snv") print}' dbSnp153.bed > dbSnp153_snv.bed
#awk '{OFS="\t"}{FS="\t"}{if($14=="del") print}' dbSnp153.bed > dbSnp153_del.bed
#awk '{OFS="\t"}{FS="\t"}{if($14=="ins") print}' dbSnp153.bed > dbSnp153_ins.bed
#awk '{OFS="\t"}{FS="\t"}{if($14=="delins") print}' dbSnp153.bed > dbSnp153_delins.bed

# SNP processing - split by ref nucleotide
awk '{FS="\t"}{OFS="\t"}{
 if ($5=="A") {
   print $0 > "dbSnp153_snv_ref-A.bed";
 } else if ($5=="T") {
   print $0 > "dbSnp153_snv_ref-T.bed";
 } else if ($5=="C") {
   print $0 > "dbSnp153_snv_ref-C.bed";
 } else if ($5=="G") {
   print $0 > "dbSnp153_snv_ref-G.bed";
 }
}' dbSnp153_snv.bed

# SNP processing - split by alt nucleotide
python $SPLITALT -i dbSnp153_snv.bed

# SNP processing - split by variant annotations
python $SPLITTYPE -i dbSnp153_snv.bed


[ -d process-snp ] || mkdir process-snp

# Rebuild each filtered BED file as a BigWig file (for BigWig-style pileup)
for FILE in dbSnp153_snv_ref-*.bed;
do
	BASE=`basename $FILE ".bed"`
	echo $BASE

	# Make a smaller BED
	bedtools intersect -sorted -a $FILE -b <(bedtools sorted -i )

	# Assume all notated for "+" strand
	awk '{FS="\t"}{OFS="\t"}{print $1,$2,$3}' $FILE \
		| grep -v '_fix' |grep -v '_random' | grep -v 'chrUn_' \
		| grep -v '_alt' | uniq -c | awk '{OFS="\t"}{print $2,$3,$4,$1}' > process-snp/$BASE.bedGraph
	$BGTOBW process-snp/$BASE.bedGraph $CHRSZ process-snp/$BASE.bw
done

# Reorganize
mv process-snp/dbSnp153*.bw $CDIR/
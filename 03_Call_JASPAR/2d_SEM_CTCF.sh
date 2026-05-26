#!/bin/bash

# data/RefPT-JASPAR-nonK562
#   |--<TFNAME>_<JASPARID>_Unbound.bed
#   |--<TFNAME>_<JASPARID>_K562-specific-Unbound.bed

### CHANGE ME
WRK=/storage/group/bfp2/default/hxc585_HainingChen/2025_Chen_TF-Nuc/03_Call_JASPAR
#WRK=/ocean/projects/see180003p/owlang/2024-Krebs_Science/03_Call_JASPAR
#WRK=/scratch/owl5022/2024-Krebs_Science/03_Call_JASPAR

# - java

#set -exo
module load bedtools
module load anaconda3
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Inputs and outputs
GENOME=$WRK/../data/hg38_files/hg38.fa
BLACKLIST=$WRK/../data/hg38_files/ENCFF356LFX_hg38_exclude.bed
WT_SEM_CTCFBAMFILE=$WRK/../data/BAM/SEM_CTCF_BX_WT_hg38.bam
CTCFKD_SEM_CTCFBAMFILE=$WRK/../data/BAM/SEM_CTCF_BX_CTCTKO_hg38.bam
CTCFBEDFOLDER=$WRK/../data/RefPT-JASPAR/
OUTDIR=$WRK/../data/RefPT-JASPAR-nonK562/


# Script shortcuts
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar

mkdir -p $OUTDIR
mkdir -p $OUTDIR/1000bp
mkdir -p $OUTDIR/150bp
## sort CTCF_MA1929.1.bed by WT_SEM_CTCF, and sort by ratio of CTCFKD_SEM_CTCFBAMFILE to WT_SEM_CTCF
cat $CTCFBEDFOLDER/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile*.bed | bedtools sort -i | uniq > $OUTDIR/CTCF_MA1929.1.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $OUTDIR/CTCF_MA1929.1.bed -o $OUTDIR/CTCF_MA1929.1_100bp.bed
java -jar "$SCRIPTMANAGER" read-analysis tag-pileup  $OUTDIR/CTCF_MA1929.1_100bp.bed $WT_SEM_CTCFBAMFILE -1 -s 6  --combined --cpu 4 -M  $OUTDIR/WT_SEM_CTCF_CTCF_MA1929.1_100bp_read1
java -jar "$SCRIPTMANAGER" read-analysis aggregate-data --sum $OUTDIR/WT_SEM_CTCF_CTCF_MA1929.1_100bp_read1_combined.cdt -o $OUTDIR/WT_SEM_CTCF_CTCF_MA1929.1_100bp_read1_combined.out
java -jar "$SCRIPTMANAGER" read-analysis tag-pileup  $OUTDIR/CTCF_MA1929.1_100bp.bed $CTCFKD_SEM_CTCFBAMFILE -1 -s 6 --combined --cpu 4 -M  $OUTDIR/CTCTKD_SEM_CTCF_CTCF_MA1929.1_100bp_read1
java -jar "$SCRIPTMANAGER" read-analysis aggregate-data --sum $OUTDIR/CTCTKD_SEM_CTCF_CTCF_MA1929.1_100bp_read1_combined.cdt -o $OUTDIR/CTCTKD_SEM_CTCF_CTCF_MA1929.1_100bp_read1_combined.out


cut -f 2  $OUTDIR/WT_SEM_CTCF_CTCF_MA1929.1_100bp_read1_combined.out | tail -n +2  | cut -f 2 | \
paste $OUTDIR/CTCF_MA1929.1.bed - > $OUTDIR/CTCF_MA1929.1_WT-SEM-CTCF.tsv
cut -f 2  $OUTDIR/CTCTKD_SEM_CTCF_CTCF_MA1929.1_100bp_read1_combined.out | tail -n +2  | cut -f 2 | \
paste $OUTDIR/CTCF_MA1929.1_WT-SEM-CTCF.tsv - > $OUTDIR/CTCF_MA1929.1_WT-SEM-CTCF_CTCTKD_SEM_CTCF.tsv

rm  $OUTDIR/WT_SEM_CTCF_CTCF_MA1929.1_100bp_read1_combined.cdt $OUTDIR/CTCTKD_SEM_CTCF_CTCF_MA1929.1_100bp_read1_combined.cdt
rm   $OUTDIR/CTCTKD_SEM_CTCF_CTCF_MA1929.1_100bp_read1_combined.out $OUTDIR/WT_SEM_CTCF_CTCF_MA1929.1_100bp_read1_combined.out
rm $OUTDIR/CTCF_MA1929.1_100bp.bed
rm  $OUTDIR/CTCF_MA1929.1_WT-SEM-CTCF.tsv


bedtools sort -i $OUTDIR/CTCF_MA1929.1_WT-SEM-CTCF_CTCTKD_SEM_CTCF.tsv  | awk -v DIR="$OUTDIR" 'BEGIN{OFS="\t";FS="\t"}{
    if ($7 >= 1) {
        system("echo \"" $0 "\" >> " DIR "/CTCF_MA1929.1_SORT-WT-SEM-CTCF.tsv")
    }
}'  

awk  'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$8/$7}' $OUTDIR/CTCF_MA1929.1_SORT-WT-SEM-CTCF.tsv | sort -k9,9nr > $OUTDIR/CTCF_MA1929.1_SORT-CTCFKD-to-WT-SEM-CTCF.tsv

rm $OUTDIR/CTCF_MA1929.1_SORT-WT-SEM-CTCF.tsv

wc -l $OUTDIR/CTCF_MA1929.1_SORT-CTCFKD-to-WT-SEM-CTCF.tsv
#  16671 /CTCF_MA1929.1_SORT-CTCFKD-to-WT-SEM-CTCF.tsv

sort -k7,7nr $OUTDIR/CTCF_MA1929.1_SORT-CTCFKD-to-WT-SEM-CTCF.tsv  | head -n 4168 | cut -f 1-6 - > $OUTDIR/CTCF_MA1929.1_1.bed

for file in  $OUTDIR/CTCF_MA1929.1_1.bed  ; do
	filename=`basename $file ".bed"`
	java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $file -o $OUTDIR/1000bp/${filename}_1000bp.bed
	java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 150 $file -o $OUTDIR/150bp/${filename}_150bp.bed
done

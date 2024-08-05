#!/bin/bash

# get 
module load anaconda
source activate bioinfo
### CHANGE ME
WRK=/storage/group/bfp2/default/hxc585_HainingChen/
MOTIF=$WRK/Fox_NFIA_CTCF/data/RefPT-Motif
SCRIPTMANAGER=$WRK/2023_Chen_PIC3/bin/ScriptManager-v0.14.jar
GENOME=$WRK/2023_Chen_PIC3/hg19_files/hg19.fa
Genome=$WRK/2023_Chen_PIC3/hg19_files/hg19.info
Nuc="$WRK/2024-Chen_Benzonase-ChIP-exo/03_RfMotif_sort_Nuc/K562_nuc_uHex_uTetra_1bp_sort_nostrand.bed"
NFIA1=$WRK/Fox_NFIA_CTCF/data/BAM/K562_NFIA_BX_rep1_hg19.bam
NFIA2=$WRK/Fox_NFIA_CTCF/data/BAM/K562_NFIA_BX_rep2_hg19.bam
P3001=$WRK/Fox_NFIA_CTCF/data/BAM/K562_EP300_BX_rep1_hg19.bam
K562IgG=$WRK/Fox_NFIA_CTCF/data/BAM/K562_IgG_BX_merge.bam
Input=$WRK/Fox_NFIA_CTCF/data/BAM/K562_-_BI_rep1_hg19.bam
H3=$WRK/2023_Chen_PIC3/data/his_BAM/K562_H3_BX_rep1_hg19.bam
H3K27ac=$WRK/2023_Chen_PIC3/data/his_BAM/K562_H3K27ac_BX_rep1_hg19.bam
H2AZ=$WRK/2023_Chen_PIC3/data/his_BAM/K562_H2A.Z_BX_rep1_hg19.bam
H3K4me3=$WRK/2023_Chen_PIC3/data/his_BAM/K562_H3K4me3_BX_rep1_hg19.bam
H3K9ac=$WRK/2023_Chen_PIC3/data/his_BAM/K562_H3K9ac_BX_rep1_hg19.bam
CTCF=$WRK/2024-Chen_Benzonase-ChIP-exo/data/RefPT-Motif/CTCF_Occupancy_1bp.bed
GABPA=$WRK/2024-Chen_Benzonase-ChIP-exo/data/RefPT-Motif/GABPA_Occupancy.bed
SP1=$WRK/2024-Chen_Benzonase-ChIP-exo/data/RefPT-Motif/SP1_Occupancy_1bp.bed
NRF1=$WRK/2024-Chen_Benzonase-ChIP-exo/data/RefPT-Motif/NRF1_Occupancy.bed
USF2=$WRK/2024-Chen_Benzonase-ChIP-exo/data/RefPT-Motif/USF2_Occupancy_1bp.bed
YY1=$WRK/2024-Chen_Benzonase-ChIP-exo/data/RefPT-Motif/YY1_Occupancy_1bp.bed
CoPROBAMFILE=$WRK/2023_Chen_PIC3/data/RNA_seq/K562_Procap_PE_ENCLB382OOO.bam
# determine the nucleosome engagment direction on each NFIA sites
wc -l $MOTIF/1bp/NFIA_Occupancy_1bp.bed
6225 /storage/group/bfp2/default/hxc585_HainingChen/Fox_NFIA_CTCF/data/RefPT-Motif/1bp/NFIA_Occupancy_1bp.bed
bedtools shift -i $MOTIF/1bp/NFIA_Occupancy_1bp.bed -g $Genome -p 120 -m -120 > $MOTIF/1bp/NFIA_Occupancy_down120bp.bed
bedtools shift -i $MOTIF/1bp/NFIA_Occupancy_1bp.bed -g $Genome -p -120 -m 120 > $MOTIF/1bp/NFIA_Occupancy_up120bp.bed

for file in $MOTIF/1bp/NFIA_Occupancy_down120bp.bed $MOTIF/1bp/NFIA_Occupancy_up120bp.bed  ; do
  filename=$(basename "$file" ".bed")
  java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 80 $file -o $MOTIF/1bp/${filename}_80bp.bed
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup $MOTIF/1bp/${filename}_80bp.bed "$NFIA1" -n 100 --cpu 4 -M $MOTIF/1bp/NFIA_${filename}_80bp_read1_MIN100
  rm $MOTIF/1bp/${filename}_80bp.bed
done

java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum $MOTIF/1bp/NFIA_NFIA_Occupancy_down120bp_80bp_read1_MIN100_anti.cdt -o $MOTIF/1bp/
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum $MOTIF/1bp/NFIA_NFIA_Occupancy_up120bp_80bp_read1_MIN100_sense.cdt -o $MOTIF/1bp/

cut -f 2 "$MOTIF/1bp/NFIA_NFIA_Occupancy_down120bp_80bp_read1_MIN100_anti_SCORES.out" | \
paste "$MOTIF/1bp/NFIA_NFIA_Occupancy_up120bp_80bp_read1_MIN100_sense_SCORES.out" - | \
tail -n +2 | cut -f 2-3 | \
paste "$MOTIF/1bp/NFIA_Occupancy_1bp.bed" - | \
awk -v motif="$MOTIF" '{ if (($7 > $8) || (($7 == $8) && ($6 == "+"))) print $0 > motif "/1bp/NFIA_Nuc_upstream_temp.bed"; else print $0 > motif "/1bp/NFIA_Nuc_downstream_temp.bed"   }'
wc -l $MOTIF/1bp/NFIA_Nuc_upstream_temp.bed
3307 /storage/group/bfp2/default/hxc585_HainingChen//Fox_NFIA_CTCF/data/RefPT-Motif/1bp/NFIA_Nuc_upstream_temp.bed
wc -l $MOTIF/1bp/NFIA_Nuc_downstream_temp.bed
2918 /storage/group/bfp2/default/hxc585_HainingChen//Fox_NFIA_CTCF/data/RefPT-Motif/1bp/NFIA_Nuc_downstream_temp.bed

rm NFIA/*_temp_*.out

awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,"engageNuddown"}' $MOTIF/1bp/NFIA_Nuc_downstream_temp.bed  > $MOTIF/1bp/NFIA_engagedNuc_downstream.bed
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,"engageNudup"}' $MOTIF/1bp/NFIA_Nuc_upstream_temp.bed > $MOTIF/1bp/NFIA_engagedNuc_upstream.bed
rm $MOTIF/1bp/NFIA_Nuc_downstream_temp.bed $MOTIF/1bp/NFIA_Nuc_upstream_temp.bed 
rm $MOTIF/1bp/NFIA_NFIA_Occupancy_down120bp_80bp_read1_*.cdt $MOTIF/1bp/NFIA_NFIA_Occupancy_up120bp_80bp_read1_*_SCORES.out
rm $MOTIF/1bp/NFIA_Occupancy_down120bp.bed $MOTIF/1bp/NFIA_Occupancy_up120bp.bed

## convert NFIA all to engage Nuc downstream

awk -v motif="$MOTIF" '{ if ($6 == "+") print $0 > motif "/1bp/NFIA_Nuc_upstream_+.bed"; else print $0 > motif "/1bp/NFIA_Nuc_upstream_-.bed" }'  $MOTIF/1bp/NFIA_engagedNuc_upstream.bed
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,"+",$8,$7,$9}' $MOTIF/1bp/NFIA_Nuc_upstream_-.bed > $MOTIF/1bp/NFIA_Nuc_upstream_-_convert.bed
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,"-",$8,$7,$9}' $MOTIF/1bp/NFIA_Nuc_upstream_+.bed > $MOTIF/1bp/NFIA_Nuc_upstream_+_convert.bed

rm $MOTIF/1bp/NFIA_Nuc_upstream_-.bed $MOTIF/1bp/NFIA_Nuc_upstream_+.bed

cat $MOTIF/1bp/NFIA_Nuc_upstream_-_convert.bed $MOTIF/1bp/NFIA_Nuc_upstream_+_convert.bed $MOTIF/1bp/NFIA_engagedNuc_downstream.bed | \
sort -k7,7nr > $MOTIF/1bp/NFIA_engagedNuc_downstream_upNuc_sort.bed

rm $MOTIF/1bp/NFIA_Nuc_upstream_-_convert.bed $MOTIF/1bp/NFIA_Nuc_upstream_+_convert.bed

wc -l $MOTIF/1bp/NFIA_engagedNuc_downstream_upNuc_sort.bed
6225 /storage/group/bfp2/default/hxc585_HainingChen//Fox_NFIA_CTCF/data/RefPT-Motif/1bp/NFIA_engagedNuc_downstream_upNuc_sort.bed

head -n 1000 $MOTIF/1bp/NFIA_engagedNuc_downstream_upNuc_sort.bed > $MOTIF/1bp/NFIA_engagedNuc_downstream_lessorientated.bed
tail -n +3001 $MOTIF/1bp/NFIA_engagedNuc_downstream_upNuc_sort.bed | head -n 1000 - > $MOTIF/1bp/NFIA_engagedNuc_downstream_moreorientated.bed

bedtools sort -i  $MOTIF/1bp/NFIA_engagedNuc_downstream_upNuc_sort.bed | uniq | bedtools closest -a - -b $Nuc -d -D a -t first | sort -k16,16n  > $MOTIF/1bp/NFIA_engageNucdown_NucSort.bed 

awk -v motif="$MOTIF" '{
   if ($16 < -73 ) {
        print $0 > motif "/1bp/NFIA_NucSort-UPSTREAM.bed"
    } else if ($16 > 73 ) {
        print $0 > motif "/1bp/NFIA_NucSort-DOWNSTREAM.bed"
    } else {
      print $0 > motif "/1bp/NFIA_NucSort-OVERLAP.bed"
    }
}' $MOTIF/1bp/NFIA_engageNucdown_NucSort.bed  


wc -l $MOTIF/1bp/NFIA_NucSort_*.bed
1259 /storage/group/bfp2/default/hxc585_HainingChen//Fox_NFIA_CTCF/data/RefPT-Motif/1bp/NFIA_NucSort-DOWNSTREAM.bed
4064 /storage/group/bfp2/default/hxc585_HainingChen//Fox_NFIA_CTCF/data/RefPT-Motif/1bp/NFIA_NucSort-OVERLAP.bed
902 /storage/group/bfp2/default/hxc585_HainingChen//Fox_NFIA_CTCF/data/RefPT-Motif/1bp/NFIA_NucSort-UPSTREAM.bed

bedtools sort -i $MOTIF/1bp/scrambledNFIA_top10k_1bp.bed  | uniq | bedtools closest -a - -b $Nuc -d -D a | sort -k13,13n  > $MOTIF/1bp/scrambledNFIA_NucSort_1bp.bed
exit

mkdir -p NFIA

#java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/1bp/NFIA_Occupancy_1bp.bed -o $MOTIF/1000bp/NFIA_Occupancy_1000bp.bed

java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 500 $MOTIF/1bp/NFIA_engagedNuc_downstream_upNuc_sort.bed -o $MOTIF/1000bp/NFIA_engagedNuc_downstream_upNuc_sort_500bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/1bp/NFIA_engageNucdown_NucSort.bed  -o $MOTIF/1000bp/NFIA_engageNucdown_NucSort_1000bp.bed

java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/1bp/NFIA_NucSort-UPSTREAM.bed -o $MOTIF/1000bp/NFIA_NucSort-UPSTREAM_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/1bp/NFIA_NucSort-DOWNSTREAM.bed -o $MOTIF/1000bp/NFIA_NucSort-DOWNSTREAM_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/1bp/NFIA_NucSort-OVERLAP.bed -o $MOTIF/1000bp/NFIA_NucSort-OVERLAP_1000bp.bed

java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/1bp/NFIA_engagedNuc_downstream_lessorientated.bed -o $MOTIF/1000bp/NFIA_engagedNuc_downstream_lessorientated_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/1bp/NFIA_engagedNuc_downstream_moreorientated.bed -o $MOTIF/1000bp/NFIA_engagedNuc_downstream_moreorientated_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/1bp/scrambledNFIA_NucSort_1bp.bed -o $MOTIF/1000bp/scrambledNFIA_NucSort_1000bp.bed

for file in $MOTIF/1000bp/NFIA_engageNucdown_Nucsort_1000bp.bed $MOTIF/1000bp/scrambledNFIA_NucSort_1000bp.bed ; do
  BED=$(basename "$file" ".bed")
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup $file "$NFIA1" --cpu 4 -M NFIA/NFIA1_${BED}_read1
  java -jar "$SCRIPTMANAGER" figure-generation heatmap -p .95 --blue NFIA/NFIA1_${BED}_read1_sense.cdt -o NFIA/NFIA1_${BED}_read1_sense_treeview.png
  java -jar "$SCRIPTMANAGER" figure-generation heatmap -p .95 --red NFIA/NFIA1_${BED}_read1_anti.cdt -o NFIA/NFIA1_${BED}_read1_anti_treeview.png
  java -jar "$SCRIPTMANAGER" figure-generation merge-heatmap NFIA/NFIA1_${BED}_read1_sense_treeview.png NFIA/NFIA1_${BED}_read1_anti_treeview.png -o NFIA/NFIA1_${BED}_read1_merge_treeview.png
  java -jar $SCRIPTMANAGER figure-generation label-heatmap NFIA/NFIA1_${BED}_read1_merge_treeview.png -f 20 -l -500 -m 0 -r 500 -o NFIA/NFIA1_${BED}_read1_merge_treeview.svg
  rm NFIA/NFIA1_${BED}_read1_*.cdt NFIA/NFIA1_${BED}_read1_*.png
done

for file in $MOTIF/1000bp/NFIA_engagedNuc_downstream_upNuc_sort_500bp.bed ; do
  BED=$(basename "$file" ".bed")
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup $file "$NFIA1" --cpu 4 -M NFIA/NFIA1_${BED}_read1
  java -jar "$SCRIPTMANAGER" figure-generation heatmap -p .95 --blue NFIA/NFIA1_${BED}_read1_sense.cdt -o NFIA/NFIA1_${BED}_read1_sense_treeview.png
  java -jar "$SCRIPTMANAGER" figure-generation heatmap -p .95 --red NFIA/NFIA1_${BED}_read1_anti.cdt -o NFIA/NFIA1_${BED}_read1_anti_treeview.png
  java -jar "$SCRIPTMANAGER" figure-generation merge-heatmap NFIA/NFIA1_${BED}_read1_sense_treeview.png NFIA/NFIA1_${BED}_read1_anti_treeview.png -o NFIA/NFIA1_${BED}_read1_merge_treeview.png
  java -jar $SCRIPTMANAGER figure-generation label-heatmap NFIA/NFIA1_${BED}_read1_merge_treeview.png -f 20 -l -250 -m 0 -r 250 -o NFIA/NFIA1_${BED}_read1_merge_treeview.svg
  rm NFIA/NFIA1_${BED}_read1_*.cdt NFIA/NFIA1_${BED}_read1_*.png
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup $file "$K562IgG" --cpu 4 -M NFIA/IgG_${BED}_read1
  java -jar "$SCRIPTMANAGER" figure-generation heatmap -p .95 --blue NFIA/IgG_${BED}_read1_sense.cdt -o NFIA/IgG_${BED}_read1_sense_treeview.png
  java -jar "$SCRIPTMANAGER" figure-generation heatmap -p .95 --red NFIA/IgG_${BED}_read1_anti.cdt -o NFIA/IgG_${BED}_read1_anti_treeview.png
  java -jar "$SCRIPTMANAGER" figure-generation merge-heatmap NFIA/IgG_${BED}_read1_sense_treeview.png NFIA/IgG_${BED}_read1_anti_treeview.png -o NFIA/IgG_${BED}_read1_merge_treeview.png
  java -jar $SCRIPTMANAGER figure-generation label-heatmap NFIA/IgG_${BED}_read1_merge_treeview.png -f 20 -l -250 -m 0 -r -250 -o NFIA/IgG_${BED}_read1_merge_treeview.svg
  rm NFIA/IgG_${BED}_read1_*.cdt NFIA/IgG_${BED}_read1_*.png
done

for file in  $MOTIF/1000bp/NFIA_engagedNuc_downstream_moreorientated_1000bp.bed $MOTIF/1000bp/NFIA_engagedNuc_downstream_lessorientated_1000bp.bed  ; do
  filename=$(basename "$file" ".bed")
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup "$file" $NFIA1  --cpu 4 -2 -n 100 -o NFIA/NFIA1_${filename}_MIN100_read2.out
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup  "$file" $NFIA1 --cpu 4 -1 -n 100 -o NFIA/NFIA1_${filename}_MIN100_read1.out
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup  "$file" $K562IgG --cpu 4 -1 -n 100 -o NFIA/IgG_${filename}_MIN100_read1.out
done

for file in $Input $H3K27ac  ; do 
  BAM=$(basename "$file" ".bam")
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup $MOTIF/1000bp/NFIA_engageNucdown_Nucsort_1000bp.bed "$file" --cpu 4 -m -M NFIA/${BAM}_NFIA_engageNucdown_Nucsort_1000bp_midpoint
  java -jar "$SCRIPTMANAGER" figure-generation heatmap -p .95 --black NFIA/${BAM}_NFIA_engageNucdown_Nucsort_1000bp_midpoint_combined.cdt -o NFIA/${BAM}_NFIA_engageNucdown_Nucsort_1000bp_midpoint_combined_treeview.png
  java -jar $SCRIPTMANAGER figure-generation label-heatmap NFIA/${BAM}_NFIA_engageNucdown_Nucsort_1000bp_midpoint_combined_treeview.png -f 20 -l -500 -m 0 -r 500 -o NFIA/${BAM}_NFIA_engageNucdown_Nucsort_1000bp_midpoint_combined_treeview.svg
  rm NFIA/${BAM}_NFIA_*.cdt  NFIA/${BAM}_NFIA_*.png
done


for file in $MOTIF/1000bp/NFIA_engageNucdown_Nuc*_1000bp.bed  ; do
  filename=$(basename "$file" ".bed")
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup "$file" $NFIA1  --cpu 4 -2 -n 100 -o NFIA/NFIA1_${filename}_MIN100_read2.out
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup  "$file" $NFIA1 --cpu 4 -1 -n 100 -o NFIA/NFIA1_${filename}_MIN100_read1.out
done



for file in $Input ; do
  BAM=$(basename "$file" ".bam")
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup $MOTIF/1000bp/scrambledNFIA_NucSort_1000bp.bed "$file" --cpu 4 -m -M NFIA/${BAM}_scrambledNFIA_NucSort_1000bp_midpoint
  java -jar "$SCRIPTMANAGER" figure-generation heatmap -p .95 --black NFIA/${BAM}_scrambledNFIA_NucSort_1000bp_midpoint_combined.cdt -o NFIA/${BAM}_scrambledNFIA_NucSort_1000bp_midpoint_combined_treeview.png
  java -jar $SCRIPTMANAGER figure-generation label-heatmap NFIA/${BAM}_scrambledNFIA_NucSort_1000bp_midpoint_combined_treeview.png -f 20 -l -500 -m 0 -r 500 -o NFIA/${BAM}_scrambledNFIA_NucSort_1000bp_midpoint_combined_treeview.svg
  rm NFIA/${BAM}_scrambledNFIA_*.cdt  NFIA/${BAM}_scrambledNFIA_*.png
done

for file in $Input ; do
  BAM=$(basename "$file" ".bam")
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup $MOTIF/1000bp/scrambledNFIA_NucSort_1000bp.bed "$file" --cpu 4 -1 -o NFIA/${BAM}_scrambledNFIA_NucSort_1000bp_read1.out
  java -jar "$SCRIPTMANAGER" read-analysis tag-pileup $MOTIF/1000bp/scrambledNFIA_NucSort_1000bp.bed "$file" --cpu 4 -2 -o NFIA/${BAM}_scrambledNFIA_NucSort_1000bp_read2.out
done




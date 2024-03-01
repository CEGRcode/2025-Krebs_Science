Non-Motif reference points

## hg19_K562_Nucleosome.bed
Nucleosome maps from Krebs et al, (in prep)

## hg19_K562_Nucleosome_1bp_SORT-genomic_nostrand.bed.gz
The reference nucleosome map used by motif sorts (`NucSort`)
```
# Rexpand to center bp
java -jar ../../bin/ScriptManager-v0.14.jar coordinate-manipulation expand-bed -c 1 K562_nuc_uHex_uTetra.bed.gz -o 1bp/hg19_K562_Nucleosome_1bp.bed
# Sort by genomic coord (for BEDTools), remove strand info, set id to genomic coord, and compress
bedtools sort -i 1bp/hg19_K562_Nucleosome_1bp.bed \
    | awk '{OFS="\t"}{FS="\t"}{print $1,$2,$3,$1"_"$2"_"$3,$5,"."}' \
    | gzip > 1bp/hg19_K562_Nucleosome_1bp_SORT-genomic_nostrand.bed.gz
# Clean-up
rm 1bp/hg19_K562_Nucleosome_1bp.bed
```


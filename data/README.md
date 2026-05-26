# data

This directory stores large data (FASTQ & BAM files) and globally used data and reference files.

<details>
<summary> Full execution summary
</summary>

```
data
  |--RefPT-Krebs
    |--MinusOneDyad_SORT-DistToExpressedTSS.bed
    |--MinusOneDyad_SORT-DistToUnexpressedTSS.bed
    |--MinusOneDyad_SORT-Expression.bed
    |--NFR_SORT-NFRLength.bed
    |--PlusOneDyad_SORT-DistToExpressedTSS.bed
    |--PlusOneDyad_SORT-DistToUnexpressedTSS.bed
    |--PlusOneDyad_SORT-Expression.bed
    |--PlusOneDyad_SORT-Expression_GROUP-Nuc-Dyad.bed
    |--PlusOneDyad_SORT-Expression_WithUnexpressed.bed
    |--PlusOneDyad_SORT-pHex-dHex_GROUP-Nuc-Dyad.bed
    |--PlusOneDyad_SORT-pHex-dHex_GROUP-Nuc-Dyad_GROUP-BOTTOM-1K.bed
    |--PlusOneDyad_SORT-pHex-dHex_GROUP-Nuc-Dyad_GROUP-TOP-1K.bed
    |--PlusOneDyad_SORT-pHN-dHN.bed
    |--PlusOneDyad_SORT-pHN-dHN_GROUP-BOTTOM-2500.bed
    |--PlusOneDyad_SORT-pHN-dHN_GROUP-TOP-2500.bed
    |--TSS_GROUP-All_SORT-CappedExpression.bed
    |--TSS_GROUP-All_SORT-CpG.bed
    |--TSS_GROUP-Expressed_SORT-CpG.bed
    |--TSS_GROUP-Expressed_SORT-Expression.bed
    |--TSS_GROUP-Unexpressed.bed
    |--200bp
      |--PlusOneDyad_SORT-pHex-dHex_GROUP-Nuc-Dyad_200bp.bed
    |--500bp
      |--PlusOneDyad_SORT-pHex-dHex_GROUP-Nuc-Dyad_GROUP-BOTTOM-1K_500bp.bed
      |--PlusOneDyad_SORT-pHex-dHex_GROUP-Nuc-Dyad_GROUP-TOP-1K_500bp.bed
    |--2000bp
      |--PlusOneDyad_SORT-Expression_2000bp.bed
      |--PlusOneDyad_SORT-Expression_GROUP-Nuc-Dyad_2000bp.bed
      |--PlusOneDyad_SORT-Expression_WithUnexpressed_2000bp.bed
      |--PlusOneDyad_SORT-pHN-dHN_2000bp.bed
      |--PlusOneDyad_SORT-pHN-dHN_GROUP-BOTTOM-2500_2000bp.bed
      |--PlusOneDyad_SORT-pHN-dHN_GROUP-TOP-2500_2000bp.bed
      |--TSS_GROUP-All_SORT-CappedExpression_2000bp.bed
      |--TSS_GROUP-All_SORT-CpG_2000bp.bed
      |--TSS_GROUP-Expressed_SORT-CpG_2000bp.bed
      |--TSS_GROUP-Expressed_SORT-Expression_2000bp.bed
      |--TSS_GROUP-Unexpressed_2000bp.bed
  |--RefPT-Motif
    |--<TF>_SORT-Occupancy.bed
    |--...
    |--1000bp
      |--<TF>_SORT-Occupancy_1000bp.bed
      |--...
    |--1bp
      |--<TF>_SORT-Occupancy_1bp.bed
      |--...
    |--250bp
      |--<TF>_SORT-Occupancy_250bp.bed
      |--...
    |--500bp
      |--<TF>_SORT-Occupancy_500bp.bed
      |--...
  |--RefPT-JASPAR
    |--<TF>_<JASPAR>_SORT-TFnucRatio_GROUP-Quartile*.bed
    |--...
    |--1000bp
      |--<TF>_<JASPAR>_SORT-TFnucRatio_GROUP-Quartile*_1000bp.bed
      |--...
  |--JASPAR
    |--<TF>_<JASPAR>.meme
  |--RefPT-JASPAR-nonK562
    |--<TF>_<JASPAR>_K562-specific-Unbound.bed
    |--<TF>_<JASPAR>_*_KD-*depleted.bed
    |--...
    |--1000bp
      |--<TF>_<JASPAR>_K562-specific-Unbound_1000bp.bed
      |--<TF>_<JASPAR>_*_KD-*depleted_1000bp.bed
      |--...
    |--150bp
      |--<TF>_<JASPAR>_K562-specific-Unbound_150bp.bed
      |--<TF>_<JASPAR>_*_KD-*depleted_150bp.bed
      |--...
  |--Conservation-SNP
    |--hg38.phyloP30way.bw
    |--dbSnp153_snv.bw
  |--sample-MEME
    |--<Target>.meme.txt
  |--BAM
    |--<Strain>_<Target>_<Assay>_<rep>_hg38.bam
    |--NormalizationFactors
      |--<BAMFILE>_NCISb_ScalingFactors.out
      |--<BAMFILE>_Total_ScalingFactors.out
 
```

</details>

### data/BAM/

Merged and renamed BAM files go here

#### data/BAM/NormalizationFactors/

Your normalization factors stored in `.txt` files go here.

### data/MEME/

MEME motif discovery results from BX K562 go here

### data/RefPT-Other/

Nucleosome and TSS centered reference files with various sorts and expansions

### data/RefPT-Motif/

Motif-centered reference files with various sorts and expansions.

### data/RefPT-Krebs/

Nucleosome and TSS RefPTs built for this manuscript.

> [!NOTE]
> Due to file size limits on Github, we split the `BNase-Nucleosomes.bed` file into two parts:
>
> ```sh
> wc -l BNase-Nucleosomes.bed 
> # 5607896 BNase-Nucleosomes.bed
> head -n 2500000 BNase-Nucleosomes.bed | gzip > BNase-Nucleosomes_part1.bed.gz
> tail -n 3107896 BNase-Nucleosomes.bed | gzip > BNase-Nucleosomes_part2.bed.gz
> ```
> 
> So you will need to reconstruct the file with:
> 
> ```sh
> gzip -d BNase-Nucleosomes_part1.bed.gz
> gzip -d BNase-Nucleosomes_part2.bed.gz
> cat BNase-Nucleosomes_part1.bed BNase-Nucleosomes_part2.bed > BNase-Nucleosomes.bed
> ```

### data/RefPT-JASPAR/

JASPAR-based motifs for benzonase cut-site analysis.

### data/RefPT-Other/

Externally retrieved RefPTs.

### data/FASTQ/

Raw sequencing datasets go here.

### data/hg38_files

The ref genome builds and annotations from the set-up the scripts in `00_Download_and_Preprocessing` populate this directory. It's required for the alignments in `00_Download_and_Preprocessing`, strain checks in `01_Run_GenoPipe`, and for running other sequence analyses and figure generation.

### data/mm10_files

Similar to `data/hg38_files` but for mouse datasets.

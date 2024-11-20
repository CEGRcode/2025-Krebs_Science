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
    |--<TF>_Bound-LowerHalf.bed
    |--...
    |--1000bp
      |--<TF>_Bound-LowerHalf_1000bp.bed
      |--...
```

</details>

### data/BAM/
merged and renamed BAM files go here

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

### data/RefPT-JASPAR/
JASPAR-based motifs for benzonase cut-site analysis.

### data/RefPT-Other/
Externally retrieved RefPTs.


### data/FASTQ/
Raw sequencing datasets go here.

### data/hg19_files
Run the following commands in the terminal from the `data/hg19_files` directory to download and index the reference genome. (Required to run alignments in `00_Download_and_Preprocessing`, strain checks in `01_Run_GenoPipe`, and for running other sequence analyses and figure generation)
```
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz -O hg19.fa.gz
gzip -d hg19.fa.gz         # hg19.fa
bwa index hg19.fa          # hg19.fa.amb, hg19.fa.ann, hg19.fa.bwt, hg19.fa.pac, hg19.fa.sa
samtools faidx hg19.fa     # hg19.fa.fai
```

...and there are some other reference files already provided for you within this repository.
```
hg19.info                   # used by ChExMix
hg19_exclude_contig.txt     # used by ChExMix
hg19_background_model.txt   # used by ChExMix
hg19_exclude.bed            # used by 00_Download_and_Preprocessing/3_normalize_samples.pbs
                            # and 02_Call_RefPT/2_FIMO_Motifs_from_Genome.pbs
```

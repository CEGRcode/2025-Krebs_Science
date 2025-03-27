
Call "lowly bound" motif reference points for motif-centered figures.

<details>
<summary> Full execution summary
</summary>

```
data
  |--JASPAR
    |--ATF7_MA0834-1.meme
    |--BACH1_MA1633-1.meme
    |--CTCF_MA1930-1.meme
    |--ELF1_MA0473-3.meme
    |--MAX_MA0058-3.meme
    |--MEIS2_MA1640-1.meme
    |--NFIC_MA1527-1.meme
    |--REST_MA0138-2.meme
    |--SP1_MA0079-5.meme
    |--SPI1_MA0080-6.meme
    |--ZKSCAN1_MA1585-1.meme
  |--RefPT-Motifs
    |--ATF7_Bound-LowerHalf.bed
    |--ATF7_Bound-NonPromoter.bed
    |--ATF7_Bound-Promoter.bed
    |--BACH1_Bound-LowerHalf.bed
    |--CTCF_Bound-LowerHalf.bed
    |--ELF1_Bound-LowerHalf.bed
    |--MAX_Bound-LowerHalf.bed
    |--MEIS2_Bound-LowerHalf.bed
    |--NFIC_Bound-LowerHalf.bed
    |--REST_Bound-LowerHalf.bed
    |--SP1_Bound-LowerHalf.bed
    |--SPI1_Bound-LowerHalf.bed
    |--ZKSCAN1_Bound-LowerHalf.bed
    |--1000bp/
      |--...<see above with _1000bp.bed suffix>
03_Call_Motifs
  |--narrowPeak
    |--ATF7_ENCFF868QLL.bed.gz
    |--BACH1_ENCFF423EMU.bed.gz
    |--CTCF_ENCFF738TKN.bed.gz
    |--ELF1_ENCFF392MUM.bed.gz
    |--MAX_ENCFF422NGZ.bed.gz
    |--MEIS2_ENCFF613RNG.bed.gz
    |--NFIC_ENCFF370ENX.bed.gz
    |--REST_ENCFF895QLA.bed.gz
    |--SP1_ENCFF300XUA.bed.gz
    |--SPI1_ENCFF664XPS.bed.gz
    |--ZKSCAN1_ENCFF163VUK.bed.gz
```

</details>

### 0a_Download_JAPAR_and_ENCODE_data.sh
Download PWMs from JASPAR (`.meme`) and ChIP binding peaks from ENCODE (`.bed.gz`)
```
sh 0_Download_JAPAR_and_ENCODE_data.sh
```
For each TF, there should be one of each the following files:
```
data/JASPAR/TF_MAXXXX-X.meme
03_Call_Motifs/narrowPeaks/TF_ENCFFXXXXXX.meme
```

### 1_FIMO_Motifs_from_Genome.sbatch
Run FIMO for each JASPAR motif to get all instances in the genome and filter to keep motifs that are at least 500bp from all blacklisted regions.
```
sbatch 1_FIMO_Motifs_from_Genome.sbatch
```
Results in four intermediate files for each TF:
```
FIMO/TF/fimo.gff
FIMO/TF/fimo.bed
FIMO/TF/fimo_1000bp.bed
FIMO/TF/filtered.bed
```

### 2a_Intersect_Motifs_wENCODE_ChIP-seq_peaks.sbatch
Get "lowly bound" motifs by intersecting the motif occurrences with ENCODE's ChIP-seq peaks and filtering for the bottom half of motifs (sorted by ENCODE binding score),
```
sbatch 2a_Intersect_Motifs_wENCODE_ChIP-seq_peaks.sbatch
```
Along with intermediate files, two final motif RefPT files are generated for each TF:
```
../data/RefPT-Motif/TF.bed
../data/RefPT-Motif/1000bp/TF_1000bp.bed
```

### 0b_Download_nonK562_ENCODE_data.sh
Download non-K562 ChIP binding peaks from ENCODE (`.bed.gz`)
```
sh 0b_Download_nonK562_ENCODE_data.sh
```
the downloaded bed file will be stored at 
```
03_Call_Motifs/NonK562_narrowPeaks/TF_Cellline_ENCFFXXXXXX.bed.gz
```

### 2b_Intersect_Motifs_wENCODE_ChIP-seq_peaks_nonK562.sh
get unbound motif of CTCF, CREM, ZKSCAN1, ATF2, ZNF263. Here the unbound could be devided into specificly-unbound in K562 and unbound in any cell line.
```
sbatch 2b_Intersect_Motifs_wENCODE_ChIP-seq_peaks_nonK562.sh
```
Along with intermediate files, two final motif RefPT files are generated for each TF:
```
../data/RefPT-JASPAR-nonK562/TF_K562-specific-Unbound.bed
../data/RefPT-JASPAR-nonK562/TF_Unbound.bed
../data/RefPT-JASPAR-nonK562/1000bp/TF_K562-specific-Unbound_1000bp.bed
../data/RefPT-JASPAR-nonK562/1000bp/TF_Unbound_1000bp.bed
```
### 2c_mm10_CTCF.sh
Get CTCF motif binding sites from mm10,split into 4 categories based on CTCF ChIP and MPE-seq in mouse cell line J1
```
sh 2c_mm10_CTCF.sh
```
Along with intermediate files, two final motif RefPT files are generated for each TF:
```
../data/RefPT-JASPAR-nonK562/MA1929_1_mm10_intersected_MPE-seq10min_164bp_category*.bed
../data/RefPT-JASPAR-nonK562/1000bp/MA1929_1_mm10_intersected_MPE-seq10min_164bp_category*_1000bp.bed
```

```

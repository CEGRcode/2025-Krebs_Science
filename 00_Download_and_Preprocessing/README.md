# 00_Download_and_Preprocessing

This directory includes the shell and SLURM scripts for downloading and preprocessing the sequencing samples for analysis throughout the rest of this repo.

The various omic processing can can be executed in any order after the numbered set-up scripts have been run. Make sure you execute these shell scripts from the `00_Download_and_Preprocessing` directory.

```sh
# Set-up
sh 0a_Setup_hg38_reffiles.sh
sh 0b_Setup_mm10_reffiles.sh
sh 0c_Download_Conservation-SNP_data.sh

# Internal data processing
sbatch 1_align-dedup_PEGR.sh
sbatch Benzonase-ChIP_merge.sbatch
sbatch Benzonase-seq_merge.sbatch
sbatch ChIP-exo-v6_merge.sh

# External data processing
sbatch ATAC-seq_download_ENCODE.sbatch
sbatch CoPRO_download-merge_ENCODE.sbatch
sbatch CUTRUN_download-align-dedup-merge.sbatch
sbatch DNase-FLASH_download-align-dedup-filter.sbatch
sbatch DNase-seq_download_ENCODE.sbatch

sbatch GROcap_download_ENCODE.sbatch
sbatch MNase-ChIP_download-align-dedup-filter-merge.sbatch
sbatch MNase-seq_download-align-filter-merge_ENCODE.sbatch # (deep, single-end, ABI SOLiD)
sbatch MNase-seq_download-align-filter-merge_SRR.sbatch  # (titrations)
sbatch MPE-seq_download-align-dedup-filter-merge.sbatch

# Post-processing constant calculations
sbatch X_get_scaling_factors.sbatch
```

Benzonase-seq
ChIP-exo (histone)
ChIP-exo (tf, w and w/o Benzonase)


<details>
<summary> Directory Structure
</summary>

```
|--00_Download_and_Preprocessing
  |--Trimmomatic-0.36
    |--adapters
      |--TruSeq3-PE.fa
      |--...
    |--trimmomatic-0.36.jar
    |--...
  |--CUTRUN
    |--...
  |--ENCODE
    |--...
  |--...all the submission scripts
|--data
  |--BAM
    |--ATAC-seq_ENCFF077FBI_hg38.bam
    |--CoPRO_Capped_merge_hg38.bam
    |--CUTRUN_CTCF_merge_hg38.bam
    |--CUTRUN_H2AZ_merge_hg38.bam
    |--CUTRUN_H3K4me1_merge_hg38.bam
    |--CUTRUN_H3K4me3_merge_hg38.bam
    |--CUTRUN_H3K27ac_merge_hg38.bam
    |--CUTRUN_H3K27ac-Millipore_merge_hg38.bam
    |--CUTRUN_H3K27me3_merge_hg38.bam
    |--CUTRUN_IgG_merge_hg38.bam
    |--DNase-FLASH_SRR801881_hg38.bam
    |--DNase-FLASH_SRR801880_hg38.bam
    |--DNase-seq_ENCFF425WDA_rep1_hg38.bam
    |--DNase-seq_ENCFF518XTC_rep1_hg38.bam
    |--GRO-cap_ENCODE_merge_hg38.bam
    |--MNase-ChIP_H3K4me3_merge_hg38.bam
    |--MNase-seq_ENCODE_merge_hg38.bam
    |--MNase-seq_5U_rep1_hg38.bam
    |--MNase-seq_21U_rep1_hg38.bam
    |--MNase-seq_79U_rep1_hg38.bam
    |--MNase-seq_304U_rep1_hg38.bam
    |--MPE-seq_10min_merge_rep1_mm10.bam
    |--MPE-seq_10min_rep2_mm10.bam
    |--MPE-seq_20min_merge_rep1_mm10.bam
    |--MPE-seq_20min_rep2_mm10.bam
    |--MPE-seq_30min_merge_rep1_mm10.bam
    |--MPE-seq_30min_rep2_mm10.bam
  |--hg38_files
    |--hg38.fa
  |--mm10_files
    |--mm10.fa
  |--RefPT-Other
    |--CpG_Islands.bed
  |--Conservation-SNP
    |--dbSnp153_snv.bw
    |--hg38.phyloP30way.bw

```

</details>


## Set-up scripts

### 0a_Setup_hg38_reffiles.sh

Downloads the hg38.fa genome and creates all the indexes needed to run the preprocessing scripts (`.fai` index, [bwa](https://bio-bwa.sourceforge.net/bwa.shtml) index, [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer) index, and [bowtie (1.2.3)](https://bowtie-bio.sourceforge.net/manual.shtml#the-bowtie-build-indexer) colorspace index).

### 0b_Setup_mm10_reffiles.sh

Downloads the mm10.fa genome and creates all the indexes needed to run the preprocessing scripts for mouse data (i.e. MPE-seq).

### 0c_Download_Conservation-SNP-CpG_data.sh

Download conservation and SNPs from USCS browser and reformat `../data/RefPT-Other/CpGIslands.tsv.gz`

### Set-up Trimmomatic

Need Trimmomatic set-up to process CUT&RUN data below.

```sh
curl http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip > Trimmomatic-0.36.zip
unzip Trimmomatic-0.36.zip
```


## Internal/Novel data processing

### 1_align-dedup_PEGR.sbatch

Script that mimics Galaxy core pipeline to perform initial alignment off of raw FASTQ data and remove duplicates. All internal data aligned with this script.

### Merging replicates

Merge BAM files for technical replicates (and IgG biological replicates) and rename BAM files to use standard naming system that includes experiment metadata (Cell line, IP target, antibody, assay-type, replicate, and genome build).

```
# Sample ids from PEGR to use with api download from PEGR
sampleIDs_Benzonase-ChIP.txt
sampleIDs_Benzonase-seq.txt
sampleIDs_ChIP-exo-v6.txt
```

## External/Published data processing

Each external dataset is downloaded and processed to adhere to the associated preprocessing standards resulting in a collection of .bam files for pileups.

## X_get_scaling_factors.sbatch

Both TotalTag and NCIS normalization factors are calculated for each `data/BAM/*hg38.bam` and saved to `data/BAM/NormalizationFactors/` with the name `*TotalTag_ScalingFactors.out` or `*_NCISb_ScalingFactors.out`. For the NCIS, a blacklist reference and IgG control BAMs that are cell-line and assay-specific are input based on the standard BAM filename structure (parse `_` delimited tokens for assay and cell line info).



This directory includes the shell and SLURM scripts for downloading and preprocessing the sequencing samples for analysis throughout the rest of this repo.

The following scripts should be executed in numerical order. The user should update the working directory (`$WRK`) filepath variable within every script before executing.

### `sampleIDs.txt`
The set of ids for samples from PEGR to analyze

### `MEMEsampleIDs.txt`
The set of ids for samples from PEGR to pull MEME motif analysis (PWMs) for.

### 0_download_data.sh
TODESCRIBE

```
bash 0_download_data.sh
```

### 1_align_data.sbatch
Script that mimics Galaxy core pipeline to perform initial alignment off of raw FASTQ data and remove duplicates.

```
# ^change the number of BAM files samples (SBATCH --array)
# To execute, type
sbatch 1_align_data.sbatch
##
# To check status, type
sbatch -u <myusername> -t
```

TODO: mimic core pipeline MEME motif calling

### 2_Merge_Replicates.sh
Merge BAM files for technical replicates (and IgG biological replicates) and rename BAM files to use standard naming system that includes experiment metadata (Cell line, IP target, antibody, assay-type, replicate, and genome build).

```
bash 2_Merge_Replicates.sh
```

### 3_normalize_samples.sbatch
Both TotalTag and NCIS normalization factors are calculated for each `data/BAM/SAMPLE.bam` and saved to `data/BAM/NormalizationFactors/` with the name `SAMPLE_NCISb.out`. For the NCIS, a blacklist reference and IgG control BAMs that are cell-line and assay-specific are input based on the standard BAM filename structure (parse `_` delimited tokens for assay and cell line info).

```
# ^change the number of BAM files samples (SBATCH --array)
# To execute, type
sbatch 3_normalize_samples.sbatch
##
# To check status, type
sbatch -u <myusername> -t
```
### 4_download_conservation-snp.sh
Download conservation and SNPs from USCS browser.
make sure proper scripts in bin
```
bash 4_download_conservation-snp.sh
```

The following scripts can can be executed in any order after `0_Setup_hg38_reffiles.sh` has been run. The user should update the working directory (`$WRK`) filepath variable within every script before executing.

1 - Benzonase-seq
2 - ChIP-exo (w and w/o Benzonase)
3 - CoPRO
4 - MNase-seq (titrations)
5 - MNase-seq (deep, single-end, ABI SOLiD)
6 - DNase-seq
7 - MNase-ChIP
8 - CUTRUN
9 - Global calculations

## 0_Setup_hg38_reffiles.sh

Downloads the hg38.fa genome and creates all the indexes needed to run the preprocessing scripts (`.fai` index, [bwa](https://bio-bwa.sourceforge.net/bwa.shtml) index, [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer) index, and [bowtie (1.2.3)](https://bowtie-bio.sourceforge.net/manual.shtml#the-bowtie-build-indexer) colorspace index).

```
sh 0_Setup_hg38_reffiles.sh
```

Previously worked with hg19 build. This build file is kept here for records:

```
sh 0_Setup_hg19_reffiles.sh
```

## Internal/Novel data processing

Download FASTQ/BAM files from PEGR using the [EGC_utility_scripts](https://github.com/CEGRcode/EGC_utility_scripts) repo's `generate_BAM_file_from_PEGR.py`. You may modify this to also include other data downloaded from public repositories.

```
# Sample ids from PEGR to use with api download scripts
benzonase_samples.txt
chip_samples.txt
```

Then replicate can be merged and renamed to standard file names using the following:

```
# Aggregate and dedup replicates
sbatch Benzonase-ChIP_download-merge-dedup-merge.sbatch
sbatch Benzonase-seq_download-merge-dedup-merge.sbatch
```

## External/Published data processing

Describe how to align data or use scripts saved here to align data.

```
sbatch CUTRUN_download-align-dedup-merge.sbatch
sbatch CoPRO_download-ENCODE.sbatch
sbatch DNase-seq_download.sbatch
sbatch MNase-ChIP_download-align-dedup-filter-merge.sbatch
sbatch MNase-seq-ENCODE_download-align-filter-merge.sbatch
sbatch MNase-seq-Titrations_download-align-dedup-filter.sbatch
```

## X_get_scaling_factors.sbatch

Total tag normalization was performed to get appropriate scaling factors for the data.

```
# ^change the number of BAM files samples (PBS -t)
# To execute, type
# sbatch X_get_scaling_factors.sbatch
##
# To check status, type
# squeue -u <myusername> -t
```

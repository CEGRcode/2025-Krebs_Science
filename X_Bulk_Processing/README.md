# X_Bulk_Processing

This directory contains scripts to pileup select BAM files against select RefPT BED file and saves the analysis to `2026-Chen_NucleosomePhasing/Library`.

Please make sure your BAM files are downloaded and indexed in the appropriate directory (`data/BAM/`) with normalization factors calculated (`data/BAM/NormalizationFactors/`) and reference points built (`data/RefPT-XXX`).

See `00_Download_and_preprocessing` to read more about how BAM files were built.
See `02_Call_Nucleosomes`, `03_Call_JASPAR` and `04_Call_Motifs` to read more about how RefPT files were built.

Scripts were designed to run on linux systems with a Slurm scheduler. You may need to adjust the `#SBATCH` parameters according to your HPC system's configurations. You may run a subset of job arrays by setting the `#SBATCH --array 1,2-10` to the range of samples you wish to execute.

- [Read more on Slurm Job arrays](https://slurm.schedmd.com/job_array.html)
- No Slurm scheduler? `SLURM_ARRAY_TASK_ID` can be hardcoded in a for loop to run them sequentially.

## 1a_make_weblogos.sh

```sh
sh 1a_make_weblogos.sh
```

## 1b_bulk_4-color-plot.sh

```sh
sh 1b_bulk_4-color-plot.sh
```

## 2x_MotifAnalyses.sbatch

```sh
sbatch 2a_MotifAnalyses.sbatch
sbatch 2b_unboundsites_otherassay_Pileups.sbatch
sh 2c_10bp_Dinucleotide_4Q.sh
sh 2d_10bp_midflank_DNAshape_4Q.sh
```

## 3_Midpoint_Pileups.sbatch

```sh
sbatch 3_Midpoint_Pileups.sbatch
```

## 4x_Five_Read1_Pileups

```sh
sbatch 4_Five_Read1_Pileups.sbatch
sbatch 4b_10bp_pileups.sh
```

## 5_Five_Read2_Pileups.sbatch

```sh
sbatch 5_Five_Read2_Pileups.sbatch
```

## 6_CoPRO_TSS_ActiveSite_Pileups.sbatch

```sh
sbatch 5_CoPRO_TSS_ActiveSite_Pileups.sbatch
```

## 7x_Convervation/dbSnp153_pileups

```sh
sh 7a_dbSnp153_pileups.sbatch
sh 7b_Conservation-phylo_pileups.sh
```

## 8x_CTCF10phase

```sh
sbatch 8a_10phase_of_CTCF_Q1_Q4.sbatch
sh 8b_CTCF_phaseallign.sh
sbatch 8c_10phase_of_CTCF_Q1_both.sbatch
```

## 9x_FoxA10phase

```sh
sbatch 9a_10phase_of_FOXA_uHepG2.sbatch
sh 9b_FOXA_phaseallign.sh
```

## 10x_NFIA10phase

```sh
sbatch 10a_10phase_of_NFIA_unbound.sbatch
sh 10b_NFIA_phaseallign.sh
```

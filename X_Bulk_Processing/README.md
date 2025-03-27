This directory contains scripts to pileup select BAM files against select RefPT BED file and saves the analysis to `2024-Krebs_Science/Library`.

Please make sure your BAM files are downloaded and indexed in the appropriate directory (`data/BAM/`) with normalization factors calculated (`data/BAM/NormalizationFactors/`) and reference points built (`data/RefPT-XXX`).

The user should update the working directory ($WRK) filepath variable within every script before executing

See `00_Download_and_preprocessing` to read more about how BAM files were built.
See `02_Call_Nucleosomes`, `03_Call_JASPAR` and `04_Call_Motifs` to read more about how RefPT files were built.

Scripts were designed to run on linux systems with a Slurm scheduler. Simply update `WRK` to match the directory path where you cloned this repo and set the `#SBATCH --array 1,2-10` to the range you wish to execute.

- [Read more on Slurm Job arrays](https://slurm.schedmd.com/job_array.html)
- No Slurm scheduler? `SLURM_ARRAY_TASK_ID` can be hardcoded and run as a regular shell script

## 1a_make_weblogos.sh

```
sh 1a_make_weblogos.sh
```

## 1b_bulk_4-color-plot.sh

```
sh 1b_bulk_4-color-plot.sh
```

## 2x_MotifAnalyses.sbatch

```
sbatch 2a_MotifAnalyses.sbatch
sbatch 2b_unboundsites_otherassay_Pileups.sbatch
sh 2c_10bp_Dinucleotide_4Q.sh
sh 2d_10bp_midflank_DNAshape_4Q.sh
```

## 3_Midpoint_Pileups.sbatch

```
sbatch 3_Midpoint_Pileups.sbatch
```

## 4x_Five_Read1_Pileups

```
sbatch 4_Five_Read1_Pileups.sbatch
sbatch 4b_10bp_pileups.sh
```

## 5_Five_Read2_Pileups.sbatch

```
sbatch 5_Five_Read2_Pileups.sbatch
```

## 6_CoPRO_TSS_ActiveSite_Pileups.sbatch

```
sbatch 5_CoPRO_TSS_ActiveSite_Pileups.sbatch
```

## 7x_Convervation/dbSnp153_pileups

```
sh 7a_dbSnp153_pileups.sbatch
sh 7b_Conservation-phylo_pileups.sh
```

## 8x_CTCF10phase

```
sbatch 8a_10phase_of_CTCF_Q1_Q4.sbatch
sh 8b_CTCF_phaseallign.sh
sbatch 8c_10phase_of_CTCF_Q1_both.sbatch
```

## 9x_FoxA10phase

```
sbatch 9a_10phase_of_FOXA_uHepG2.sbatch
sh 9b_FOXA_phaseallign.sh
```


## 10x_NFIA10phase

```
sbatch 10a_10phase_of_NFIA_unbound.sbatch
sh 10b_NFIA_phaseallign.sh
```

## 11_cut-site_nucleotide-content.sh

```
sh 11_cut-site_nucleotide-content.sh
```

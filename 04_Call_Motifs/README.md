The following scripts should be executed in numerical order. The user should update the working directory ($WRK) filepath variable within every script before executing

### 1_Get_PWM.sh
Create Motif RefPT files for figure pileups and weblogos and save the PWM for identifying motif sites from MEME.
```
bash 1_Get_PWM.sh
```
copies motif PWM files into the `PWM` directory with `<NAME>.meme.txt` format.

### 2_FIMO_Motifs_from_Genome.sbatch
Call motifs in genome from PWM and get blacklist filtered & exlusion zone filtered bed coordinates.
```
# Update submission info (e.g. replace NSAMPLES with integer)
sbatch 2_FIMO_Motifs_from_Genome.sbatch
```

### 3_Filter_and_Sort_by_occupancy.sbatch
Filter motifs to include "bound" sets and sort by occupancy.
```
# Update submission info (e.g. replace NSAMPLES with integer)
sbatch 3_Filter_and_Sort_by_occupancy.sbatch
```
For every FIMO result for each TF, reformat into BED and call "bound" sets by BAM file metadata parsing.

### 4_NFIA_motif.sh
Perform NFIA-specific transformations, groups, and sorts for figure RefPT building.
```
# Update submission info (e.g. add WRK path)
sh 4_NFIA_motif.sh
```

### 5_FoxA_motif.sh
Perform FoxA-specific transformations, groups, and sorts for figure RefPT building.
```
# Update submission info (e.g. add WRK path)
sh 5_FoxA_motif.sh
```

### 6_CTCF_motif.sh
Perform CTCF-specific transformations, groups, and sorts for figure RefPT building.
```
# Update submission info (e.g. add WRK path)
sh 6_CTCF_motif.sh
```

```

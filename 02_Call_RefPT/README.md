The following scripts should be executed in numerical order. The user should update the working directory ($WRK) filepath variable within every script before executing

### 1_Get_PWM.sh
Create Motif RefPT files for figure pileups and save the PWM and scrambled_PWM files for identifying motif sites from MEME.
```
bash 1_Get_PWM.sh
```
copies motif PWM files into the `PWM` directory with `<NAME>.meme.txt` format.

### 2_FIMO_Motifs_from_Genome.pbs
Call motifs (K562-BX and scrambled) in genome from PWM and get filtered bed coordinate. 
```
# Update submission info (e.g. replace NSAMPLES with integer)
sbatch 2_FIMO_Motifs_from_Genome.sbatch
```
Generates filtered, BED-formatted motif files under `FIMO/<TF>_motif1_unsorted.bed` to be proccessed further (see `3_Filter_and_Sort_by_occupancy.sbatch` and `4_scrambled_top10K.sh`). 

### 3_Filter_and_Sort_by_occupancy.sbatch
Filter motifs to include "bound" sets and sort by occupancy (K562-BX-rep1).
```
# Update submission info (e.g. replace NSAMPLES with integer)
qsub 4_Filter_and_Sort_by_occupancy.pbs
```
For every TF of interest, the above scripts should make
```
../data/RefPT-Motif/<TF>_Occupancy.bed
../data/RefPT-Motif/1000bp/<TF>_Occupancy_1000bp.bed
../data/RefPT-Motif/1bp/<TF>_Occupancy_1bp.bed
```

### 4_scrambled_top10K.sh
take top 10K scrambled sites based on fimo score
```
bash 4_scrambled_top10K.sh
```
For each scrambled motif, the above scripts should make
``` 
../data/RefPT-Motif/1000bp/scrambled<TF>_top10K_1000bp.bed
../data/RefPT-Motif/1bp/scrambled<TF>_top10K_1bp.bed
../data/RefPT-Motif/scrambled<TF>_top10K.bed
```

### 5_WDR5TSS_sort.sh
specific cutouff for WDR5 binding motif or TSS
```
# Update submission info (e.g. determine the location of tss_all_k562.bed from Core et all and WDR5_MOTIF1_sorted.bed)
sh 5_WDR5TSS_sort.sh
```
the above scripts should make
```
../data/BED/TSS95_WDR5_sort_1000bp.bed
../data/BED/TSS_WDR5_same_1000bp.bed
../data/RefPT-Motif/1000bp/WDR5_Occupancy_1000bp.bed
../data/RefPT-Motif/1000bp/WDR5_TSS_oppo_1000bp.bed
../data/RefPT-Motif/1bp/WDR5_Occupancy_1bp.bed
../data/RefPT-Motif/WDR5_Occupancy.bed
../data/RefPT-Motif/WDR5_TSS_same_32bp.bed
../data/RefPT-Motif/WDR5_TSS_oppo_32bp.bed
```

### 6_CTCF_ThirdsByOccupancy.sh
break CTCF sites into 3 class based on occypancy level
```
sh 6_CTCF_ThirdsByOccupancy.sh
```
the above script should make thre new RefPTs
``` 
../data/RefPT-Motif/1000bp/CTCF_Occupancy_TOP_1000bp.bed
../data/RefPT-Motif/1000bp/CTCF_Occupancy_MIDDLE_1000bp.bed
../data/RefPT-Motif/1000bp/CTCF_Occupancy_BOTTOM_1000bp.bed
```

### 7_Motif_Sort_by_Nucleosome.sbatch
Have all the RefPT_1bp.bed resorted by their distance to the nearest nucleosome
```
sbatch 7_Motif_Sort_by_Nucleosome.sbatch
```
For every TF of interest, the above scripts should make
```
../data/RefPT-Motif/1bp/<TF>_NucSort_1bp.bed
../data/RefPT-Motif/1000bp/<TF>_NucSort_1000bp.bed
```

### 8_NFIA_ThirdsByNucPosition.sh
Break NFIA into 3 classes based on their relative location to neareast nucleosome
```
# Update submission info (e.g. determin location of NFIA_Occupancy_1bp_Nuc_sort.bed)
bash 8_NFIA_ThirdsByNucPosition.sh
```
For NFIA motif sites, the above scripts should make
```
../data/RefPT-Motif/1000bp/NFIA_NucSort-DOWNSTREAM_1000bp.bed
../data/RefPT-Motif/1000bp/NFIA_NucSort-OVERLAP_1000bp.bed
../data/RefPT-Motif/1000bp/NFIA_NucSort-UPSTREAM_1000bp.bed
../data/RefPT-Motif/100bp/NFIA_NucSort-DOWNSTREAM_100bp.bed
../data/RefPT-Motif/100bp/NFIA_NucSort-OVERLAP_100bp.bed
../data/RefPT-Motif/100bp/NFIA_NucSort-UPSTREAM_100bp.bed
../data/RefPT-Motif/1bp/NFIA_NucSort-DOWNSTREAM_1bp.bed
../data/RefPT-Motif/1bp/NFIA_NucSort-OVERLAP_1bp.bed
../data/RefPT-Motif/1bp/NFIA_NucSort-UPSTREAM_1bp.bed
```


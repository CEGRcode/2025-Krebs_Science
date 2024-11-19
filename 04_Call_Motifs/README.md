The following scripts should be executed in numerical order. The user should update the working directory ($WRK) filepath variable within every script before executing

### 1_Get_PWM.sh
Create Motif RefPT files for figure pileups and save the PWM and scrambled_PWM files for identifying motif sites from MEME.
```
bash 1_Get_PWM.sh
```
copies motif PWM files into the `PWM` directory with `<NAME>.meme.txt` format.

### 2_FIMO_Motifs_from_Genome.sbatch
Call motifs (K562-BX and scrambled) in genome from PWM and get blacklist filtered & exlusion zone filtered bed coordinates.
```
# Update submission info (e.g. replace NSAMPLES with integer)
sbatch 2_FIMO_Motifs_from_Genome.sbatch
```

### 3_Filter_and_Sort_by_occupancy.sbatch
Filter motifs to include "bound" sets and sort by occupancy (K562, BX).
```
# Update submission info (e.g. replace NSAMPLES with integer)
sbatch 3_Filter_and_Sort_by_occupancy.sbatch
```
For every TF of interest (except WDR5), the above script should make:
```
../data/RefPT-Motif/TF_Occupancy.bed
../data/RefPT-Motif/1bp/TF_Occupancy_1bp.bed
../data/RefPT-Motif/1000bp/TF_Occupancy_1000bp.bed
```

### 4_scrambled_top10k.sbatch
Filter scrambled motifs to keep top 10k motif matches based on MEME score (K562, BX).
```
# Update submission info (e.g. add WRK path)
sh 4_scrambled_top10k.sh
```
For every scrambledTF of interest, the above script should make:
```
../data/RefPT-Motif/scrambledTF_Occupancy.bed
../data/RefPT-Motif/1bp/scrambledTF_Occupancy_1bp.bed
../data/RefPT-Motif/1000bp/scrambledTF_Occupancy_1000bp.bed
```

### 5_WDR5_and_TSS-sort.sh
Call both WDR5_Occupancy and TSS sorted by distance to WDR5 reference points.
```
# Update submission info (e.g. add WRK path)
sh 5_WDR5TSS_sort.sh
```
The above script should make:
```
../data/RefPT-Motif/WDR5_Occupancy.bed
../data/RefPT-Motif/1bp/WDR5_Occupancy_1bp.bed
../data/RefPT-Motif/1000bp/WDR5_Occupancy_1000bp.bed
../data/RefPT-Other/TSS_DistWDR5.bed
../data/RefPT-Other/TSS_DistWDR5_1000bp.bed
../data/RefPT-Other/TSS_DistWDR5_Filter-SameStrand_1000bp.bed
```

### 6_CTCF_ThirdsByOccupancy.sh
Split CTCF sites into thirds based on occupancy (TOP, MIDDLE, BOTTOM).
```
# Update submission info (e.g. add WRK path)
sh 6_CTCF_ThirdsByOccupancy.sh
```
The above script should make:
```
../data/RefPT-Motif/1000bp/CTCF_Occupancy_TOP_1000bp.bed
../data/RefPT-Motif/1000bp/CTCF_Occupancy_MIDDLE_1000bp.bed
../data/RefPT-Motif/1000bp/CTCF_Occupancy_BOTTOM_1000bp.bed
```

### 7_Motif_Sort_by_Nucleosome.sbatch
For each `TF.meme.txt` in `PWM/`, make a RefPT file sorted by distance to closes Nucleosome coordinates (`../data/RefPT-Other/1bp/hg19_K562_Nucleosome_1bp_SORT-genomic_nostrand.bed`)
```
# Update submission info (e.g. replace NSAMPLES with integer)
sbatch 7_Motif_Sort_by_Nucleosome.sbatch
```
For every PWM, the above script should make:
```
../data/RefPT-Motif/1bp/TF_NucSort_1bp.bed
../data/RefPT-Motif/1000bp/TF_NucSort_1000bp.bed
```

### 8_NFIA_ThirdsByNucPosition.sh
Split Nucleosome distance-sorted NFIA motifs (`NFIA_NucSort_1bp.bed`) into three sections by how far upstream, downstream, or whether the nucleosome is overlapping the motif or not. (-73 <= overlap <= 73)
```
# Update submission info (e.g. add WRK path)
sh 8_NFIA_ThirdsByNucPosition.sh
```
The above script should make:
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
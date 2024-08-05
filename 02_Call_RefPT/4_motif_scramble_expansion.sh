#!/bin/bash
module load anaconda
source activate bioinfo
# Sort each FIMO'd scrambled motif by FIMO score and take the top 10K for each
# motif analysis 

### CHANGE ME
WRK=/path/to/2024-Chen_Nature/02_Call_RefPT
###

# Dependencies
# - java

# Script shortcuts
WRK=/storage/group/bfp2/default/hxc585_HainingChen/
FIMO=$WRK/Fox_NFIA_CTCF/02_Call_RefPT/FIMO/
SCRIPTMANAGER=$WRK/2023_Chen_PIC3/bin/ScriptManager-v0.14.jar
GENOME=$WRK/2023_Chen_PIC3/data/hg19_files/hg19.fa
Genome=$WRK/2023_Chen_PIC3/hg19_files/hg19.info
Nuc="$WRK/2024-Chen_Benzonase-ChIP-exo/03_RfMotif_sort_Nuc/K562_nuc_uHex_uTetra_1bp_sort_nostrand.bed"

cd $FIMO
# Inputs and outputs
BLACKLIST=$WRK/2023_Chen_PIC3/data/hg19_files/hg19_exclude.bed
MOTIF=$WRK/Fox_NFIA_CTCF/data/RefPT-Motif
[ -d $MOTIF ] || mkdir $MOTIF
[ -d $MOTIF/1000bp ] || mkdir $MOTIF/1000bp
[ -d $MOTIF/1bp ] || mkdir $MOTIF/1bp


# Loop through each scrambled* motif
for PWM in ../PWM/scrambled*.meme.txt;
do
    # Parse TF from PWM filename
    TF=`basename $PWM ".meme.txt"`
    # Sort by FIMO score and get top 10K
    sort -k5,5nr $TF/$TF\_motif1_unsorted.bed | head -n 10000 >  $TF/$TF\_top10k.bed
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1 $TF/$TF\_top10k.bed -o $MOTIF/1bp/$TF\_top10k_1bp.bed
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $TF/$TF\_top10k.bed -o $MOTIF/1000bp/$TF\_top10k_1000bp.bed
done

for FILE in */*_Occupancy.bed ;
do
    filename=`basename $FILE ".bed"`
    # Expand by 1bp and 1000bp
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1 $FILE -o $MOTIF/1bp/$filename\_1bp.bed
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $FILE -o $MOTIF/1000bp/$filename\_1000bp.bed
done


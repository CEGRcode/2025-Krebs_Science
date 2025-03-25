#!/bin/bash

# Make weblogos for every JASPAR MEME file in `data/JASPAR/*.meme` and MEME file in `02_Call_RefPT/PWM/*.meme.txt` (both strands)

### CHANGE ME
WRK=/storage/group/bfp2/default/hxc585_HainingChen/2025_Chen_TF-Nuc/X_Bulk_Processing
#WRK=/storage/group/bfp2/default/owl5022-OliviaLang/2024-Krebs_Science/X_Bulk_Processing
#WRK=/scratch/owl5022/2024-Krebs_Science/X_Bulk_Processing
###

# Dependencies
# - ceqlogo

set -exo
module load anaconda3
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Define the source files as LIST
LIST=`ls $WRK/../04_Call_Motifs/PWM/*_M1.meme.txt $WRK/../data/JASPAR/*.meme`

# Inputs and outputs
OUTDIR=$WRK/Library/WebLogos

# Set up output directories
[ -d Library ] || mkdir Library
[ -d $OUTDIR ] || mkdir $OUTDIR

# Loop through the PWM files
for PWMFILE in ${LIST[*]} ;
do
    BASE=`basename $PWMFILE ".txt"`
    BASE=`basename $BASE ".meme"`

    # Generate logo (reg and rc)
    ceqlogo -i $PWMFILE -m 1    -o $OUTDIR/$BASE\_logo.eps -f EPS
    ceqlogo -i $PWMFILE -m 1 -r -o $OUTDIR/$BASE\_logoRC.eps -f EPS
done

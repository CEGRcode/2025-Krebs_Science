#!/bin/bash

# Make weblogos for every MEME file in `02_Call_RefPT/PWM/*.meme.txt` (both strands)

### CHANGE ME
WRK=/path/to/2024-Krebs_Science/X_Bulk_Processing
WRK=/storage/home/owl5022/scratch/2024-Krebs_Science/X_Bulk_Processing
###

# Dependencies
# - ceqlogo

set -exo
module load anaconda3
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Define the source files as LIST
LIST=(
    "$WRK/../02_Call_RefPT/PWM/CTCF_M1.meme.txt"
    "$WRK/../02_Call_RefPT/PWM/FOXA2_M1.meme.txt"
    "$WRK/../02_Call_RefPT/PWM/NFIA_M1.meme.txt"
    "$WRK/../02_Call_RefPT/PWM/ZKSCAN1_M1.meme.txt"
)

# Inputs and outputs
OUTDIR=$WRK/Library/WebLogos

# Set up output directories
[ -d Library ] || mkdir Library
[ -d $OUTDIR ] || mkdir $OUTDIR

# Loop through the PWM files
for PWMFILE in ${LIST[*]} ;
do
    BASE=`basename $PWMFILE ".meme.txt"`
    # Generate logo (reg and rc)
    ceqlogo -i $PWMFILE -m 1 -r -o $OUTDIR/$BASE\_logo.eps -f EPS
    ceqlogo -i $PWMFILE -m 1 -r -o $OUTDIR/$BASE\_logoRC.eps -f EPS
done


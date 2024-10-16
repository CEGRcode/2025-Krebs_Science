#!/bin/bash

# Make weblogos for every MEME file in `02_Call_RefPT/PWM/*.meme.txt`

### CHANGE ME
WRK=/PATH/to/Title/
Library=$WRK/Library
[ -d $Library ] || mkdir -p $Library

###

# Dependencies
# - ceqlogo

set -exo
module load anaconda3
source activate meme

# Fill in placeholder constants with your directories
PWM=$WRK/02_Call_RefPT/PWM
OUTDIR=$Library/WebLogos

# Set up output directories
[ -d $OUTDIR ] || mkdir -p $OUTDIR

# Loop through the PWM files
for PWMFILE in $PWM/*.meme.txt;
do
    BASE=`basename $PWMFILE ".meme.txt"`
    # Generate logo
    ceqlogo -i $PWMFILE -m 1 -o $OUTDIR/$BASE\_logo.eps -f EPS
done

# make logo for ZKSCAN1
mkdir -p $OUTDIR/ZKSCAN1_MEME 
meme $WRK/Library/ZKSCAN1_Occupancy_flip_shift2_1000bp/ZKSCAN1_Occupancy_flip_shift2_1000bp_31bp.fa -oc $OUTDIR/ZKSCAN1_MEME -nmotifs 1 -maxw 31 -dna


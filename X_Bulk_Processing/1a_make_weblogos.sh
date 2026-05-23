#!/bin/bash

# Make weblogos for every JASPAR MEME file in `data/JASPAR/*.meme` and MEME file in `02_Call_RefPT/PWM/*.meme.txt` (both strands)

# Dependencies
# - ceqlogo

set -exo
source activate bx

# Define the source files as LIST
LIST=`ls ../04_Call_Motifs/PWM/*_M1.meme.txt ../data/JASPAR/*.meme`

# Inputs and outputs
OUTDIR=Library/WebLogos

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

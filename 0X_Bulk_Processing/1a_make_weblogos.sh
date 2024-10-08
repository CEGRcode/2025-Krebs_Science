#!/bin/bash

# Make weblogos for every MEME file in `02_Call_RefPT/PWM/*.meme.txt`

### CHANGE ME
WRK=/storage/group/bfp2/default/hxc585_HainingChen/
Library=$WRK/Fox_NFIA_CTCF/Library
[ -d $Library ] || mkdir -p $Library

###

# Dependencies
# - ceqlogo

set -exo
module load anaconda3
source activate meme

# Fill in placeholder constants with your directories
PWM=$WRK/Fox_NFIA_CTCF/02_Call_RefPT/PWM
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


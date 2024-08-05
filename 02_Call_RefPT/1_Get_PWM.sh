#!/bin/bash

# Organize select MEME reference files for RefPT building into the PWM directory

### CHANGE ME
WRK=/path/to/2024-Chen_Nature/02_Call_RefPT
###

# Inputs and outputs
MDIR=$WRK/../data/sample-MEME

# Create output directories if they don't exist
[ -d PWM ] || mkdir PWM

# Hardcode move of MEME files from PEGR workflow to raname as TF.meme.txt file
cp $MDIR/30478_CTCF_07-729_K562_-_-_Lysis-50Unuclease-4cycSonic-NoPapainDigestion-bothAmpureSteps_BX.meme.txt PWM/CTCF.meme.txt
cp $MDIR/32116_NFIA_HPA008884_K562_-_-_-_BX.meme.txt PWM/NFIA.meme.txt
cp $MDIR/32354_FoxA1_ab23738_HepG2_-_-_-_BX.meme.txt PWM/FoxA1.meme.txt


# Move scrambled motifs
cp $MDIR/handshuffled-30478_CTCF.meme.txt PWM/scrambledCTCF.meme.txt
cp $MDIR/handshuffled-32116_NFIA.meme.txt PWM/scrambledNFIA.meme.txt

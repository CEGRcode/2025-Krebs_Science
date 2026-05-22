#!/bin/bash

# Organize select MEME reference files for RefPT building into the PWM directory

# Dependencies
# - MEME suite (MEME)

set -exo
source activate bx

# Inputs and outputs
MDIR=../data/sample-MEME

# Create output directories if they don't exist
[ -d PWM ] || mkdir PWM

# # Hardcode move of MEME files from PEGR workflow to raname as TF.meme.txt file
# cp $MDIR/33924_CTCF_07-729_K562_-_IMDM_-_BX.meme.txt PWM/CTCF.meme.txt
# cp $MDIR/38480_FOXA2_ab256493_HepG2_-_-_-_BX.meme.txt PWM/FOXA2.meme.txt
# cp $MDIR/32116_NFIA_HPA008884_K562_-_-_-_BX.meme.txt PWM/NFIA.meme.txt
# cp $MDIR/34158_HNF4A_HPA004712_HepG2_-_-_-_BX.meme.txt PWM/HNF4A.meme.txt

# Extract the first motif
meme-get-motif -id 1 PWM/CTCF.meme.txt > PWM/CTCF_M1.meme.txt
meme-get-motif -id 1 PWM/FOXA2.meme.txt > PWM/FOXA2_M1.meme.txt
meme-get-motif -id 1 PWM/NFIA.meme.txt > PWM/NFIA_M1.meme.txt
meme-get-motif -id 1 PWM/HNF4A.meme.txt > PWM/HNF4A_M1.meme.txt
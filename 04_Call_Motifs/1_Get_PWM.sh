#!/bin/bash

# Organize select MEME reference files for RefPT building into the PWM directory

### CHANGE ME
WRK=/path/to/2024-Krebs_Science/04_Call_Motifs
WRK=/storage/home/owl5022/scratch/2024-Krebs_Science/04_Call_Motifs
###

# Dependencies
# - MEME suite (MEME)

set -exo
module load anaconda3
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Inputs and outputs
MDIR=$WRK/data/sample-MEME

# Create output directories if they don't exist
[ -d PWM ] || mkdir PWM

# Hardcode move of MEME files from PEGR workflow to raname as TF.meme.txt file
cp $MDIR/33924_CTCF_07-729_K562_-_IMDM_-_BX.meme.txt PWM/CTCF.meme.txt
cp $MDIR/38480_FOXA2_ab256493_HepG2_-_-_-_BX.meme.txt PWM/FOXA2.meme.txt
cp $MDIR/32116_NFIA_HPA008884_K562_-_-_-_BX.meme.txt PWM/NFIA.meme.txt
cp $MDIR/34048_ZKSCAN1_HPA006672_K562_-_IMDM_-_BX.meme.txt PWM/ZKSCAN1.meme.txt

meme-get-motif -id 1 PWM/CTCF.meme.txt > PWM/CTCF_M1.meme.txt
meme-get-motif -id 1 PWM/FOXA2.meme.txt > PWM/FOXA2_M1.meme.txt
meme-get-motif -id 1 PWM/NFIA.meme.txt > PWM/NFIA_M1.meme.txt
meme-get-motif -id 1 PWM/ZKSCAN1.meme.txt > PWM/ZKSCAN1_M1.meme.txt
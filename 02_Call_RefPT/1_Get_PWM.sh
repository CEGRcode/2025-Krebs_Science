#!/bin/bash

# Organize select MEME reference files for RefPT building into the PWM directory

### CHANGE ME
WRK=/path/to/2023-Chen_Benzonase-ChIP-exo/02_Call_RefPT
###

# Inputs and outputs
MDIR=$WRK/../data/sample-MEME

# Create output directories if they don't exist
[ -d PWM ] || mkdir PWM

# Hardcode move of MEME files from PEGR workflow to raname as TF.meme.txt file
cp $MDIR/30478_CTCF_07-729_K562_-_-_Lysis-50Unuclease-4cycSonic-NoPapainDigestion-bothAmpureSteps_BX.meme.txt PWM/CTCF.meme.txt
cp $MDIR/32241_E4F1_A300-832A_K562_-_-_-_BX.meme.txt PWM/E4F1.meme.txt ### Eeeek!
cp $MDIR/32114_GABPA_HPA003258_K562_-_-_-_BX.meme.txt PWM/GABPA.meme.txt
cp $MDIR/32096_MEIS2_1A11_K562_-_-_-_BX.meme.txt PWM/MEIS2.meme.txt
cp $MDIR/32116_NFIA_HPA008884_K562_-_-_-_BX.meme.txt PWM/NFIA.meme.txt
cp $MDIR/31970_NRF1_3H1-s_K562_-_-_-_BX.meme.txt PWM/NRF1.meme.txt
cp $MDIR/32120_Sp1_HPA001853_K562_-_-_-_BX.meme.txt PWM/SP1.meme.txt
cp $MDIR/32849_SRF_A303-172A_K562_-_-_-_BX.meme.txt PWM/SRF.meme.txt
cp $MDIR/29261_USF1_170427_K562_-_-_-_BX.meme.txt PWM/USF1.meme.txt
cp $MDIR/32201_USF2_1A11_K562_-_-_-_BX.meme.txt PWM/USF2.meme.txt
cp $MDIR/33918_WDR5_HPA047182_K562_-_IMDM_-_BX.meme.txt PWM/WDR5.meme.txt
cp $MDIR/32097_YY1_1B2_K562_-_-_-_BX.meme.txt PWM/YY1.meme.txt

# Move scrambled motifs
cp $MDIR/handshuffled-30478_CTCF.meme.txt PWM/scrambledCTCF.meme.txt
cp $MDIR/handshuffled-32116_NFIA.meme.txt PWM/scrambledNFIA.meme.txt
cp $MDIR/handshuffled-_WDR5.meme.txt PWM/scrambledWDR5.meme.txt

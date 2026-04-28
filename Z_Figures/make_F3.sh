#!/bin/bash

# Organize data from X_Bulk_Processing into Z_Figures for F4

### CHANGE ME
WRK=/storage/group/bfp2/default/hxc585_HainingChen/2025_Chen_TF-Nuc/Z_Figures
#WRK=/storage/home/owl5022/scratch/2024-Krebs_Science/Z_Figures
###

set -exo
module load anaconda3
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

LIBRARY=$WRK/../X_Bulk_Processing/Library
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar

[ -d F3 ] || mkdir F3

[ -d F3/A ] || mkdir F3/A
# ===============================================================================================================================


# heatmap
cp $LIBRARY/10phase/CTCF_Q4/10xplot/*read1_original_phase_5.out F3/A
cp $LIBRARY/10phase/CTCF_Q4/10xplot/*read1_shuffle_phase_5.out F3/A
cp $LIBRARY/10phase/CTCF_Q4/10xplot/*read1_original_phase_1.out F3/A
cp $LIBRARY/10phase/CTCF_Q4/10xplot/*read1_shuffle_phase_1.out F3/A
cp $LIBRARY/10phase/CTCF_Q4/*original_250bp_phase_sort_merge.svg F3/A
cp $LIBRARY/10phase/CTCF_Q4/*shuffle_250bp_phase_sort_merge.svg  F3/A

# ===============================================================================================================================

[ -d F3/B ] || mkdir F3/B
cp $LIBRARY/10phase/CTCF_phase_aligned/*Q4_phase_prefered.out F3/B
cp $LIBRARY/10phase/CTCF_phase_aligned/*Q1_phase_prefered.out F3/B




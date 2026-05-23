#!/bin/bash

set -exo
source activate bx

LIBRARY=../X_Bulk_Processing/Library
SCRIPTMANAGER=../bin/ScriptManager-v0.15.jar

[ -d F3 ] || mkdir F3

# ===============================================================================================================================

[ -d F3/A ] || mkdir F3/A

# Heatmaps
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

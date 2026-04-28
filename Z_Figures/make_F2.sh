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

[ -d F2 ] || mkdir F2

[ -d F2/E ] || mkdir F2/E

# ===============================================================================================================================

cp $LIBRARY/WebLogos/CTCF_MA1929.1_logo.eps F2/E
cp $LIBRARY/WebLogos/ZKSCAN1_MA1585.1_logo.eps F2/E
cp $LIBRARY/WebLogos/REST_MA0138.3_logo.eps F2/E
cp $LIBRARY/WebLogos/ATF7_MA0834.1_logo.eps F2/E

# Composites
cp $LIBRARY/CTCF_Q*_1000bp/Composites/*.out F2/E
cp $LIBRARY/BI_Pileups/MA1929/Composites/*.out F2/E
cp $LIBRARY/BI_Pileups/MA1585/Composites/*.out F2/E
cp $LIBRARY/BI_Pileups/MA0138/Composites/*.out F2/E
cp $LIBRARY/BI_Pileups/MA0834/Composites/*.out F2/E


[ -d F2B ] || mkdir F2/B
# ===============================================================================================================================

cp $LIBRARY/CTCF_Q4_500bp/SVG/*.svg  F2/B
cp $LIBRARY/CTCF_Q4_1000bp/Composites/*.out F2/B


[ -d F2/C ] || mkdir F2/C
# ===============================================================================================================================

cp $LIBRARY/CTCF_Q1-ClosestDyad_500bp/SVG/*.svg  F2/C
cp $LIBRARY/CTCF_Q1_1000bp/Composites/*.out F2/C


[ -d F2/D ] || mkdir F2/D
# ===============================================================================================================================
cp $LIBRARY/10phase/CTCF_Q1/10xplot/*read1_original_phase_3.out F2/D
cp $LIBRARY/10phase/CTCF_Q1/10xplot/*read1_original_phase_8.out F2/D
cp $LIBRARY/10phase/CTCF_Q1/10xplot/*read1_shuffle_phase_3.out F2/D
cp $LIBRARY/10phase/CTCF_Q1/10xplot/*read1_shuffle_phase_8.out F2/D
cp $LIBRARY/10phase/CTCF_Q1/*original_250bp_phase_sort_merge.svg F2/D
cp $LIBRARY/10phase/CTCF_Q1/*shuffle_250bp_phase_sort_merge.svg  F2/D


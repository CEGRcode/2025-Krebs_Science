#!/bin/bash

# Organize data from X_Bulk_Processing into Z_Figures for S3c

### CHANGE ME
WRK=/storage/group/bfp2/default/hxc585_HainingChen/2025_Chen_TF-Nuc/Z_Figures
#WRK=/storage/home/owl5022/scratch/2024-Krebs_Science/Z_Figures
###

set -exo
module load anaconda3
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

LIBRARY=$WRK/../X_Bulk_Processing/Library
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar

[ -d S3 ] || mkdir S3

# ===============================================================================================================================

[ -d S3/A ] || mkdir S3/A

# Composites
cp $LIBRARY/10phase/CTCF_Q4/10xplot/*_-1Nuc_read1_original_phase_*.out S3/A/
cp $LIBRARY/10phase/CTCF_Q4/10xplot/*_+1Nuc_read1_original_phase_*.out S3/A/
cp $LIBRARY/10phase/CTCF_Q4/10xplot/RR*.out S3/A/
cp $LIBRARY/10phase/CTCF_Q4/10xplot/YY*.out S3/A/
cp $LIBRARY/10phase/CTCF_Q4/10xplot/WW*.out S3/A/
cp $LIBRARY/10phase/CTCF_Q4/10xplot/SS*.out S3/A/



[ -d S3/B ] || mkdir S3/B
cp $LIBRARY/DNAshape/ATF7_MA0834.1/Composites/*.out S3/B
cp $LIBRARY/DNAshape/CTCF_MA1929.1/Composites/*.out S3/B
cp $LIBRARY/DNAshape/ZKSCAN1_MA1585.1/Composites/*.out S3/B
cp $LIBRARY/DNAshape/ESRRA_MA0592.1/Composites/*.out S3/B
cp $LIBRARY/DNAshape/JUND_MA0492.2/Composites/*.out S3/B
cp $LIBRARY/DNAshape/NRF1_MA0506.3/Composites/*.out S3/B
cp $LIBRARY/DNAshape/REST_MA0138.3/Composites/*.out S3/B
cp $LIBRARY/DNAshape/SP1_MA0079.5/Composites/*.out S3/B
cp $LIBRARY/DNAshape/SPI1_MA0080.7/Composites/*.out S3/B
cp $LIBRARY/DNAshape/ZNF263_MA0528.1/Composites/*.out S3/B

[ -d S3/C ] || mkdir S3/C

cp $LIBRARY/10phase/Random_Q1/*.svg S3/C/
cp $LIBRARY/10phase/Random_phase_aligned/*.out S3/C/





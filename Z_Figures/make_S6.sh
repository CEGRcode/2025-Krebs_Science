#!/bin/bash

set -exo
source activate bx

LIBRARY=../X_Bulk_Processing/Library
SCRIPTMANAGER=../bin/ScriptManager-v0.15.jar

[ -d S6 ] || mkdir S6

# ===============================================================================================================================

[ -d S6/A ] || mkdir S6/A

cp $LIBRARY/NFIA_SORT-Occupancy_GROUP-Q4/Composites/BNase-seq_50U-10min_merge_hg38_NFIA_SORT-Occupancy_GROUP-Q4_5read1_TotalTag.out S6/A/

# ===============================================================================================================================

[ -d S6/B ] || mkdir S6/B

# Heatmaps
cp $LIBRARY/10phase/NFIA_SORT-Occupancy_GROUP-Q4/BNase-seq_50U-10min_merge_hg38_NFIA_SORT-Occupancy_GROUP-Q4_nearNuc_read1_shuffle_250bp_phase_sort_merge.svg S6/B/
cp $LIBRARY/10phase/NFIA_SORT-Occupancy_GROUP-Q4/BNase-seq_50U-10min_merge_hg38_NFIA_SORT-Occupancy_GROUP-Q4_nearNuc_read1_original_250bp_phase_sort_merge.svg S6/B/

# Composites
cp $LIBRARY/10phase/NFIA_SORT-Occupancy_GROUP-Q4/10xplot/BNase-seq_50U-10min_merge_hg38_NFIA_SORT-Occupancy_GROUP-Q4_nearNuc_read1_original_phase_4.out S6/B/
cp $LIBRARY/10phase/NFIA_SORT-Occupancy_GROUP-Q4/10xplot/BNase-seq_50U-10min_merge_hg38_NFIA_SORT-Occupancy_GROUP-Q4_nearNuc_read1_shuffle_phase_4.out S6/B/
cp $LIBRARY/10phase/NFIA_SORT-Occupancy_GROUP-Q4/10xplot/BNase-seq_50U-10min_merge_hg38_NFIA_SORT-Occupancy_GROUP-Q4_nearNuc_read1_original_phase_9.out S6/B/
cp $LIBRARY/10phase/NFIA_SORT-Occupancy_GROUP-Q4/10xplot/BNase-seq_50U-10min_merge_hg38_NFIA_SORT-Occupancy_GROUP-Q4_nearNuc_read1_shuffle_phase_9.out S6/B/
cp $LIBRARY/10phase/NFIA_SORT-Occupancy_GROUP-Q4/10xplot/WW_NFIA_nearNuc_original_phase_9.out S6/B/
cp $LIBRARY/10phase/NFIA_SORT-Occupancy_GROUP-Q4/10xplot/RR_NFIA_nearNuc_original_phase_9.out S6/B/
cp $LIBRARY/10phase/NFIA_SORT-Occupancy_GROUP-Q4/10xplot/YY_NFIA_nearNuc_original_phase_9.out S6/B/
cp $LIBRARY/10phase/NFIA_SORT-Occupancy_GROUP-Q4/10xplot/SS_NFIA_nearNuc_original_phase_9.out S6/B/

#!/bin/bash

set -exo
source activate bx

LIBRARY=../X_Bulk_Processing/Library
SCRIPTMANAGER=../bin/ScriptManager-v0.15.jar
BED=CTCF_SORT-Occupancy_1000bp
BED1=CTCF_SORT-SMC3Engagement_1000
BED2=CTCF_SORT-SMC3Engagement_GROUP-High_1000bp
BED3=CTCF_SORT-SMC3Engagement_GROUP-Low_1000bp

[ -d S4 ] || mkdir S4

# ===============================================================================================================================

[ -d S4/A ] || mkdir S4/A

# Composites
cp $LIBRARY/$BED/Composites/K562_CTCF_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS.out S4/A/
## tag pile up by nonxl CTCF
java -jar $SCRIPTMANAGER read-analysis tag-pileup ../data/RefPT-Motif/1000bp/CTCF_SORT-Occupancy_1000bp.bed ../data/BAM/K562_CTCF_NonXLBX_rep1_hg38.bam -n 100 --cpu 4 -o S4/A/K562_CTCF_NonXLBX_rep1_hg38_CTCF_SORT-Occupancy_1000bp_5read1-MIN100.out

# ===============================================================================================================================

[ -d S4/B ] || mkdir S4/B

# Heatmaps
cp $LIBRARY/$BED1/SVG/K562_CTCF_BX_rep1_hg38_${BED1}_5read1-MIN100_NCIS_merge_label.svg S4/B/
cp $LIBRARY/$BED1/SVG/K562_RAD21_BX_rep1_hg38_${BED1}_5read1-MIN100_NCIS_merge_label.svg S4/B/
cp $LIBRARY/$BED1/SVG/K562_SMC3_BX_rep1_hg38_${BED1}_5read1-MIN100_NCIS_merge_label.svg S4/B/

# Composites
cp $LIBRARY/$BED2/Composites/K562_CTCF_BX_rep1_hg38_${BED2}_5read1-MIN100_NCIS.out S4/B/
cp $LIBRARY/$BED2/Composites/K562_RAD21_BX_rep1_hg38_${BED2}_5read1-MIN100_NCIS.out S4/B/
cp $LIBRARY/$BED2/Composites/K562_SMC3_BX_rep1_hg38_${BED2}_5read1-MIN100_NCIS.out S4/B/
cp $LIBRARY/$BED3/Composites/K562_CTCF_BX_rep1_hg38_${BED3}_5read1-MIN100_NCIS.out S4/B/
cp $LIBRARY/$BED3/Composites/K562_RAD21_BX_rep1_hg38_${BED3}_5read1-MIN100_NCIS.out S4/B/
cp $LIBRARY/$BED3/Composites/K562_SMC3_BX_rep1_hg38_${BED3}_5read1-MIN100_NCIS.out S4/B/

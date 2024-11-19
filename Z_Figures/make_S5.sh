#!/bin/bash

# Organize data from X_Bulk_Processing into Z_Figures for S5

### CHANGE ME
WRK=/path/to/2024-Krebs_Science/Z_Figures
WRK=/storage/home/owl5022/scratch/2024-Krebs_Science/Z_Figures
###

set -exo
module load anaconda3
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

LIBRARY=$WRK/../X_Bulk_Processing/Library
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar

[ -d S5 ] || mkdir S5

# ===============================================================================================================================

[ -d S5/a ] || mkdir S5/a

# Composites
BED=FOXA_LABEL-K562_SORT-ClosestDyad_GROUP-NFR_1000bp
cp $LIBRARY/$BED/Composites/K562_FOXA2_BX_rep1_hg38_${BED}_5read1_NCIS.out S5/a/
cp $LIBRARY/$BED/Composites/K562_FOXA2_BX_rep1_hg38_${BED}_5read2_NCIS.out S5/a/
BED=FOXA_LABEL-K562_SORT-ClosestDyad_GROUP-NoOverlap_1000bp
cp $LIBRARY/$BED/Composites/K562_FOXA2_BX_rep1_hg38_${BED}_5read1_NCIS.out S5/a/
cp $LIBRARY/$BED/Composites/K562_FOXA2_BX_rep1_hg38_${BED}_5read2_NCIS.out S5/a/


# ===============================================================================================================================

[ -d S5/b ] || mkdir S5/b

# Heatmaps
BED=FOXA_LABEL-uHepG2_SORT-NucleosomeEngagement_500bp
cp $LIBRARY/$BED/SVG/HepG2_FOXA1_BX_rep1_hg38_${BED}_5read1_NCIS_merge_label.svg S5/b/
cp $LIBRARY/$BED/SVG/HepG2_FOXA1_BX_rep1_hg38_${BED}_5read2_NCIS_merge_label.svg S5/b/

# Composites
BED=FOXA_LABEL-uHepG2_SORT-NucleosomeEngagement_GROUP-Engaged_500bp
cp $LIBRARY/$BED/Composites/HepG2_FOXA1_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS.out S5/b/
cp $LIBRARY/$BED/Composites/HepG2_FOXA1_BX_rep1_hg38_${BED}_5read2-MIN100_NCIS.out S5/b/
BED=FOXA_LABEL-uHepG2_SORT-NucleosomeEngagement_GROUP-LessEngaged_500bp
cp $LIBRARY/$BED/Composites/HepG2_FOXA1_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS.out S5/b/
cp $LIBRARY/$BED/Composites/HepG2_FOXA1_BX_rep1_hg38_${BED}_5read2-MIN100_NCIS.out S5/b/

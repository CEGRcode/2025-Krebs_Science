#!/bin/bash

# Organize data from X_Bulk_Processing into Z_Figures for F5

### CHANGE ME
WRK=/path/to/2024-Krebs_Science/Z_Figures
WRK=/storage/home/owl5022/scratch/2024-Krebs_Science/Z_Figures
###

set -exo
module load anaconda3
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

LIBRARY=$WRK/../X_Bulk_Processing/Library
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar

[ -d F5 ] || mkdir F5

# ===============================================================================================================================

[ -d F5/a ] || mkdir F5/a

BED=FOXA_SORT-ClosestDyad_STACK-K562-uHepG2_1000bp

# Sequence figs
cp $LIBRARY/WebLogos/FOXA2_M1_logo.eps F5/a/
cp $LIBRARY/$BED/FourColor/${BED}_31bp.svg F5/a/

# Heatmaps
cp $LIBRARY/$BED/SVG/K562_IgG_BX_merge_hg38_${BED}_5read1_Raw_merge_label.svg F5/a/
cp $LIBRARY/$BED/SVG/K562_FOXA2_BX_rep1_hg38_${BED}_5read1_NCIS_merge_label.svg F5/a/
cp $LIBRARY/$BED/SVG/K562_FOXA1_BX_rep1_hg38_${BED}_5read1_NCIS_merge_label.svg F5/a/
cp $LIBRARY/$BED/SVG/HepG2_IgG_BX_merge_hg38_${BED}_5read1_Raw_merge_label.svg F5/a/
cp $LIBRARY/$BED/SVG/HepG2_FOXA1_BX_rep1_hg38_${BED}_5read1_NCIS_merge_label.svg F5/a/
cp $LIBRARY/$BED/SVG/HepG2_FOXA2_BX_rep1_hg38_${BED}_5read1_NCIS_merge_label.svg F5/a/
cp $LIBRARY/$BED/SVG/BNase-seq_50U-10min_merge_hg38_${BED}_midpoint_TotalTag_combined.svg F5/a/


# ===============================================================================================================================

[ -d F5/b ] || mkdir F5/b

# Composites
BED=FOXA_LABEL-K562_SORT-ClosestDyad_1000bp
cp $LIBRARY/$BED/Composites/BNase-seq_50U-10min_merge_hg38_${BED}_5read2_TotalTag.out F5/b/
BED=FOXA_LABEL-uHepG2_SORT-ClosestDyad_1000bp
cp $LIBRARY/$BED/Composites/BNase-seq_50U-10min_merge_hg38_${BED}_5read2_TotalTag.out F5/b/

# ===============================================================================================================================

[ -d F5/c ] || mkdir F5/c

# Composites
BED=FOXA_LABEL-K562_SORT-ClosestDyad_GROUP-NFR_1000bp
cp $LIBRARY/$BED/Composites/K562_FOXA1_BX_rep1_hg38_${BED}_5read1_NCIS.out F5/c/
cp $LIBRARY/$BED/Composites/K562_FOXA1_BX_rep1_hg38_${BED}_5read2_NCIS.out F5/c/
BED=FOXA_LABEL-K562_SORT-ClosestDyad_GROUP-NoOverlap_1000bp
cp $LIBRARY/$BED/Composites/K562_FOXA1_BX_rep1_hg38_${BED}_5read1_NCIS.out F5/c/
cp $LIBRARY/$BED/Composites/K562_FOXA1_BX_rep1_hg38_${BED}_5read2_NCIS.out F5/c/


# ===============================================================================================================================

[ -d F5/d ] || mkdir F5/d

# Heatmaps
BED=FOXA_LABEL-K562_SORT-NucleosomeEngagement_500bp
cp $LIBRARY/$BED/SVG/K562_FOXA1_BX_rep1_hg38_${BED}_5read1_NCIS_merge_label.svg F5/d
cp $LIBRARY/$BED/SVG/K562_FOXA1_BX_rep1_hg38_${BED}_5read2_NCIS_merge_label.svg F5/d

# Composites
BED=FOXA_LABEL-K562_SORT-NucleosomeEngagement_GROUP-Engaged_500bp
cp $LIBRARY/$BED/Composites/K562_FOXA1_BX_rep1_hg38_${BED}_5read2-MIN100_NCIS.out F5/d/
cp $LIBRARY/$BED/Composites/K562_FOXA1_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS.out F5/d/
cp $LIBRARY/$BED/Composites/K562_FOXA1_BX_rep1_hg38_${BED}_5read2_NCIS.out F5/d/
cp $LIBRARY/$BED/Composites/K562_FOXA1_BX_rep1_hg38_${BED}_5read1_NCIS.out F5/d/
cp $LIBRARY/$BED/Composites/BNase-seq_50U-10min_merge_hg38_${BED}_5read2_TotalTag.out F5/d/
BED=FOXA_LABEL-K562_SORT-NucleosomeEngagement_GROUP-LessEngaged_500bp
cp $LIBRARY/$BED/Composites/K562_FOXA1_BX_rep1_hg38_${BED}_5read2-MIN100_NCIS.out F5/d/
cp $LIBRARY/$BED/Composites/K562_FOXA1_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS.out F5/d/
cp $LIBRARY/$BED/Composites/K562_FOXA1_BX_rep1_hg38_${BED}_5read1_NCIS.out F5/d/
cp $LIBRARY/$BED/Composites/K562_FOXA1_BX_rep1_hg38_${BED}_5read1_NCIS.out F5/d/
cp $LIBRARY/$BED/Composites/BNase-seq_50U-10min_merge_hg38_${BED}_5read2_TotalTag.out F5/d/

#!/bin/bash

set -exo
source activate bx

LIBRARY=../X_Bulk_Processing/Library
SCRIPTMANAGER=../bin/ScriptManager-v0.15.jar

[ -d F6 ] || mkdir F6

# ===============================================================================================================================

[ -d F6/A ] || mkdir F6/A

cp $LIBRARY/WebLogos/NFIA_M1_logo.eps F6/A/

# Composites
BED=NFIA_SORT-Occupancy_500bp
cp $LIBRARY/$BED/Composites/dbSnp153_snv_${BED}.out F6/A
cp $LIBRARY/$BED/Composites/hg38.phyloP30way_${BED}.out F6/A

# ===============================================================================================================================

[ -d F6/B ] || mkdir F6/B

# Heatmaps
BED=NFIA_SORT-DistClosestDyad_1000bp
cp $LIBRARY/$BED/SVG/BNase-seq_50U-10min_merge_hg38_${BED}_midpoint-MIN100_TotalTag_combined.svg F6/B
cp $LIBRARY/$BED/SVG/K562_NFIA_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS_merge_label.svg F6/B

# ===============================================================================================================================

[ -d F6/C ] || mkdir F6/C

BED=NFIA_SORT-Occupancy_500bp
cp $LIBRARY/$BED/FourColor/${BED}_31bp.svg F6/C/

# Heatmaps
BED=NFIA_SORT-Occupancy_500bp
cp $LIBRARY/$BED/SVG/K562_NFIA_BX_rep1_hg38_${BED}_5read1_NCIS_merge_label.svg F6/C/
cp $LIBRARY/$BED/SVG/K562_IgG_BX_merge_hg38_${BED}_5read1_Raw_merge_label.svg F6/C/

# Fat heatmap (threshold should match Bulk Processing config file)
THRESH=1
BED=NFIA-d250bp_SORT-Occupancy_250bp
CDT=$LIBRARY/$BED/CDT/K562_NFIA_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS
BASE=F6/C/K562_NFIA_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation heatmap --blue -x 600 -y 600 -a $THRESH ${CDT}_sense.cdt -o ${BASE}_sense_treeviewFat.png
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation heatmap --red  -x 600 -y 600 -a $THRESH ${CDT}_anti.cdt -o ${BASE}_anti_treeviewFat.png
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation merge-heatmap ${BASE}_sense_treeviewFat.png  ${BASE}_anti_treeviewFat.png -o ${BASE}_merge_treeviewFat.png
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation label-heatmap ${BASE}_merge_treeviewFat.png \
    -l "-250" -m "0" -r "+250" -w 1 -f 20 \
    -o ${BASE}_merge_treeviewFat_label.svg

# ===============================================================================================================================

[ -d F6/D ] || mkdir F6/D

# Composites
BED=NFIA_SORT-DistClosestDyad_GROUP-Downstream_1000bp
cp $LIBRARY/$BED/Composites/K562_NFIA_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS.out F6/D
cp $LIBRARY/$BED/Composites/K562_NFIA_BX_rep1_hg38_${BED}_5read2-MIN100_NCIS.out F6/D

BED=NFIA_SORT-DistClosestDyad_GROUP-Overlap_1000bp
cp $LIBRARY/$BED/Composites/K562_NFIA_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS.out F6/D
cp $LIBRARY/$BED/Composites/K562_NFIA_BX_rep1_hg38_${BED}_5read2-MIN100_NCIS.out F6/D

BED=NFIA_SORT-DistClosestDyad_GROUP-Upstream_1000bp
cp $LIBRARY/$BED/Composites/K562_NFIA_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS.out F6/D
cp $LIBRARY/$BED/Composites/K562_NFIA_BX_rep1_hg38_${BED}_5read2-MIN100_NCIS.out F6/D

# ===============================================================================================================================

# Figure 6E is a graphic

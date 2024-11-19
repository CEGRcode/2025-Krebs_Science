#!/bin/bash

# Organize data from X_Bulk_Processing into Z_Figures for F6

### CHANGE ME
WRK=/path/to/2024-Krebs_Science/Z_Figures
WRK=/storage/home/owl5022/scratch/2024-Krebs_Science/Z_Figures
###

set -exo
module load anaconda3
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

LIBRARY=$WRK/../X_Bulk_Processing/Library
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar
VIOLIN=$WRK/../bin/make_violin_plot.py

[ -d F6 ] || mkdir F6

# ===============================================================================================================================

[ -d F6/a ] || mkdir F6/a

cp $LIBRARY/WebLogos/NFIA_M1_logo.eps F6/a/

# Composites
BED=NFIA_SORT-Occupancy_500bp
cp $LIBRARY/$BED/Composites/dbSnp153_snv_${BED}.out F6/a
cp $LIBRARY/$BED/Composites/hg38.phyloP30way_${BED}.out F6/a

# ===============================================================================================================================

[ -d F6/b ] || mkdir F6/b

BED=NFIA_SORT-Occupancy_500bp
cp $LIBRARY/$BED/FourColor/${BED}_31bp.svg F6/b/

# Heatmaps
BED=NFIA_SORT-Occupancy_500bp
cp $LIBRARY/$BED/SVG/K562_NFIA_BX_rep1_hg38_${BED}_5read1_NCIS_merge_label.svg F6/b/
cp $LIBRARY/$BED/SVG/K562_IgG_BX_merge_hg38_${BED}_5read1_Raw_merge_label.svg F6/b/

# Fat heatmap (threshold should match Bulk Processing config file)
THRESH=1
BED=NFIA-d250bp_SORT-Occupancy_250bp
CDT=$LIBRARY/$BED/CDT/K562_NFIA_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS
BASE=F6/b/K562_NFIA_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation heatmap --blue -x 600 -y 600 -a $THRESH ${CDT}_sense.cdt -o ${BASE}_sense_treeviewFat.png
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation heatmap --red  -x 600 -y 600 -a $THRESH ${CDT}_anti.cdt -o ${BASE}_anti_treeviewFat.png
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation merge-heatmap ${BASE}_sense_treeviewFat.png  ${BASE}_anti_treeviewFat.png -o ${BASE}_merge_treeviewFat.png
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation label-heatmap ${BASE}_merge_treeviewFat.png \
    -l "-250" -m "0" -r "+250" -w 1 -f 20 \
    -o ${BASE}_merge_treeviewFat_label.svg

# ===============================================================================================================================

[ -d F6/c ] || mkdir F6/c

# Heatmaps
BED=NFIA_SORT-DistClosestDyad_1000bp
cp $LIBRARY/$BED/SVG/BNase-seq_50U-10min_merge_hg38_${BED}_midpoint-MIN100_TotalTag_combined.svg F6/c
cp $LIBRARY/$BED/SVG/K562_NFIA_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS_merge_label.svg F6/c

# Composites
BED=NFIA_SORT-DistClosestDyad_GROUP-Downstream_1000bp
cp $LIBRARY/$BED/Composites/K562_NFIA_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS.out F6/c
cp $LIBRARY/$BED/Composites/K562_NFIA_BX_rep1_hg38_${BED}_5read2-MIN100_NCIS.out F6/c

BED=NFIA_SORT-DistClosestDyad_GROUP-Overlap_1000bp
cp $LIBRARY/$BED/Composites/K562_NFIA_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS.out F6/c
cp $LIBRARY/$BED/Composites/K562_NFIA_BX_rep1_hg38_${BED}_5read2-MIN100_NCIS.out F6/c

BED=NFIA_SORT-DistClosestDyad_GROUP-Upstream_1000bp
cp $LIBRARY/$BED/Composites/K562_NFIA_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS.out F6/c
cp $LIBRARY/$BED/Composites/K562_NFIA_BX_rep1_hg38_${BED}_5read2-MIN100_NCIS.out F6/c

# ===============================================================================================================================

[ -d F6/d ] || mkdir F6/d

DATAFILE=F6/d/ViolinData.txt

# Aggregate data file for violin
BED=NFIA-u95_SORT-Occupancy_150bp
BASE=$LIBRARY/$BED/CDT/K562_NFIA_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m $BASE\_sense.cdt $BASE\_anti.cdt -o F6/d/Oriented_Upstream.out
tail -n +2 F6/d/Oriented_Upstream.out | awk 'BEGIN{OFS="\t";FS="\t"}{print $2+$3, "Oriented_Upstream"}' > $DATAFILE

BED=NFIA-d95_SORT-Occupancy_150bp
BASE=$LIBRARY/$BED/CDT/K562_NFIA_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m $BASE\_sense.cdt $BASE\_anti.cdt -o F6/d/Oriented_Downstream.out
tail -n +2 F6/d/Oriented_Downstream.out | awk 'BEGIN{OFS="\t";FS="\t"}{print $2+$3, "Oriented_Downstream"}' >> $DATAFILE

BED=NFIA-u95_REORIENT-Random_SORT-Occupancy_150bp
BASE=$LIBRARY/$BED/CDT/K562_NFIA_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m $BASE\_sense.cdt $BASE\_anti.cdt -o F6/d/Random_Upstream.out
tail -n +2 F6/d/Random_Upstream.out | awk 'BEGIN{OFS="\t";FS="\t"}{print $2+$3, "Random_Upstream"}' >> $DATAFILE

BED=NFIA-d95_REORIENT-Random_SORT-Occupancy_150bp
BASE=$LIBRARY/$BED/CDT/K562_NFIA_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m $BASE\_sense.cdt $BASE\_anti.cdt -o F6/d/Random_Downstream.out
tail -n +2 F6/d/Random_Downstream.out | awk 'BEGIN{OFS="\t";FS="\t"}{print $2+$3, "Random_Downstream"}' >> $DATAFILE

# Compress for storage savings
gzip -f $DATAFILE

# Plot violin data
python $VIOLIN -i $DATAFILE.gz -o F6/d/ViolinData.svg

# Clean-up
rm F6/d/*.out
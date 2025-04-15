#!/bin/bash

# Organize data from X_Bulk_Processing into Z_Figures for F6

### CHANGE ME
WRK=/storage/group/bfp2/default/hxc585_HainingChen/2025_Chen_TF-Nuc/Z_Figures
#WRK=/storage/home/owl5022/scratch/2024-Krebs_Science/Z_Figures
###

set -exo
module load anaconda3
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

LIBRARY=$WRK/../X_Bulk_Processing/Library
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar
VIOLIN=$WRK/../bin/NFIA_violin.py

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

BED=NFIA-d250bp_SORT-Occupancy_250bp
mv $LIBRARY/$BED/SVG/K562_NFIA_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS_merge_label.svg  F6/b/

[ -d F6/c ] || mkdir F6/c

# Heatmaps
BED=NFIA_SORT-DistClosestDyad_1000bp
cp $LIBRARY/$BED/SVG/BNase-seq_50U-10min_merge_hg38_${BED}_midpoint_TotalTag_combined.svg F6/c
cp $LIBRARY/$BED/SVG/K562_NFIA_BX_rep1_hg38_${BED}_5read1_NCIS_merge_label.svg F6/c

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

DATAFILE=F6/d/ViolinData_original.txt

BED=NFIA-u95_SORT-Occupancy_150bp
BASE=$LIBRARY/$BED/CDT/K562_NFIA_BX_rep1_hg38_${BED}_5read1-MIN100
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m $BASE\_sense.cdt $BASE\_anti.cdt -o F6/d/Oriented_Upstream.out
tail -n +2 F6/d/Oriented_Upstream.out | awk '{OFS="\t"} {print $1,"oriented_upstream", $2+$3}'  > $DATAFILE

BED=NFIA-d95_SORT-Occupancy_150bp
BASE=$LIBRARY/$BED/CDT/K562_NFIA_BX_rep1_hg38_${BED}_5read1-MIN100
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m $BASE\_sense.cdt $BASE\_anti.cdt -o F6/d/Oriented_Downstream.out
tail -n +2 F6/d/Oriented_Downstream.out | awk '{OFS="\t"} {print $1,"oriented_downstream", $2+$3}'  >> $DATAFILE
awk '{OFS = ","} {print $1, $2, $3}' $DATAFILE > F6/d/temp.csv
awk 'BEGIN {OFS = ","; print "Location", "NFIA-Nuc_engagement"} { 
  if ($3 != 0) {
    print $1, $2, log10($3) 
  }
}' F6/d/temp.csv > F6/d/NFIA_Nuc_Engagement.csv

rm F6/d/temp.csv


DATAFILE2=F6/d/ViolinData_random.txt

BED=NFIA-u95_REORIENT-Random_SORT-Occupancy_150bp
BASE=$LIBRARY/$BED/CDT/K562_NFIA_BX_rep1_hg38_${BED}_5read1-MIN100
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m $BASE\_sense.cdt $BASE\_anti.cdt -o F6/d/Random_Upstream.out
tail -n +2 F6/d/Random_Upstream.out | awk '{OFS="\t"} {print $1, "random_upstream", $2+$3}'  >> $DATAFILE2

BED=NFIA-d95_REORIENT-Random_SORT-Occupancy_150bp
BASE=$LIBRARY/$BED/CDT/K562_NFIA_BX_rep1_hg38_${BED}_5read1-MIN100
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m $BASE\_sense.cdt $BASE\_anti.cdt -o F6/d/Random_Downstream.out
tail -n +2 F6/d/Random_Downstream.out | awk '{OFS="\t"} {print $1,"random_downstream", $2+$3}' >> $DATAFILE2

awk '{OFS = ","} {print $1, $2, $3}' $DATAFILE2 > F6/d/temp.csv
awk 'BEGIN {OFS = ","; print "Location", "NFIA-Nuc_engagement"} { 
  if ($3 != 0) {
    print $1, $2, log10($3) 
  }
}' F6/d/temp.csv > F6/d/Random_NFIA_Nuc_Engagement.csv

rm F6/d/temp.csv
# Plot violin data
python $VIOLIN -i F6/d/NFIA_Nuc_Engagement.csv -o F6/d/ViolinData1.svg
python $VIOLIN -i F6/d/Random_NFIA_Nuc_Engagement.csv -o F6/d/ViolinData2.svg
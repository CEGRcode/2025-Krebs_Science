#!/bin/bash

# NFIA plots need custom color thresholding scale for visualization and IgG scales should be matching.

### CHANGE ME
WRK=/path/to/2024-Chen_Nature/Z_Figures
WRK=/storage/home/owl5022/scratch/2024-Chen_Nature/Z_Figures
###

# Dependencies
# - java
# - samtools

set -exo
module load samtools

# Fill in placeholder constants with your directories
LIBRARY=../0X_Bulk_Processing/Library
BAMDIR=$WRK/../data/BAM
MOTIF=$WRK/../data/RefPT-Motif

# Setup ScriptManager for job array
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.14.jar

# Set up output directories
[ -d logs ] || mkdir logs
[ -d F2/b ] || mkdir -p F2/b

# NFIA pileup on NFIA_Occupancy_1000bp
CDTBASE="$LIBRARY/NFIA_Occupancy_1000bp/CDT/K562_NFIA_BX_rep1_hg19_NFIA_Occupancy_1000bp_read1";
CDT=`basename $CDTBASE`
NSITES=`wc -l $CDTBASE\_sense.cdt | awk '{print $1-1}'`

# Heatmap match and mismatch (different colors)
java -jar $SCRIPTMANAGER figure-generation heatmap -a 1 --blue $CDTBASE\_sense.cdt -o SENSE.png
java -jar $SCRIPTMANAGER figure-generation heatmap -a 1 --red  $CDTBASE\_anti.cdt  -o ANTI.png
java -jar $SCRIPTMANAGER figure-generation merge-heatmap SENSE.png ANTI.png -o $CDT\_merge.png
java -jar $SCRIPTMANAGER figure-generation label-heatmap $CDT\_merge.png \
    -l "-500" -m "0" -r "+500" -w 1 -f 20 \
    -x "NFIA_Occupancy_1000bp" -y "NFIA_Occupancy_1000bp (${NSITES} sites)" \
    -o F2/b/$CDT\_merge_label.svg

# IgG pileup on NFIA_Occupancy_1000bp
CDTBASE="$LIBRARY/NFIA_Occupancy_1000bp/CDT/K562_IgG_BX_merge_hg19_NFIA_Occupancy_1000bp_read1";
CDT=`basename $CDTBASE`
NSITES=`wc -l $CDTBASE\_sense.cdt | awk '{print $1-1}'`

# Heatmap match and mismatch (different colors)
java -jar $SCRIPTMANAGER figure-generation heatmap -a 1 --blue $CDTBASE\_sense.cdt -o SENSE.png
java -jar $SCRIPTMANAGER figure-generation heatmap -a 1 --red  $CDTBASE\_anti.cdt  -o ANTI.png
java -jar $SCRIPTMANAGER figure-generation merge-heatmap SENSE.png ANTI.png -o $CDT\_merge.png
java -jar $SCRIPTMANAGER figure-generation label-heatmap $CDT\_merge.png \
    -l "-500" -m "0" -r "+500" -w 1 -f 20 \
    -x "NFIA_Occupancy_1000bp" -y "NFIA_Occupancy_1000bp (${NSITES} sites)" \
    -o F2/b/$CDT\_merge_label.svg


[ -d F2/f ] || mkdir F2/f

# NFIA pileup on NFIA_NucSort_1000bp
CDTBASE="$LIBRARY/NFIA_NucSort_1000bp/CDT/K562_NFIA_BX_rep1_hg19_NFIA_NucSort_1000bp_read1";
CDT=`basename $CDTBASE`
NSITES=`wc -l $CDTBASE\_sense.cdt | awk '{print $1-1}'`

# Heatmap match and mismatch (different colors)
java -jar $SCRIPTMANAGER figure-generation heatmap -a 1 --blue $CDTBASE\_sense.cdt -o SENSE.png
java -jar $SCRIPTMANAGER figure-generation heatmap -a 1 --red  $CDTBASE\_anti.cdt  -o ANTI.png
java -jar $SCRIPTMANAGER figure-generation merge-heatmap SENSE.png ANTI.png -o $CDT\_merge.png
java -jar $SCRIPTMANAGER figure-generation label-heatmap $CDT\_merge.png \
    -l "-500" -m "0" -r "+500" -w 1 -f 20 \
    -x "NFIA_NucSort_1000bp" -y "NFIA_NucSort_1000bp (${NSITES} sites)" \
    -o F2/f/$CDT\_merge_label.svg

# Clean-up
rm SENSE.png ANTI.png

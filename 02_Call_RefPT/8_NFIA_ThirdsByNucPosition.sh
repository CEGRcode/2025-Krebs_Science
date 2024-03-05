#!/bin/bash

# Break NFIA sorted by nucleosome position into UPSTREAM, OVERLAP, and DOWNSTREAM distance filter groups

### CHANGE ME
WRK=/path/to/2024-Chen_Nature/02_Call_RefPT
###

# Dependencies
# - java
# - awk

set -exo

# Inputs and outputs
MOTIF=$WRK/../data/RefPT-Motif
BEDFILE=$MOTIF/1bp/NFIA_NucSort_1bp.bed

[ -d $MOTIF/100bp ] || mkdir $MOTIF/100bp
[ -d $MOTIF/1000bp ] || mkdir $MOTIF/1000bp

# Script shortcuts
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.14.jar

# Split by sorted distance
awk '{
  if ($13 > 73) {
    print $0 > "downstream.bed";
  } else if ($13 < -73) {
    print $0 > "upstream.bed";
  } else {
    print $0 > "overlap.bed";
  }
}' $BEDFILE

# Expand by 1000bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 overlap.bed -o $MOTIF/1000bp/NFIA_NucSort-OVERLAP_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 upstream.bed -o $MOTIF/1000bp/NFIA_NucSort-UPSTREAM_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 downstream.bed -o $MOTIF/1000bp/NFIA_NucSort-DOWNSTREAM_1000bp.bed

# Expand by 100bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 overlap.bed -o $MOTIF/100bp/NFIA_NucSort-OVERLAP_100bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 upstream.bed -o $MOTIF/100bp/NFIA_NucSort-UPSTREAM_100bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 downstream.bed -o $MOTIF/100bp/NFIA_NucSort-DOWNSTREAM_100bp.bed

# Rename three files
mv downstream.bed $MOTIF/1bp/NFIA_NucSort-DOWNSTREAM_1bp.bed
mv upstream.bed $MOTIF/1bp/NFIA_NucSort-UPSTREAM_1bp.bed
mv overlap.bed $MOTIF/1bp/NFIA_NucSort-OVERLAP_1bp.bed

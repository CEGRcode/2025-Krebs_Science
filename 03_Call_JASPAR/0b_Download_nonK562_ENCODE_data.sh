#!/bin/bash

# Download JASPAR PWMs (.meme) and ENCODE narrowPeak (.bed) files for building
# "lowly bound" TF Motif RefPTs. Some of these downloads will not be used in
# larger analysis so they have been omitted from being saved in the global
# `../data/` directory.

### CHANGE ME
WRK=WRK=/storage/group/bfp2/default/hxc585_HainingChen/2025_Chen_TF-Nuc/03_Call_JASPAR
#WRK=/ocean/projects/see180003p/owlang/2024-Krebs_Science/03_Call_JASPAR
#WRK=/scratch/owl5022/2024-Krebs_Science/03_Call_JASPAR
METADATA=nonK562TF_JASPAR_ENCODE_config.txt
###

# Dependencies
# - wget

JDIR=../data/JASPAR

set -exo
[ -d $JDIR ] || mkdir $JDIR
[ -d NonK562_narrowPeaks ] || mkdir NonK562_narrowPeaks

while read line; do
	# Format input variables
	TF=`echo $line | awk '{print $1}'`
	Cell=`echo $line | awk '{print $2}'`
	ENCFF=`echo $line | awk '{print $3}'`

	# Download ENCODE peaks
	wget -c -O NonK562_narrowPeaks/$TF\_$Cell\_$ENCFF.bed.gz https://www.encodeproject.org/files/$ENCFF/@@download/$ENCFF.bed.gz

done < $METADATA

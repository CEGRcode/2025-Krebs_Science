#!/bin/bash

# Download JASPAR PWMs (.meme) and ENCODE narrowPeak (.bed) files for building
# "lowly bound" TF Motif RefPTs. Some of these downloads will not be used in
# larger analysis so they have been omitted from being saved in the global
# `../data/` directory.

METADATA=nonK562TF_JASPAR_ENCODE_config.txt

# Dependencies
# - wget

set -exo
source activate bx

# Inputs and outputs
JDIR=../data/JASPAR

[ -d $JDIR ] || mkdir $JDIR
[ -d NonK562_narrowPeaks ] || mkdir NonK562_narrowPeaks

while read line; do
	# Format input variables
	TF=`echo $line | awk '{print $1}'`
	Cell=`echo $line | awk '{print $2}'`
	ENCFF=`echo $line | awk '{print $3}'`

	# Download ENCODE peaks
	wget -c -O NonK562_narrowPeaks/${TF}_${Cell}_$ENCFF.bed.gz https://www.encodeproject.org/files/$ENCFF/@@download/$ENCFF.bed.gz

done < $METADATA

#!/bin/bash

# Get TOP, MIDDLE, and BOTTOM thirds of CTCF sites by occupancy

### CHANGE ME
WRK=/path/to/2023-Chen_Benzonase-ChIP-exo/02_Call_RefPT
###

set -exo

# Inputs and outputs
MOTIF=$WRK/../data/RefPT-Motif/1000bp
REFPT=$MOTIF/CTCF_Occupancy_1000bp.bed

# Get line counts and calculate one third
NSITES=`wc -l $REFPT | awk '{print $1}'`
NTHIRD=$(($NSITES/3))
echo "${NSITES} / 3 = ${NTHIRD}"

# Get top N CTCF sites (TOP)
head -n $NTHIRD $REFPT > $MOTIF/CTCF_Occupancy-TOP_1000bp.bed

# Get next N CTCF sites (MIDDLE)
head -n $((2 * $NTHIRD)) $REFPT | tail -n $NTHIRD > $MOTIF/CTCF_Occupancy-MIDDLE_1000bp.bed

# Get remaining third CTCF sites (BOTTOM)
tail -n $(($NSITES - 2*$NTHIRD )) $REFPT > $MOTIF/CTCF_Occupancy-BOTTOM_1000bp.bed
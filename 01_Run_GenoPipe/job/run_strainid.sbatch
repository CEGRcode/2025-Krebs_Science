#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 01:00:00
#SBATCH -A open
#SBATCH -o logs/run_strainid.log.out
#SBATCH -e logs/run_strainid.log.err
#SBATCH --array 1-20

module load samtools
module load anaconda3
source activate genopipe # load conda environment with GenoPipe dependencies

### CHANGE ME

WRK=/storage/group/bfp2/default/hxc585_HainingChen/Fox_NFIA_CTCF/
STRAINID=$WRK/bin/StrainID
GENOME=$WRK/hg38_files/hg38.fa
###

INPUTBAM=$WRK/data/BAM
OUTPUT=$WRK/01_Run_GenoPipe/strainid

cd $WRK
[ -d logs ] || mkdir logs
[ -d $OUTPUT ] || mkdir $OUTPUT

cd $STRAINID
bash identify-Strain.sh -i $INPUTBAM -v $STRAINID/hg38_VCF/ -g $GENOME -o $OUTPUT

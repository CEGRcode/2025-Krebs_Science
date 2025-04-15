#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=24gb
#SBATCH -t 00:40:00
#SBATCH -A open
#SBATCH -o logs/2_exo_read_pileups.log.out-%a
#SBATCH -e logs/2_exo_read_pileups.log.err-%a
#SBATCH --array 1-33

# Pileup read1 (exo cut sites) for a custom combinations of BAM x BED files.
# Heatmaps included and all with scaling. 

### CHANGE ME
WRK=/storage/group/bfp2/default/hxc585_HainingChen/2025_Chen_TF-Nuc/X_Bulk_Processing
METADATA=$WRK/10bp_region.txt
GINFO=$WRK/../data/hg38_files/hg38.chrom.sizes
GENOME=$WRK/../data/hg38_files/hg38.fa
###
module load anaconda
source activate bioinfo
# Dependencies
# - java
# - perl
# - samtools

set -exo
#module load samtools

# Fill in placeholder constants with your directories
BAMDIR=$WRK/data/BAM
OUTDIR=$WRK/Library/10bp

# Setup ScriptManager for job array
#SCRIPTMANAGER=$WRK/../2023_Chen_PIC3/bin/ScriptManager-v0.15.jar
ORIGINAL_SCRIPTMANAGER=../bin/ScriptManager-v0.15.jar
SCRIPTMANAGER=../bin/ScriptManager-v0.15-$SLURM_ARRAY_TASK_ID.jar
cp $ORIGINAL_SCRIPTMANAGER $SCRIPTMANAGER

# Script shortcuts
COMPOSITE=$WRK/../bin/sum_Col_CDT.pl
chisquare=$WRK/../bin/chisquare.py

# Set up output directories
[ -d logs ] || mkdir logs
[ -d $OUTDIR ] || mkdir $OUTDIR

cd $WRK/
# Determine BED file for the current job array index
BEDFILE=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $1}'`
BED=`basename $BEDFILE "_1bp.bed"`

# Determine BAM file for the current job array index
BAMFILE=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $2}'`
BAM=`basename $BAMFILE ".bam"`
#[ -f $BAMFILE.bai ] || samtools index $BAMFILE

# strand
Strand=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $3}'`

# fragment size
Size=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $4}'`

# extract
Extract=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $5}'`

# window
Window=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $6}'`

DIR=$OUTDIR/$BED
[ -d $DIR ] || mkdir $DIR
[[ -d $DIR/CDT ]] || mkdir $DIR/CDT
[[ -d $DIR/SCORES ]] || mkdir $DIR/SCORES

mkdir -p $OUTDIR/BX_10bp


BASE=${BAM}_${BED}_read1-MIN${Size}
echo $BASE

# expand

java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $BEDFILE -o $DIR/${BED}_1000bp.bed

# Pileup (read 1)
java -jar $SCRIPTMANAGER read-analysis tag-pileup $DIR/${BED}_1000bp.bed $BAMFILE --cpu 4 -5 -1 -n ${Size} -M $DIR/CDT/${BASE}

# take sum score of extract region
cat $DIR/CDT/${BASE}_${Strand}.cdt | cut -f 1-2,$Extract > $DIR/CDT/${BASE}_${Extract}_${Strand}.cdt
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum $DIR/CDT/${BASE}_${Extract}_${Strand}.cdt -o $DIR/SCORES/${BASE}_${Extract}_${Strand}_SCORES.out
tail -n +2 $DIR/SCORES/${BASE}_${Extract}_${Strand}_SCORES.out > $DIR/${BASE}_${Extract}_${Strand}_score.out
# calculate ratio of each phase
cut -f 3- $DIR/CDT/${BASE}_${Extract}_${Strand}.cdt | tail -n +2  | \
        awk -v max_col=${Window} '{
            OFS="\t";
            # Initialize sums for all 10 columns
            sum1 = 0; sum2 = 0; sum3 = 0; sum4 = 0; sum5 = 0;
            sum6 = 0; sum7 = 0; sum8 = 0; sum9 = 0; sum10 = 0;

            # Sum the columns dynamically in steps of 10
            for (i = 1; i <= max_col; i += 10) sum1 += $i;
            for (i = 2; i <= (max_col + 1); i += 10) sum2 += $i;
            for (i = 3; i <= (max_col + 2); i += 10) sum3 += $i;
            for (i = 4; i <= (max_col + 3); i += 10) sum4 += $i;
            for (i = 5; i <= (max_col + 4); i += 10) sum5 += $i;
            for (i = 6; i <= (max_col + 5); i += 10) sum6 += $i;
            for (i = 7; i <= (max_col + 6); i += 10) sum7 += $i;
            for (i = 8; i <= (max_col + 7); i += 10) sum8 += $i;
            for (i = 9; i <= (max_col + 8); i += 10) sum9 += $i;
            for (i = 10; i <= (max_col + 9); i += 10) sum10 += $i;

            print sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9, sum10;
        }' | paste $DIR/${BASE}_${Extract}_${Strand}_score.out  - |  awk '{
    OFS="\t";
    print $1,$2,$3/($2+1),$4/($2+1),$5/($2+1),$6/($2+1),$7/($2+1),$8/($2+1),$9/($2+1),$10/($2+1),$11/($2+1),$12/($2+1)}' > $DIR/${BASE}_${Extract}_${Strand}_score_peak.out

rm $DIR/${BASE}_${Extract}_${Strand}_score.out 
rm $DIR/${BED}_1000bp.bed


for file in  $DIR/${BASE}_${Extract}_${Strand}_score_peak.out ; do
            filename=$(basename "$file" ".out")
            
            # Loop through phases 0 to 9 to avoid repetitive code
            for phase in {0..9}; do
                # Each phase corresponds to columns 3 to 12
                awk -v phase="$phase" -v filename="$filename" 'BEGIN {OFS=","} {
                    print $1, phase, $((phase+3)) > (filename"_"phase".csv");
                }' $file
            done

            # Create the final CSV with header and concatenate all phase CSVs
            echo -e "Region,Nucleosomephase,enrichment" >  $DIR/${filename}.csv
            cat ${filename}_*.csv >>  $DIR/${filename}.csv

            # Remove individual phase files
            rm ${filename}_*.csv
done
mkdir -p $WRK/Library/BX_10bp
for folder in $OUTDIR/${BED} ; do
    Target=$(basename "$folder")
    output_file="${Target}.out"
    python $chisquare $folder "$output_file"
    mv "$output_file" $WRK/Library/BX_10bp
done



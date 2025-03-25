#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 00:40:00
#SBATCH -A open
#SBATCH -o logs/2_MotifAnalyses.log.out-%a
#SBATCH -e logs/2_MotifAnalyses.log.err-%a
#SBATCH --array 1-34

# Loop through each BEDFILE in data/RefPT-Motif/1000bp and perform motif analyses:
# - Need BI_Pileups.txt file
# - Four color plots (QC)
# - BI pileups
# - DNA shape
# - Dinucleotide

### CHANGE ME
WRK=/storage/group/bfp2/default/hxc585_HainingChen/2025_Chen_TF-Nuc/X_Bulk_Processing
#WRK=/ocean/projects/see180003p/owlang/2024-Krebs_Science/X_Bulk_Processing
#WRK=/storage/group/bfp2/default/owl5022-OliviaLang/2024-Krebs_Science/X_Bulk_Processing
METADATA=$WRK/BI_Pileups.txt
THREADS=4
###

# Dependencies
# - java
# - opencv
# - perl
# - python

set -exo
module load anaconda
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Setup ScriptManager for job array
ORIGINAL_SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15-$SLURM_ARRAY_TASK_ID.jar
cp $ORIGINAL_SCRIPTMANAGER $SCRIPTMANAGER

# Script shortcuts
RESIZE=../bin/resize_png.py
DINUCLEOTIDE=../bin/scan_FASTA_for_motif_as_binary_string.py
COMPOSITE=../bin/sum_Col_CDT.pl

# Inputs and outputs
GENOME=$WRK/../data/hg38_files/hg38.fa
MOTIF=$WRK/../data/RefPT-JASPAR
OUTDIR=$WRK/Library/BI_Pileups

# Set up output directories
[ -d logs ] || mkdir logs
[ -d $OUTDIR ] || mkdir $OUTDIR

BAMFILE=$WRK/../data/BAM/BNase-seq_50U-10min_merge_hg38.bam
BAM=`basename $BAMFILE ".bam"`

# Determine BED file for the current job array index
RefID=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $1}'`

BEDFILE1=../data/RefPT-JASPAR/1000bp/*_${RefID}*GROUP-Quartile1_1000bp.bed
BEDFILE2=../data/RefPT-JASPAR/1000bp/*_${RefID}*GROUP-Quartile2_1000bp.bed
BEDFILE3=../data/RefPT-JASPAR/1000bp/*_${RefID}*GROUP-Quartile3_1000bp.bed
BEDFILE4=../data/RefPT-JASPAR/1000bp/*_${RefID}*GROUP-Quartile4_1000bp.bed
BED1=`basename $BEDFILE1 ".bed"`
BED2=`basename $BEDFILE2 ".bed"`
BED3=`basename $BEDFILE3 ".bed"`
BED4=`basename $BEDFILE4 ".bed"`

loc1=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $2}'`
column1=$((loc1+502))
loc2=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $3}'`
column2=$((loc2+502))
loc3=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $5}'`
column3=$((loc3+502))
loc4=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $6}'`
column4=$((loc4+502))

Region_sense=$((loc2 - loc1))
Region_anti=$((loc4 - loc3))
Rounded_Region_sense=$((Region_sense / 10 * 10)) 
Rounded_Region_anti=$((Region_anti / 10 * 10))
sense_end_column=$((Rounded_Region_sense + 2))
anti_end_column=$((Rounded_Region_anti + 2))

# Count sites
NSITES=`wc -l $BEDFILE | awk '{print $1}'`

echo "(${SLURM_ARRAY_TASK_ID}) ${RefID} "x" ${BAMFILE} "

DIR=$OUTDIR/${RefID}
[ -d $DIR ] || mkdir $DIR
[[ -d $DIR/CDT ]] || mkdir $DIR/CDT
[[ -d $DIR/Composites ]] || mkdir $DIR/Composites


# ===============================================================================================================================

echo "Run Motif endo cuts BI pileup"
BASE1=$BAM\_$BED1\_5both
BASE2=$BAM\_$BED2\_5both
BASE3=$BAM\_$BED3\_5both
BASE4=$BAM\_$BED4\_5both

# Pileup (endo)
java -jar $SCRIPTMANAGER read-analysis tag-pileup $BEDFILE1 $BAMFILE --cpu $THREADS -5 -a -o $DIR/Composites/$BASE1.out -M $DIR/CDT/$BASE1
java -jar $SCRIPTMANAGER read-analysis tag-pileup $BEDFILE2 $BAMFILE --cpu $THREADS -5 -a -o $DIR/Composites/$BASE2.out -M $DIR/CDT/$BASE2
java -jar $SCRIPTMANAGER read-analysis tag-pileup $BEDFILE3 $BAMFILE --cpu $THREADS -5 -a -o $DIR/Composites/$BASE3.out -M $DIR/CDT/$BASE3
java -jar $SCRIPTMANAGER read-analysis tag-pileup $BEDFILE4 $BAMFILE --cpu $THREADS -5 -a -o $DIR/Composites/$BASE4.out -M $DIR/CDT/$BASE4
# Take region 
cat $DIR/CDT/${BASE1}_sense.cdt | cut -f 1-2,${column1}-${column2}  > $DIR/${RefID}_Q1_sense_final.cdt
cat $DIR/CDT/${BASE1}_anti.cdt | cut -f 1-2,${column3}-${column4}  > $DIR/${RefID}_Q1_anti_final.cdt
cat $DIR/CDT/${BASE2}_sense.cdt | cut -f 1-2,${column1}-${column2}  > $DIR/${RefID}_Q2_sense_final.cdt
cat $DIR/CDT/${BASE2}_anti.cdt | cut -f 1-2,${column3}-${column4}  > $DIR/${RefID}_Q2_anti_final.cdt
cat $DIR/CDT/${BASE3}_sense.cdt | cut -f 1-2,${column1}-${column2}  > $DIR/${RefID}_Q3_sense_final.cdt
cat $DIR/CDT/${BASE3}_anti.cdt | cut -f 1-2,${column3}-${column4}  > $DIR/${RefID}_Q3_anti_final.cdt
cat $DIR/CDT/${BASE4}_sense.cdt | cut -f 1-2,${column1}-${column2}  > $DIR/${RefID}_Q4_sense_final.cdt
cat $DIR/CDT/${BASE4}_anti.cdt | cut -f 1-2,${column3}-${column4}  > $DIR/${RefID}_Q4_anti_final.cdt

  


for file in $DIR/*.cdt ; do
    filename=$(basename "$file" "_final.cdt")
    Ref=`basename $file ".cdt" | cut -d "_" -f 1`
    Sort=`basename $file ".cdt" | cut -d "_" -f 2`
    Strand=`basename $file ".cdt" | cut -d "_" -f 3`
    java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum $DIR/${Ref}_${Sort}_${Strand}_final.cdt -o $DIR/${Ref}_${Sort}_${Strand}_final_SCORES.out
    tail -n +2 $DIR/${Ref}_${Sort}_${Strand}_final_SCORES.out > $DIR/${Ref}_${Sort}_${Strand}_score.out  
    rm $DIR/${Ref}_${Sort}_${Strand}_final_SCORES.out

    if [[ "$Strand" == "sense" ]]; then  
        tail -n +2 ${Ref}_${Sort}_sense_final.cdt | cut -f 3-${sense_end_column} | \
        awk -v max_col=$Rounded_Region_sense '{
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
        }' | paste $DIR/${Ref}_${Sort}_sense_score.out - | \
        awk '{
            OFS="\t";
            print $1, $2, $3/($2+1), $4/($2+1), $5/($2+1), $6/($2+1), $7/($2+1), $8/($2+1), $9/($2+1), $10/($2+1), $11/($2+1), $12/($2+1);
        }' > $DIR/${Ref}_${Sort}_sense_score_peak.out

    elif [[ "$Strand" == "anti" ]]; then 
        tail -n +2 ${Ref}_${Sort}_anti_final.cdt | cut -f 3-${anti_end_column} | \
        awk -v max_col=$Rounded_Region_anti '{
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
        }' | paste $DIR/${Ref}_${Sort}_anti_score.out - | \
        awk '{
            OFS="\t";
            print $1, $2, $3/($2+1), $4/($2+1), $5/($2+1), $6/($2+1), $7/($2+1), $8/($2+1), $9/($2+1), $10/($2+1), $11/($2+1), $12/($2+1);
        }' > $DIR/${Ref}_${Sort}_anti_score_peak.out 
    fi 
    
    rm $DIR/${Ref}_${Sort}_${Strand}_score.out
done


for folder in $OURDIR/$DIR ; do
    if [ -d "$folder" ]; then
        X=$folder
        for file in $X/*_score_peak.out ; do
            filename=$(basename "$file" ".out")
            
            # Loop through phases 0 to 9 to avoid repetitive code
            for phase in {0..9}; do
                # Each phase corresponds to columns 3 to 12
                awk -v phase="$phase" -v filename="$filename" 'BEGIN {OFS=","} {
                    print $1, phase, $((phase+3)) > (filename"_"phase".csv");
                }' $file
            done

            # Create the final CSV with header and concatenate all phase CSVs
            echo -e "Region,Nucleosomephase,enrichment" > $X/${filename}.csv
            cat ${filename}_*.csv >> $X/${filename}.csv

            # Remove individual phase files
            rm ${filename}_*.csv
        done
    fi
done

mkdir -p $WRK/Library/BI_10bp 

for folder in $OURDIR/$DIR ; do
    
    if [ -d "$folder" ]; then
        X=$folder
        output_file="${X}.out"
        python $chisquare "$X" $output_file
		mv $output_file $WRK/Library/BI_10bp
    fi
done
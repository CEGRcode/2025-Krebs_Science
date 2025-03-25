#!/bin/bash

### CHANGE ME
WRK=/storage/group/bfp2/default/hxc585_HainingChen/2025_Chen_TF-Nuc/X_Bulk_Processing
#WRK=/ocean/projects/see180003p/owlang/2024-Krebs_Science/0X_Bulk_Processing
#WRK=/scratch/owl5022/2024-Krebs_Science/0X_Bulk_Processing
METADATA=$WRK/10phase_of_CTCF_Q1_Q4.txt
###

# Dependencies
# - bedtools
# - java

module load bedtools
module load anaconda3
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Inputs and outputs
GENOME=$WRK/../data/hg38_files/hg38.fa
BAMFILE=$WRK/../data/BAM/BNase-seq_50U-10min_merge_hg38.bam
BLACKLIST=$WRK/../data/hg38_files/ENCFF356LFX_hg38_exclude.bed
MOTIF=$WRK/../data/RefPT-JASPAR
OUTDIR=$WRK/../Library/10phase

# Create output directories if they don't exist
[ -d $OUTDIR ] || mkdir $OUTDIR
# cd to data file
cd OUTDIR=$WRK/Library/
# Script shortcuts
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar
COMPOSITE=$WRK/../bin/sum_Col_CDT.pl
SHUFFLE=$WRK/../bin/shuffle_script.py
chisquare=$WRK/../bin/chisquare.py

# Determine Input
DIRQ1=${OUTDIR}/FOXA_LABEL-uHepG2_SORT-ClosestDyad
Ref=FOXA
# Determine Output
DIR=${OUTDIR}/FOXA_phase_aligned
mkdir -p $DIR


for  XX in RR YY WW SS ; do
        cat $DIRQ1/10xCDT/${XX}_${Ref}_nearNuc_original_phase_0_sense.cdt | cut -f  1-2,253-752 > $DIR/${XX}_${Ref}_nearNuc_phase_0_sense.cdt
        cat $DIRQ1/10xCDT/${XX}_${Ref}_nearNuc_original_phase_1_sense.cdt | cut -f  1-2,254-753 | tail -n +2 > $DIR/${XX}_${Ref}_nearNuc_phase_1_sense.cdt
        cat $DIRQ1/10xCDT/${XX}_${Ref}_nearNuc_original_phase_2_sense.cdt | cut -f  1-2,255-754 | tail -n +2 >  $DIR/${XX}_${Ref}_nearNuc_phase_2_sense.cdt
        cat $DIRQ1/10xCDT/${XX}_${Ref}_nearNuc_original_phase_3_sense.cdt | cut -f  1-2,256-755 | tail -n +2 >  $DIR/${XX}_${Ref}_nearNuc_phase_3_sense.cdt
        cat $DIRQ1/10xCDT/${XX}_${Ref}_nearNuc_original_phase_9_sense.cdt | cut -f  1-2,252-751 | tail -n +2 >  $DIR/${XX}_${Ref}_nearNuc_phase_9_sense.cdt
        cat $DIR/${XX}_${Ref}_nearNuc_phase_0_sense.cdt $DIR/${XX}_${Ref}_nearNuc_phase_1_sense.cdt $DIR/${XX}_${Ref}_nearNuc_phase_2_sense.cdt $DIR/${XX}_${Ref}_nearNuc_phase_3_sense.cdt $DIR/${XX}_${Ref}_nearNuc_phase_9_sense.cdt > $DIR/${XX}_${Ref}_nearNuc_phase_prefered_sense.cdt
        rm $DIR/${XX}_${Ref}_nearNuc_phase_0_sense.cdt $DIR/${XX}_${Ref}_nearNuc_phase_1_sense.cdt $DIR/${XX}_${Ref}_nearNuc_phase_2_sense.cdt $DIR/${XX}_${Ref}_nearNuc_phase_3_sense.cdt $DIR/${XX}_${Ref}_nearNuc_phase_9_sense.cdt
done


for  XX in RR YY WW SS ; do
    perl $COMPOSITE $DIR/${XX}_${Ref}_nearNuc_phase_prefered_sense.cdt $DIR/${XX}_${Ref}_nearNuc_phase_prefered.out  
done


for file in $DIR/*_${Ref}_*_phase_prefered_sense.cdt ; do
    filename=$(basename "$file" "_sense.cdt")
    cat $DIR/${filename}_sense.cdt | cut -f  1-2  > $DIR/${filename}_sense_ref.cdt
    cat $DIR/${filename}_sense.cdt | cut -f  273-372,133-232 | paste $DIR/${filename}_sense_ref.cdt - > $DIR/${filename}_273-372-133-232_sense.cdt
    rm $DIR/${filename}_sense_ref.cdt
done



mkdir -p $DIR/FOXA_phase_prefered
mv $DIR/*232_sense.cdt  $DIR/FOXA_phase_prefered


for file in $DIR/FOXA_phase_prefered/*_phase_prefered_*_sense.cdt  ; do
      filename=$(basename "$file" ".cdt")
      java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum $file -o $DIR/FOXA_phase_prefered/${filename}_SCORES.out 
      tail -n +2 $DIR/FOXA_phase_prefered/${filename}_SCORES.out  > $DIR/FOXA_phase_prefered/${filename}_score.out
      
      tail -n +2 $file | cut -f 3-202 | \
      awk  -v max_col=200 '{
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
        }'  | paste $DIR/FOXA_phase_prefered/${filename}_score.out  - |  awk '{
    OFS="\t";
    print $1,$2,$3/($2+1),$4/($2+1),$5/($2+1),$6/($2+1),$7/($2+1),$8/($2+1),$9/($2+1),$10/($2+1),$11/($2+1),$12/($2+1)}' > $DIR/FOXA_phase_prefered/${filename}_score_peak.out
    rm $DIR/FOXA_phase_prefered/${filename}_SCORES.out 
    rm $DIR/FOXA_phase_prefered/${filename}_score.out
done

for folder in $DIR/FOXA_phase_prefered ; do
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

mkdir -p $WRK/Library/10phase_10bp
for folder in $DIR/FOXA_phase_prefered ; do
    
    if [ -d "$folder" ]; then
        X=$folder
        output_file="${X}.out"
        python $chisquare "$X" "$output_file" 
        mv "$output_file" $WRK/Library/10phase_10bp/
    fi
done




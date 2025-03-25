#!/bin/bash

WRK=/storage/group/bfp2/default/hxc585_HainingChen/2025_Chen_TF-Nuc/X_Bulk_Processing

set -exo
#module load samtools

# Fill in placeholder constants with your directories

# Fill in placeholder constants with your directories
GENOME="$WRK/../data/hg38_files/hg38.fa"
OURDIR=$WRK/Library/Dinucleotide/
mkdir -p $OURDIR
cd $OURDIR

# Script shortcuts
SCRIPTMANAGER="$WRK/../bin/ScriptManager-v0.15.jar"
COMPOSITE=$WRK/../bin/sum_Col_CDT.pl
chisquare=$WRK/../bin/chisquare.py
MOTIFSCAN=$WRK/../bin/scan_FASTA_for_motif_as_binary_string.py

## run Dinucleotode scan, make composite plot of dinucleotide
for file in $WRK/../data/RefPT-JASPAR/1000bp/*1_1000bp.bed $WRK/../data/RefPT-JASPAR/1000bp/*4_1000bp.bed ; do 
    filename=$(basename "$file" ".bed")
     Ref=`basename $file ".cdt" | cut -d "_" -f 1`
     ID=`basename $file ".cdt" | cut -d "_" -f 2`
     Sort=$(basename "$file" ".cdt" | cut -d "_" -f 4 | sed 's/[^0-9]*//g')
     DIR=$OURDIR/${Ref}_${ID}
     mkdir -p $DIR
     [[ -d $DIR/CDT ]] || mkdir $DIR/CDT
     [[ -d $DIR/Composites ]] || mkdir $DIR/Composites
     java -jar $SCRIPTMANAGER sequence-analysis fasta-extract $GENOME $WRK/data/RefPT-JASPAR/1000bp/${filename}.bed -o $DIR/${filename}_1000bp.fa
     python $MOTIFSCAN -i  $DIR/${filename}_1000bp.fa -m AA -o $DIR/CDT/AA_${Ref}_${ID}_Q${Sort}_1000bp
     python $MOTIFSCAN -i  $DIR/${filename}_1000bp.fa -m AT -o $DIR/CDT/AT_${Ref}_${ID}_Q${Sort}_1000bp
     python $MOTIFSCAN -i  $DIR/${filename}_1000bp.fa -m AC -o $DIR/CDT/AC_${Ref}_${ID}_Q${Sort}_1000bp
     python $MOTIFSCAN -i  $DIR/${filename}_1000bp.fa -m AG -o $DIR/CDT/AG_${Ref}_${ID}_Q${Sort}_1000bp
     python $MOTIFSCAN -i  $DIR/${filename}_1000bp.fa -m TA -o $DIR/CDT/TA_${Ref}_${ID}_Q${Sort}_1000bp
     python $MOTIFSCAN -i  $DIR/${filename}_1000bp.fa -m TT -o $DIR/CDT/TT_${Ref}_${ID}_Q${Sort}_1000bp
     python $MOTIFSCAN -i  $DIR/${filename}_1000bp.fa -m TC -o $DIR/CDT/TC_${Ref}_${ID}_Q${Sort}_1000bp
     python $MOTIFSCAN -i  $DIR/${filename}_1000bp.fa -m TG -o $DIR/CDT/TG_${Ref}_${ID}_Q${Sort}_1000bp
     python $MOTIFSCAN -i  $DIR/${filename}_1000bp.fa -m CA -o $DIR/CDT/CA_${Ref}_${ID}_Q${Sort}_1000bp
     python $MOTIFSCAN -i  $DIR/${filename}_1000bp.fa -m CT -o $DIR/CDT/CT_${Ref}_${ID}_Q${Sort}_1000bp
     python $MOTIFSCAN -i  $DIR/${filename}_1000bp.fa -m CC -o $DIR/CDT/CC_${Ref}_${ID}_Q${Sort}_1000bp
     python $MOTIFSCAN -i  $DIR/${filename}_1000bp.fa -m CG -o $DIR/CDT/CG_${Ref}_${ID}_Q${Sort}_1000bp
     python $MOTIFSCAN -i  $DIR/${filename}_1000bp.fa -m GA -o $DIR/CDT/GA_${Ref}_${ID}_Q${Sort}_1000bp
     python $MOTIFSCAN -i  $DIR/${filename}_1000bp.fa -m GT -o $DIR/CDT/GT_${Ref}_${ID}_Q${Sort}_1000bp
     python $MOTIFSCAN -i  $DIR/${filename}_1000bp.fa -m GC -o $DIR/CDT/GC_${Ref}_${ID}_Q${Sort}_1000bp
     python $MOTIFSCAN -i  $DIR/${filename}_1000bp.fa -m GG -o $DIR/CDT/GG_${Ref}_${ID}_Q${Sort}_1000bp
     rm $DIR/CDT/*_${Ref}_${ID}_Q${Sort}_1000bp_anti.cdt
     perl $COMPOSITE $DIR/CDT/AA_${Ref}_${ID}_Q${Sort}_1000bp_sense.cdt $DIR/Composites/AA_${Ref}_${ID}_Q${Sort}_1000bp.out
     perl $COMPOSITE $DIR/CDT/AT_${Ref}_${ID}_Q${Sort}_1000bp_sense.cdt $DIR/Composites/AT_${Ref}_${ID}_Q${Sort}_1000bp.out
     perl $COMPOSITE $DIR/CDT/AC_${Ref}_${ID}_Q${Sort}_1000bp_sense.cdt $DIR/Composites/AC_${Ref}_${ID}_Q${Sort}_1000bp.out
     perl $COMPOSITE $DIR/CDT/AG_${Ref}_${ID}_Q${Sort}_1000bp_sense.cdt $DIR/Composites/AG_${Ref}_${ID}_Q${Sort}_1000bp.out
     perl $COMPOSITE $DIR/CDT/TA_${Ref}_${ID}_Q${Sort}_1000bp_sense.cdt $DIR/Composites/TA_${Ref}_${ID}_Q${Sort}_1000bp.out
     perl $COMPOSITE $DIR/CDT/TT_${Ref}_${ID}_Q${Sort}_1000bp_sense.cdt $DIR/Composites/TT_${Ref}_${ID}_Q${Sort}_1000bp.out
     perl $COMPOSITE $DIR/CDT/TC_${Ref}_${ID}_Q${Sort}_1000bp_sense.cdt $DIR/Composites/TC_${Ref}_${ID}_Q${Sort}_1000bp.out
     perl $COMPOSITE $DIR/CDT/TG_${Ref}_${ID}_Q${Sort}_1000bp_sense.cdt $DIR/Composites/TG_${Ref}_${ID}_Q${Sort}_1000bp.out
     perl $COMPOSITE $DIR/CDT/CA_${Ref}_${ID}_Q${Sort}_1000bp_sense.cdt $DIR/Composites/CA_${Ref}_${ID}_Q${Sort}_1000bp.out
     perl $COMPOSITE $DIR/CDT/CT_${Ref}_${ID}_Q${Sort}_1000bp_sense.cdt $DIR/Composites/CT_${Ref}_${ID}_Q${Sort}_1000bp.out
     perl $COMPOSITE $DIR/CDT/CC_${Ref}_${ID}_Q${Sort}_1000bp_sense.cdt $DIR/Composites/CC_${Ref}_${ID}_Q${Sort}_1000bp.out
     perl $COMPOSITE $DIR/CDT/CG_${Ref}_${ID}_Q${Sort}_1000bp_sense.cdt $DIR/Composites/CG_${Ref}_${ID}_Q${Sort}_1000bp.out
     perl $COMPOSITE $DIR/CDT/GA_${Ref}_${ID}_Q${Sort}_1000bp_sense.cdt $DIR/Composites/GA_${Ref}_${ID}_Q${Sort}_1000bp.out
     perl $COMPOSITE $DIR/CDT/GT_${Ref}_${ID}_Q${Sort}_1000bp_sense.cdt $DIR/Composites/GT_${Ref}_${ID}_Q${Sort}_1000bp.out
     perl $COMPOSITE $DIR/CDT/GC_${Ref}_${ID}_Q${Sort}_1000bp_sense.cdt $DIR/Composites/GC_${Ref}_${ID}_Q${Sort}_1000bp.out
     perl $COMPOSITE $DIR/CDT/GG_${Ref}_${ID}_Q${Sort}_1000bp_sense.cdt $DIR/Composites/GG_${Ref}_${ID}_Q${Sort}_1000bp.out
done

for file in $WRK/data/RefPT-JASPAR/1000bp/*1_1000bp.bed $WRK/data/RefPT-JASPAR/1000bp/*4_1000bp.bed ; do 
    filename=$(basename "$file" ".bed")
     Ref=`basename $file ".cdt" | cut -d "_" -f 1`
     ID=`basename $file ".cdt" | cut -d "_" -f 2`
     Sort=$(basename "$file" ".cdt" | cut -d "_" -f 4 | sed 's/[^0-9]*//g')
     DIR=$OURDIR/${Ref}_${ID}
     for CDT in $DIR/*_${Ref}_${ID}_Q${Sort}_1000bp_sense.cdt ; do
        Dinuc=`basename $CDT ".cdt" | cut -d "_" -f 1`
        cat $DIR/CDT/${Dinuc}_${Ref}_${ID}_Q${Sort}_1000bp_sense.cdt | cut -f  1-2  > $DIR/${Dinuc}_${Ref}_${ID}_Q${Sort}_1000bp_sense_ref.cdt
        cat $DIR/CDT/${Dinuc}_${Ref}_${ID}_Q${Sort}_1000bp_sense.cdt | cut -f  283-482 | paste  $DIR/${Dinuc}_${Ref}_${ID}_Q${Sort}_1000bp_sense_ref.cdt - > $DIR/${Dinuc}_${Ref}_${ID}_Q${Sort}_1000bp_sense_up.cdt
        cat $DIR/CDT/${Dinuc}_${Ref}_${ID}_Q${Sort}_1000bp_sense.cdt | cut -f  523-722 | paste  $DIR/${Dinuc}_${Ref}_${ID}_Q${Sort}_1000bp_sense_ref.cdt - > $DIR/${Dinuc}_${Ref}_${ID}_Q${Sort}_1000bp_sense_down.cdt
        rm $DIR/${Dinuc}_${Ref}_${ID}_Q${Sort}_1000bp_sense_ref.cdt
done
     
done

for file in  $WRK/data/RefPT-JASPAR/1000bp/*1_1000bp.bed $WRK/data/RefPT-JASPAR/1000bp/*4_1000bp.bed  ; do
    filename=$(basename "$file" ".bed")
     Ref=`basename $file ".cdt" | cut -d "_" -f 1`
     ID=`basename $file ".cdt" | cut -d "_" -f 2`
     Sort=$(basename "$file" ".cdt" | cut -d "_" -f 4 | sed 's/[^0-9]*//g')
     DIR=$OURDIR/${Ref}_${ID}
     for CDT in $DIR/*_${Ref}_${ID}_Q${Sort}_1000bp_sense_*.cdt ; do
        Dinuc=`basename $CDT ".cdt" | cut -d "_" -f 1`
        Region=`basename $CDT ".cdt" | cut -d "_" -f 7`
        java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum $DIR/${Dinuc}_${Ref}_${ID}_Q${Sort}_1000bp_sense_${Region}.cdt -o $DIR/${Dinuc}_${Ref}_${ID}_Q${Sort}_1000bp_sense_${Region}_SCORES.out
        tail -n +2 $DIR/${Dinuc}_${Ref}_${ID}_Q${Sort}_1000bp_sense_${Region}_SCORES.out > $DIR/${Dinuc}_${Ref}_${ID}_Q${Sort}_1000bp_sense_${Region}_score.out
        rm $DIR/${Dinuc}_${Ref}_${ID}_Q${Sort}_1000bp_sense_${Region}_SCORES.out
        tail -n +2 $DIR/${Dinuc}_${Ref}_${ID}_Q${Sort}_1000bp_sense_${Region}.cdt | cut -f 3-202 | \
        awk -v max_col=200 '{
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
        }' | paste $DIR/${Dinuc}_${Ref}_${ID}_Q${Sort}_1000bp_sense_${Region}_score.out   - | \
        awk '{
            OFS="\t";
            print $1, $2, $3/($2+1), $4/($2+1), $5/($2+1), $6/($2+1), $7/($2+1), $8/($2+1), $9/($2+1), $10/($2+1), $11/($2+1), $12/($2+1);
        }' > $DIR/${Dinuc}_${Ref}_${ID}_Q${Sort}_1000bp_sense_${Region}_score_peak.out
        rm $DIR/${Dinuc}_${Ref}_${ID}_Q${Sort}_1000bp_sense_${Region}_score.out
    done
done
 

for folder in $OURDIR/*_M  ; do
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

mkdir -p $WRK/Library/Dinucleotide_10bp

for folder in $OURDIR/*_M ; do
    
    if [ -d "$folder" ]; then
        X=$folder
        output_file="${X}.out"
        python $chisquare "$X" "$output_file"
        mv "$output_file" $WRK/Library/Dinucleotide_10bp/
    fi
done






#!/bin/bash

WRK=/storage/group/bfp2/default/hxc585_HainingChen/2025_Chen_TF-Nuc/X_Bulk_Processing
set -exo
#module load samtools

# Fill in placeholder constants with your directories
GENOME="$WRK/../data/hg38_files/hg38.fa"
OURDIR=$WRK/Library/DNAshape/
mkdir -p $OURDIR
cd $OURDIR

# Script shortcuts
SCRIPTMANAGER="$WRK/../bin/ScriptManager-v0.15.jar"
COMPOSITE=$WRK/../bin/sum_Col_CDT.pl
chisquare=$WRK/../bin/chisquare.py
mid_flank_ttest=$WRK/../bin/mid_flank_t-test.py

# run DNA shape analysis, make composite plot of DNA shape
for file in $WRK/../data/RefPT-JASPAR/1000bp/*1_1000bp.bed $WRK/../data/RefPT-JASPAR/1000bp/*4_1000bp.bed ; do 
    filename=$(basename "$file" ".bed")
     Ref=`basename $file ".cdt" | cut -d "_" -f 1`
     ID=`basename $file ".cdt" | cut -d "_" -f 2`
     Sort=$(basename "$file" ".cdt" | cut -d "_" -f 4 | sed 's/[^0-9]*//g')
     DIR=$OURDIR/${Ref}_${ID}
     mkdir -p $DIR
     [[ -d $DIR/CDT ]] || mkdir $DIR/CDT
     [[ -d $DIR/Composites ]] || mkdir $DIR/Composites
     java -jar $SCRIPTMANAGER sequence-analysis dna-shape-bed -p -g -r -l -o $DIR/CDT/${Ref}_${ID}_Q${Sort}_1000bp $GENOME  $WRK/data/RefPT-JASPAR/1000bp/${filename}.bed
     perl $COMPOSITE $DIR/CDT/${Ref}_${ID}_Q${Sort}_1000bp_HelT.cdt $DIR/Composites/${Ref}_${ID}_Q${Sort}_1000bp_HelT.out
     perl $COMPOSITE $DIR/CDT/${Ref}_${ID}_Q${Sort}_1000bp_Roll.cdt $DIR/Composites/${Ref}_${ID}_Q${Sort}_1000bp_Roll.out
     perl $COMPOSITE $DIR/CDT/${Ref}_${ID}_Q${Sort}_1000bp_PropT.cdt $DIR/Composites/${Ref}_${ID}_Q${Sort}_1000bp_PropT.out
     perl $COMPOSITE $DIR/CDT/${Ref}_${ID}_Q${Sort}_1000bp_MGW.cdt $DIR/Composites/${Ref}_${ID}_Q${Sort}_1000bp_MGW.out
done


for file in $WRK/data/RefPT-JASPAR/1000bp/*1_1000bp.bed $WRK/data/RefPT-JASPAR/1000bp/*4_1000bp.bed ; do 
    filename=$(basename "$file" ".bed")
     Ref=`basename $file ".cdt" | cut -d "_" -f 1`
     ID=`basename $file ".cdt" | cut -d "_" -f 2`
     Sort=$(basename "$file" ".cdt" | cut -d "_" -f 4 | sed 's/[^0-9]*//g')
     DIR=$OURDIR/${Ref}_${ID}

     cat $DIR/CDT/${Ref}_${ID}_Q${Sort}_1000bp_HelT.cdt | cut -f  1-2  > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_ref.cdt
     cat $DIR/CDT/${Ref}_${ID}_Q${Sort}_1000bp_HelT.cdt | cut -f  522-721  | paste $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_ref.cdt - > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_down.cdt
     cat $DIR/CDT/${Ref}_${ID}_Q${Sort}_1000bp_HelT.cdt | cut -f  282-481  | paste $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_ref.cdt - > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_up.cdt
     rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_ref.cdt
     cat $DIR/CDT/${Ref}_${ID}_Q${Sort}_1000bp_Roll.cdt | cut -f  1-2  > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_ref.cdt
     cat $DIR/CDT/${Ref}_${ID}_Q${Sort}_1000bp_Roll.cdt | cut -f  522-721  | paste $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_ref.cdt - > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_down.cdt
     cat $DIR/CDT/${Ref}_${ID}_Q${Sort}_1000bp_Roll.cdt | cut -f  282-481  | paste $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_ref.cdt - > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_up.cdt
     rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_ref.cdt
     cat $DIR/CDT/${Ref}_${ID}_Q${Sort}_1000bp_PropT.cdt | cut -f  1-2  > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_ref.cdt
     cat $DIR/CDT/${Ref}_${ID}_Q${Sort}_1000bp_PropT.cdt | cut -f  522-721  | paste $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_ref.cdt - > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_down.cdt
     cat $DIR/CDT/${Ref}_${ID}_Q${Sort}_1000bp_PropT.cdt | cut -f  282-481  | paste $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_ref.cdt - > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_up.cdt
     rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_ref.cdt
     cat $DIR/CDT/${Ref}_${ID}_Q${Sort}_1000bp_MGW.cdt | cut -f  1-2  > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_ref.cdt
     cat $DIR/CDT/${Ref}_${ID}_Q${Sort}_1000bp_MGW.cdt | cut -f  522-721  | paste $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_ref.cdt - > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_down.cdt
     cat $DIR/CDT/${Ref}_${ID}_Q${Sort}_1000bp_MGW.cdt | cut -f  282-481  | paste $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_ref.cdt - > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_up.cdt
     rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_ref.cdt
     
done

for file in $WRK/../data/RefPT-JASPAR/1000bp/*1_1000bp.bed $WRK/../data/RefPT-JASPAR/1000bp/*4_1000bp.bed ; do
    filename=$(basename "$file" ".bed")
     Ref=`basename $file ".cdt" | cut -d "_" -f 1`
     ID=`basename $file ".cdt" | cut -d "_" -f 2`
     Sort=$(basename "$file" ".cdt" | cut -d "_" -f 4 | sed 's/[^0-9]*//g')
     DIR=$OURDIR/${Ref}_${ID}
    java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_down.cdt -o $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_down_SCORES.out
    tail -n +2 $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_down_SCORES.out > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_down_score.out    
    rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_down_SCORES.out

    tail -n +2 $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_down.cdt | cut -f 3-202 | \
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
        }' | paste $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_down_score.out   - | \
        awk '{
            OFS="\t";
            print $1, $2, $3/($2), $4/($2), $5/($2), $6/($2), $7/($2), $8/($2), $9/($2), $10/($2), $11/($2), $12/($2);
        }' > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_down_score_peak.out
        rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_down_score.out



        java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_up.cdt -o $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_up_SCORES.out
        tail -n +2 $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_up_SCORES.out > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_up_score.out   
        rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_up_SCORES.out
        tail -n +2 $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_up.cdt | cut -f 3-202 | \
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
        }' | paste $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_up_score.out   - | \
        awk '{
            OFS="\t";
            print $1, $2, $3/($2), $4/($2), $5/($2), $6/($2), $7/($2), $8/($2), $9/($2), $10/($2), $11/($2), $12/($2);
        }' > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_up_score_peak.out
        rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_up_score.out
done 

for file in $WRK/data/RefPT-JASPAR/1000bp/*1_1000bp.bed $WRK/data/RefPT-JASPAR/1000bp/*4_1000bp.bed ; do
    filename=$(basename "$file" ".bed")
     Ref=`basename $file ".cdt" | cut -d "_" -f 1`
     ID=`basename $file ".cdt" | cut -d "_" -f 2`
     Sort=$(basename "$file" ".cdt" | cut -d "_" -f 4 | sed 's/[^0-9]*//g')
     DIR=$OURDIR/${Ref}_${ID}
    java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_down.cdt -o $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_down_SCORES.out
    tail -n +2 $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_down_SCORES.out > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_down_score.out    
    rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_down_SCORES.out

    tail -n +2 $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_down.cdt | cut -f 3-202  | \
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
        }' | paste $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_down_score.out   - | \
        awk '{
            OFS="\t";
            print $1, $2, $3/($2), $4/($2), $5/($2), $6/($2), $7/($2), $8/($2), $9/($2), $10/($2), $11/($2), $12/($2);
        }' > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_down_score_peak.out
        rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_down_score.out



        java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_up.cdt -o $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_up_SCORES.out
        tail -n +2 $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_up_SCORES.out > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_up_score.out   
        rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_up_SCORES.out
        tail -n +2 $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_up.cdt | cut -f 3-202 | \
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
        }' | paste $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_up_score.out   - | \
        awk '{
            OFS="\t";
            print $1, $2, $3/($2), $4/($2), $5/($2), $6/($2), $7/($2), $8/($2), $9/($2), $10/($2), $11/($2), $12/($2);
        }' > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_up_score_peak.out
        rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_up_score.out
done 

for file in $WRK/data/RefPT-JASPAR/1000bp/*1_1000bp.bed $WRK/data/RefPT-JASPAR/1000bp/*4_1000bp.bed ; do
    filename=$(basename "$file" ".bed")
     Ref=`basename $file ".cdt" | cut -d "_" -f 1`
     ID=`basename $file ".cdt" | cut -d "_" -f 2`
     Sort=$(basename "$file" ".cdt" | cut -d "_" -f 4 | sed 's/[^0-9]*//g')
     DIR=$OURDIR/${Ref}_${ID}
    java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_down.cdt -o $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_down_SCORES.out
    tail -n +2 $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_down_SCORES.out > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_down_score.out    
    rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_down_SCORES.out

    tail -n +2 $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_down.cdt | cut -f 3-202 | \
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
        }' | paste $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_down_score.out   - | \
        awk '{
            OFS="\t";
            print $1, $2, $3/($2), $4/($2), $5/($2), $6/($2), $7/($2), $8/($2), $9/($2), $10/($2), $11/($2), $12/($2);
        }' > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_down_score_peak.out
        rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_down_score.out



        java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_up.cdt -o $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_up_SCORES.out
        tail -n +2 $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_up_SCORES.out > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_up_score.out   
        rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_up_SCORES.out
        tail -n +2 $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_up.cdt | cut -f 3-202 | \
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
        }' | paste $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_up_score.out   - | \
        awk '{
            OFS="\t";
            print $1, $2, $3/($2), $4/($2), $5/($2), $6/($2), $7/($2), $8/($2), $9/($2), $10/($2), $11/($2), $12/($2);
        }' > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_up_score_peak.out
        rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_up_score.out
done

for file in $WRK/data/RefPT-JASPAR/1000bp/*1_1000bp.bed $WRK/data/RefPT-JASPAR/1000bp/*4_1000bp.bed ; do
    filename=$(basename "$file" ".bed")
     Ref=`basename $file ".cdt" | cut -d "_" -f 1`
     ID=`basename $file ".cdt" | cut -d "_" -f 2`
     Sort=$(basename "$file" ".cdt" | cut -d "_" -f 4 | sed 's/[^0-9]*//g')
     DIR=$OURDIR/${Ref}_${ID}
    java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_down.cdt -o $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_down_SCORES.out
    tail -n +2 $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_down_SCORES.out > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_down_score.out    
    rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_down_SCORES.out

    tail -n +2 $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_down.cdt | cut -f 3-202 | \
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
        }' | paste $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_down_score.out   - | \
        awk '{
            OFS="\t";
            print $1, $2, $3/($2), $4/($2), $5/($2), $6/($2), $7/($2), $8/($2), $9/($2), $10/($2), $11/($2), $12/($2);
        }' > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_down_score_peak.out
        rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_down_score.out



        java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_up.cdt -o $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_up_SCORES.out
        tail -n +2 $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_up_SCORES.out > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_up_score.out   
        rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_up_SCORES.out
        tail -n +2 $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_up.cdt | cut -f 3-202 | \
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
        }' | paste $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_up_score.out   - | \
        awk '{
            OFS="\t";
            print $1, $2, $3/($2), $4/($2), $5/($2), $6/($2), $7/($2), $8/($2), $9/($2), $10/($2), $11/($2), $12/($2);
        }' > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_up_score_peak.out
        rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_up_score.out
done

for folder in $OURDIR/*_M* ; do
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
mkdir -p $OURDIR/DNAshape_10bp
for folder in $OURDIR/*_M* ; do
    
    if [ -d "$folder" ]; then
        X=$folder
        output_file="${X}.out"
        python $chisquare "$X" "$output_file" 
        mv $output_file $OURDIR/DNAshape_10bp
    fi
done




for file in $WRK/data/RefPT-JASPAR/1000bp/*1_1000bp.bed $WRK/data/RefPT-JASPAR/1000bp/*4_1000bp.bed ; do 
    filename=$(basename "$file" ".bed")
     Ref=`basename $file ".cdt" | cut -d "_" -f 1`
     ID=`basename $file ".cdt" | cut -d "_" -f 2`
     Sort=$(basename "$file" ".cdt" | cut -d "_" -f 4 | sed 's/[^0-9]*//g')
     DIR=$OURDIR/${Ref}_${ID}

     cat $DIR/CDT/${Ref}_${ID}_Q${Sort}_1000bp_HelT.cdt | cut -f  1-2  > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_ref.cdt
     cat $DIR/CDT/${Ref}_${ID}_Q${Sort}_1000bp_HelT.cdt | cut -f  382-481,522-621  | paste $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_ref.cdt - > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_mid.cdt
     cat $DIR/CDT/${Ref}_${ID}_Q${Sort}_1000bp_HelT.cdt | cut -f  252-351,652-751  | paste $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_ref.cdt - > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_flank.cdt
     cat $DIR/CDT/${Ref}_${ID}_Q${Sort}_1000bp_HelT.cdt | cut -f  252-481,522-751  | paste $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_ref.cdt - > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_total.cdt
     rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_ref.cdt
     cat $DIR/CDT/${Ref}_${ID}_Q${Sort}_1000bp_Roll.cdt | cut -f  1-2  > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_ref.cdt
     cat $DIR/CDT/${Ref}_${ID}_Q${Sort}_1000bp_Roll.cdt | cut -f  382-481,522-621  | paste $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_ref.cdt - > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_mid.cdt
     cat $DIR/CDT/${Ref}_${ID}_Q${Sort}_1000bp_Roll.cdt | cut -f  252-351,652-751   | paste $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_ref.cdt - > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_flank.cdt
     cat $DIR/CDT/${Ref}_${ID}_Q${Sort}_1000bp_Roll.cdt | cut -f  252-481,522-751  | paste $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_ref.cdt - > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_total.cdt
     rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_ref.cdt
     cat $DIR/CDT/${Ref}_${ID}_Q${Sort}_1000bp_PropT.cdt | cut -f  1-2  > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_ref.cdt
     cat $DIR/CDT/${Ref}_${ID}_Q${Sort}_1000bp_PropT.cdt | cut -f  382-481,522-621   | paste $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_ref.cdt - > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_mid.cdt
     cat $DIR/CDT/${Ref}_${ID}_Q${Sort}_1000bp_PropT.cdt | cut -f  252-351,652-751 | paste $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_ref.cdt - > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_flank.cdt
     cat $DIR/CDT/${Ref}_${ID}_Q${Sort}_1000bp_PropT.cdt | cut -f  252-481,522-751  | paste $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_ref.cdt - > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_total.cdt
     rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_ref.cdt
     cat $DIR/CDT/${Ref}_${ID}_Q${Sort}_1000bp_MGW.cdt | cut -f  1-2  > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_ref.cdt
     cat $DIR/CDT/${Ref}_${ID}_Q${Sort}_1000bp_MGW.cdt | cut -f  382-481,522-621  | paste $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_ref.cdt - > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_mid.cdt
     cat $DIR/CDT/${Ref}_${ID}_Q${Sort}_1000bp_MGW.cdt | cut -f  252-351,652-751  | paste $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_ref.cdt - > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_flank.cdt
     cat $DIR/CDT/${Ref}_${ID}_Q${Sort}_1000bp_MGW.cdt | cut -f  252-481,522-751  | paste $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_ref.cdt - > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_total.cdt
     rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_ref.cdt  
done


for file in $WRK/data/RefPT-JASPAR/1000bp/*1_1000bp.bed $WRK/data/RefPT-JASPAR/1000bp/*4_1000bp.bed ; do
    filename=$(basename "$file" ".bed")
     Ref=`basename $file ".cdt" | cut -d "_" -f 1`
     ID=`basename $file ".cdt" | cut -d "_" -f 2`
     Sort=$(basename "$file" ".cdt" | cut -d "_" -f 4 | sed 's/[^0-9]*//g')
     DIR=$OURDIR/${Ref}_${ID}
    java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_mid.cdt -o $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_mid_SCORES.out
    java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_flank.cdt -o $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_flank_SCORES.out
    java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_total.cdt -o $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_total_SCORES.out
    tail -n +2 $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_total_SCORES.out > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_total_sum.out
    rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_total_SCORES.out
    tail -n +2 $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_flank_SCORES.out > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_flank_sum.out
    rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_flank_SCORES.out
    tail -n +2 $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_mid_SCORES.out > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_mid_sum.out
    rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_mid_SCORES.out

    cut -f 2 $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_total_sum.out | paste  $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_mid_sum.out -  | awk -v mid="mid" -v filename="${Ref}_${ID}_Q${Sort}" -v DIR=$DIR 'BEGIN {OFS=","} {print $1, mid, $2/$3 > (DIR "/" filename "_1000bp_MGW_mid_sum.csv")}'

    cut -f 2 $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_total_sum.out | paste  $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_flank_sum.out -  | awk -v flank="flank" -v filename="${Ref}_${ID}_Q${Sort}" -v DIR=$DIR 'BEGIN {OFS=","} {print $1, flank ,$2/$3 > (DIR "/" filename "_1000bp_MGW_flank_sum.csv") }'  

    echo -e "Region,Location,Sum" > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_sum.csv
    cat $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_mid_sum.csv $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_flank_sum.csv >> $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_sum.csv
    rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_mid_sum.csv $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_flank_sum.csv 
    rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_total_sum.out $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_mid_sum.out $DIR/${Ref}_${ID}_Q${Sort}_1000bp_MGW_flank_sum.out

    java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_mid.cdt -o $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_mid_SCORES.out
    java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_flank.cdt -o $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_flank_SCORES.out
    java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_total.cdt -o $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_total_SCORES.out
    tail -n +2 $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_total_SCORES.out > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_total_sum.out
    rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_total_SCORES.out
    tail -n +2 $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_flank_SCORES.out > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_flank_sum.out
    rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_flank_SCORES.out
    tail -n +2 $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_mid_SCORES.out > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_mid_sum.out
    rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_mid_SCORES.out

    cut -f 2 $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_total_sum.out | paste  $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_mid_sum.out -  | awk -v mid="mid" -v filename="${Ref}_${ID}_Q${Sort}" -v DIR=$DIR 'BEGIN {OFS=","} {print $1, mid, $2/$3 > (DIR "/" filename "_1000bp_PropT_mid_sum.csv")}'

    cut -f 2 $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_total_sum.out | paste  $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_flank_sum.out -  | awk -v flank="flank" -v filename="${Ref}_${ID}_Q${Sort}" -v DIR=$DIR 'BEGIN {OFS=","} {print $1, flank ,$2/$3 > (DIR "/" filename "_1000bp_PropT_flank_sum.csv") }'  

    echo -e "Region,Location,Sum" > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_sum.csv
    cat $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_mid_sum.csv $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_flank_sum.csv >> $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_sum.csv
    rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_mid_sum.csv $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_flank_sum.csv 
    rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_total_sum.out $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_mid_sum.out $DIR/${Ref}_${ID}_Q${Sort}_1000bp_PropT_flank_sum.out

    java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_mid.cdt -o $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_mid_SCORES.out
    java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_flank.cdt -o $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_flank_SCORES.out
    java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_total.cdt -o $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_total_SCORES.out
    tail -n +2 $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_total_SCORES.out > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_total_sum.out
    rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_total_SCORES.out
    tail -n +2 $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_flank_SCORES.out > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_flank_sum.out
    rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_flank_SCORES.out
    tail -n +2 $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_mid_SCORES.out > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_mid_sum.out
    rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_mid_SCORES.out

    cut -f 2 $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_total_sum.out | paste  $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_mid_sum.out -  | awk -v mid="mid" -v filename="${Ref}_${ID}_Q${Sort}" -v DIR=$DIR 'BEGIN {OFS=","} {print $1, mid, $2/$3 > (DIR "/" filename "_1000bp_Roll_mid_sum.csv")}'

    cut -f 2 $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_total_sum.out | paste  $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_flank_sum.out -  | awk -v flank="flank" -v filename="${Ref}_${ID}_Q${Sort}" -v DIR=$DIR 'BEGIN {OFS=","} {print $1, flank ,$2/$3 > (DIR "/" filename "_1000bp_Roll_flank_sum.csv") }'  

    echo -e "Region,Location,Sum" > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_sum.csv
    cat $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_mid_sum.csv $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_flank_sum.csv >> $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_sum.csv
    rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_mid_sum.csv $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_flank_sum.csv 
    rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_total_sum.out $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_mid_sum.out $DIR/${Ref}_${ID}_Q${Sort}_1000bp_Roll_flank_sum.out


    java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_mid.cdt -o $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_mid_SCORES.out
    java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_flank.cdt -o $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_flank_SCORES.out
    java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_total.cdt -o $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_total_SCORES.out
    tail -n +2 $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_total_SCORES.out > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_total_sum.out
    rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_total_SCORES.out
    tail -n +2 $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_flank_SCORES.out > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_flank_sum.out
    rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_flank_SCORES.out
    tail -n +2 $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_mid_SCORES.out > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_mid_sum.out
    rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_mid_SCORES.out

    cut -f 2 $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_total_sum.out | paste  $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_mid_sum.out -  | awk -v mid="mid" -v filename="${Ref}_${ID}_Q${Sort}" -v DIR=$DIR 'BEGIN {OFS=","} {print $1, mid, $2/$3 > (DIR "/" filename "_1000bp_HelT_mid_sum.csv")}'

    cut -f 2 $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_total_sum.out | paste  $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_flank_sum.out -  | awk -v flank="flank" -v filename="${Ref}_${ID}_Q${Sort}" -v DIR=$DIR 'BEGIN {OFS=","} {print $1, flank ,$2/$3 > (DIR "/" filename "_1000bp_HelT_flank_sum.csv") }'  

    echo -e "Region,Location,Sum" > $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_sum.csv
    cat $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_mid_sum.csv $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_flank_sum.csv >> $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_sum.csv
    rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_mid_sum.csv $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_flank_sum.csv 
    rm $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_total_sum.out $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_mid_sum.out $DIR/${Ref}_${ID}_Q${Sort}_1000bp_HelT_flank_sum.out

done 

mkdir -p  $WRK/Library/DNAshape_mid_flank 
for folder in $OURDIR/*_M* ; do
    if [ -d "$folder" ]; then
        X=$folder
        output_file="${X}.out"
        python $mid_flank_ttest "$X" "$output_file" 
        mv "$output_file"  $WRK/Library/DNAshape_mid_flank/
    fi
done






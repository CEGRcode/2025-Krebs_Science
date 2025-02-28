#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 00:10:00
#SBATCH -A open
#SBATCH -o logs/3_10phase_of_CTCF_Q1_Q4.log.out-%a
#SBATCH -e logs/3_10phase_of_CTCF_Q1_Q4.log.err-%a
#SBATCH --array 1-4

# data/RefPT-JASPAR
#   |--<TFNAME>_<JASPARID>_SORT-TFnucRatio.bed
#   |--<TFNAME>_<JASPARID>_SORT-TFnucRatio_GROUP-Quartile1.bed
#   |--<TFNAME>_<JASPARID>_SORT-TFnucRatio_GROUP-Quartile2.bed
#   |--<TFNAME>_<JASPARID>_SORT-TFnucRatio_GROUP-Quartile3.bed
#   |--<TFNAME>_<JASPARID>_SORT-TFnucRatio_GROUP-Quartile4.bed

### CHANGE ME
WRK=/path/to/2024-Krebs_Science/03_Call_JASPAR
WRK=/ocean/projects/see180003p/owlang/2024-Krebs_Science/03_Call_JASPAR
WRK=/scratch/owl5022/2024-Krebs_Science/03_Call_JASPAR
METADATA=3_10phase_of_CTCF_Q1_Q4.txt
###

# Dependencies
# - bedtools
# - java

set -exo
module load bedtools
module load anaconda3
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Inputs and outputs
GENOME=$WRK/../data/hg38_files/hg38.fa
BAMFILE=$WRK/../data/BAM/BNase-seq_50U-10min_merge_hg38.bam
BLACKLIST=$WRK/../data/hg38_files/ENCFF356LFX_hg38_exclude.bed
MOTIF=$WRK/../data/RefPT-JASPAR
OUTDIR=$WRK/../data/RefPT-FOXA10phase

# Create output directories if they don't exist
[ -d logs ] || mkdir logs
[ -d $OUTDIR ] || mkdir $OUTDIR
# cd to data file
cd $WRK/../data/
# Script shortcuts
ORIGINAL_SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15-$SLURM_ARRAY_TASK_ID.jar
cp $ORIGINAL_SCRIPTMANAGER $SCRIPTMANAGER
COMPOSITE=$WRK/../bin/sum_Col_CDT.pl
SHUFFLE=$WRK/../bin/shuffle_script.py
chisquare=$WRK/../bin/chisquare.py

# Determine BED file for the current job array index
BEDFILE=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $1}'`
BED=`basename $BEDFILE ".bed"`
Sort=`basename $BEDFILE".bed" | cut -d "_" -f 2 | cut -d "-" -f 2`
Ref=`basename $BEDFILE".bed" | cut -d "_" -f 1`

# Determine BAM file for the current job array index
BAMFILE=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $3}'`
BAM=`basename $BAMFILE ".bam"`

# window
Window=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $2}'`
# Nuc related 
Nuc=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $4}'`
# extract
Extract=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $5}'`
Extract2=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $6}'`
Extract3=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $7}'`
Extract4=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $8}'`
DIR=$OUTDIR/$BED
## Create output directories with ENCODE file name
[ -d $DIR ] || mkdir $DIR
[[ -d $DIR/SCORES ]] || mkdir $DIR/SCORES
[[ -d $DIR/temp ]] || mkdir $DIR/temp

BASE=${BAM}_${BED}_${Nuc}_1000bp_read1
echo $BASE

## add Q group to Ref
awk -v Ref=${Ref} -v Sort=${Sort} '{ OFS="\t"; print $1, $2, $3, $4, $5, $6, Ref"_"Sort }' $BEDFILE | bedtools sort -i | uniq > $DIR/temp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1 $DIR/temp.bed -o $DIR/${BED}_1bp.bed
rm $DIR/temp.bed
# expand
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $DIR/${BED}_1bp.bed -o $DIR/${BED}_${Nuc}_1000bp.bed

# Pileup (read 1)
java -jar $SCRIPTMANAGER read-analysis tag-pileup $DIR/${BED}_${Nuc}_1000bp.bed $BAMFILE --cpu 4 -5 -1 -M $DIR/SCORES/${BASE}_original

## clean 
rm $DIR/${BED}_${Nuc}_1000bp.bed

## negative control by shuffle
if [[ "$Nuc" == "+1Nuc" ]]; then  
    python $SHUFFLE $DIR/SCORES/${BASE}_original_anti.cdt $DIR/SCORES/${BASE}_shuffle_anti.cdt
    python $SHUFFLE $DIR/SCORES/${BASE}_original_sense.cdt $DIR/SCORES/${BASE}_shuffle_sense.cdt    
elif [[ "$Nuc" == "-1Nuc" ]]; then 
    :
fi

for file in $DIR/SCORES/${BASE}_*_sense.cdt ; do
    filename=`basename $file "_sense.cdt"`
    tail -n +2 $DIR/SCORES/${filename}_sense.cdt | cut -f ${Extract} | \
    awk  -v max_col=${Window} '{
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
        }'  | paste $DIR/${BED}_1bp.bed  - >  $DIR/${filename}_sense_1_peak.bed
    tail -n +2 $DIR/SCORES/${filename}_anti.cdt | cut -f ${Extract} | \
    awk  -v max_col=${Window} '{
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
        }'  |  paste $DIR/${filename}_sense_1_peak.bed  - | awk '{ OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,+$21,$22,$23,$24,$25,$26,$27,$8+$26+$9+$27+$10+$18+$11+$19+$12+$20+$13+$21+$14+$22+$15+$23+$16+$24+$17+$25 }' | \
        awk '{ OFS="\t"; print $1, $2, $3, $4, $5, $6, $7, ($8+$26)/($28+1) , ($9+$27)/($28+1) , ($10+$18)/($28+1), ($11+$19)/($28+1), ($12+$20)/($28+1), ($13+$21)/($28+1), ($14+$22)/($28+1), ($15+$23)/($28+1) ,($16+$24)/($28+1),($17+$25)/($28+1) }' > $DIR/${filename}_senseanti_1_peak.bed

        rm $DIR/${filename}_sense_1_peak.bed

        tail -n +2 $DIR/SCORES/${filename}_sense.cdt  | cut -f ${Extract2} | \
        awk  -v Window=$Window '{max_col = Window / 2;
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
        }'  | paste $DIR/${BED}_1bp.bed  - > $DIR/${filename}_sense_2_peak.bed
        tail -n +2 $DIR/SCORES/${filename}_anti.cdt | cut -f ${Extract2} | \
        awk  -v Window=$Window '{max_col = Window / 2;
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
        }'  | paste  $DIR/${filename}_sense_2_peak.bed  - | awk '{ OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,+$21,$22,$23,$24,$25,$26,$27,$8+$26+$9+$27+$10+$18+$11+$19+$12+$20+$13+$21+$14+$22+$15+$23+$16+$24+$17+$25 }' | \
        awk '{ OFS="\t"; print $1, $2, $3, $4, $5, $6, $7, ($8+$26)/($28+1) , ($9+$27)/($28+1) , ($10+$18)/($28+1), ($11+$19)/($28+1), ($12+$20)/($28+1), ($13+$21)/($28+1), ($14+$22)/($28+1), ($15+$23)/($28+1) ,($16+$24)/($28+1),($17+$25)/($28+1) } '  > $DIR/${filename}_senseanti_2_peak.bed
        rm $DIR/${filename}_sense_2_peak.bed
        
        tail -n +2 $DIR/SCORES/${filename}_sense.cdt  | cut -f ${Extract3} | \
        awk  -v Window=$Window '{max_col = Window / 2;
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
        }'  | paste $DIR/${BED}_1bp.bed  - > $DIR/${filename}_sense_3_peak.bed
        
        tail -n +2 $DIR/SCORES/${filename}_anti.cdt | cut -f ${Extract3} | \
        awk  -v Window=$Window '{max_col = Window / 2;
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
        }'  | paste  $DIR/${filename}_sense_3_peak.bed  - | awk '{ OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,+$21,$22,$23,$24,$25,$26,$27,$8+$26+$9+$27+$10+$18+$11+$19+$12+$20+$13+$21+$14+$22+$15+$23+$16+$24+$17+$25 }' | \
        awk '{ OFS="\t"; print $1, $2, $3, $4, $5, $6, $7, ($8+$26)/($28+1) , ($9+$27)/($28+1) , ($10+$18)/($28+1), ($11+$19)/($28+1), ($12+$20)/($28+1), ($13+$21)/($28+1), ($14+$22)/($28+1), ($15+$23)/($28+1) ,($16+$24)/($28+1),($17+$25)/($28+1) }' > $DIR/${filename}_senseanti_3_peak.bed
        rm $DIR/${filename}_sense_3_peak.bed

        tail -n +2 $DIR/SCORES/${filename}_sense.cdt  | cut -f ${Extract4} | \
        awk  -v Window=$Window '{max_col = Window / 2;
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
        }'  | paste $DIR/${BED}_1bp.bed  - > $DIR/${filename}_sense_4_peak.bed
        tail -n +2 $DIR/SCORES/${filename}_anti.cdt | cut -f ${Extract4} | \
        awk  -v Window=$Window '{max_col = Window / 2;
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
        }'  | paste  $DIR/${filename}_sense_4_peak.bed  - | awk '{ OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,+$21,$22,$23,$24,$25,$26,$27,$8+$26+$9+$27+$10+$18+$11+$19+$12+$20+$13+$21+$14+$22+$15+$23+$16+$24+$17+$25 }' | \
         awk '{ OFS="\t"; print $1, $2, $3, $4, $5, $6, $7, ($8+$26)/($28+1) , ($9+$27)/($28+1) , ($10+$18)/($28+1), ($11+$19)/($28+1), ($12+$20)/($28+1), ($13+$21)/($28+1), ($14+$22)/($28+1), ($15+$23)/($28+1) ,($16+$24)/($28+1),($17+$25)/($28+1) } '  >  $DIR/${filename}_senseanti_4_peak.bed
        rm $DIR/${filename}_sense_4_peak.bed
done


cat $DIR/${BASE}_original_senseanti_1_peak.bed  | awk -v DIR=$DIR -v TFNUC=${Ref}_${Nuc} '{
            max_val = $8;
            max_col = 8;
            max_count = 1;  # Initially assume max_val occurs once (in column 8)
            for (i = 9; i <= 17; i++) {
                if ($i > max_val) {
                    max_val = $i;
                    max_col = i;
                    max_count = 1;  # Reset count to 1 for the new maximum
                } else if ($i == max_val) {
                    max_count++;  # Increment count if we find another occurrence of the same max value
                }
            }
            if (max_count == 1) {
                if (max_col == 8) {
                    printf "%s\n", $0 > (DIR "/" TFNUC "_1_10x_+1.bed");
                } 
            else if (max_col == 9) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_1_10x_+2.bed");
                }
            else if (max_col == 10) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_1_10x_+3.bed");
                }
            else if (max_col == 11) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_1_10x_+4.bed");
                }
            else if (max_col == 12) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_1_10x_+5.bed");
                }
            else if (max_col == 13) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_1_10x_+6.bed");
                }
            else if (max_col == 14) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_1_10x_+7.bed");
                }
            else if (max_col == 15) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_1_10x_+8.bed");
                }
            else if (max_col == 16) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_1_10x_+9.bed");
                }
            else if (max_col == 17) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_1_10x_+0.bed");
                } 
            }
        }' 
   
cat $DIR/${BASE}_original_senseanti_2_peak.bed  |  bedtools intersect -v -a - -b $DIR/${Ref}_${Nuc}_1_10x_*.bed |  awk -v DIR=$DIR -v TFNUC=${Ref}_${Nuc} '{
            max_val = $8;
            max_col = 8;
            max_count = 1;  # Initially assume max_val occurs once (in column 8)
            for (i = 9; i <= 17; i++) {
                if ($i > max_val) {
                    max_val = $i;
                    max_col = i;
                    max_count = 1;  # Reset count to 1 for the new maximum
                } else if ($i == max_val) {
                    max_count++;  # Increment count if we find another occurrence of the same max value
                }
            }
            if (max_count == 1) {
                if (max_col == 8) {
                    printf "%s\n", $0 > (DIR "/" TFNUC "_2_10x_+1.bed");
                } 
            else if (max_col == 9) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_2_10x_+2.bed");
                }
            else if (max_col == 10) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_2_10x_+3.bed");
                }
            else if (max_col == 11) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_2_10x_+4.bed");
                }
            else if (max_col == 12) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_2_10x_+5.bed");
                }
            else if (max_col == 13) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_2_10x_+6.bed");
                }
            else if (max_col == 14) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_2_10x_+7.bed");
                }
            else if (max_col == 15) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_2_10x_+8.bed");
                }
            else if (max_col == 16) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_2_10x_+9.bed");
                }
            else if (max_col == 17) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_2_10x_+0.bed");
                } 
            }
        }' 

cat $DIR/${BASE}_original_senseanti_3_peak.bed  |  bedtools intersect -v -a - -b $DIR/${Ref}_${Nuc}_1_10x_*.bed $DIR/${Ref}_${Nuc}_2_10x_*.bed |  awk -v DIR=$DIR -v TFNUC=${Ref}_${Nuc} '{
            max_val = $8;
            max_col = 8;
            max_count = 1;  # Initially assume max_val occurs once (in column 8)
            for (i = 9; i <= 17; i++) {
                if ($i > max_val) {
                    max_val = $i;
                    max_col = i;
                    max_count = 1;  # Reset count to 1 for the new maximum
                } else if ($i == max_val) {
                    max_count++;  # Increment count if we find another occurrence of the same max value
                }
            }
            if (max_count == 1) {
                if (max_col == 8) {
                    printf "%s\n", $0 > (DIR "/" TFNUC "_3_10x_+1.bed");
                } 
            else if (max_col == 9) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_3_10x_+2.bed");
                }
            else if (max_col == 10) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_3_10x_+3.bed");
                }
            else if (max_col == 11) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_3_10x_+4.bed");
                }
            else if (max_col == 12) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_3_10x_+5.bed");
                }
            else if (max_col == 13) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_3_10x_+6.bed");
                }
            else if (max_col == 14) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_3_10x_+7.bed");
                }
            else if (max_col == 15) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_3_10x_+8.bed");
                }
            else if (max_col == 16) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_3_10x_+9.bed");
                }
            else if (max_col == 17) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_3_10x_+0.bed");
                } 
            }
        }' 
cat $DIR/${BASE}_original_senseanti_4_peak.bed  |  bedtools intersect -v -a - -b $DIR/${Ref}_${Nuc}_1_10x_*.bed $DIR/${Ref}_${Nuc}_2_10x_*.bed $DIR/${Ref}_${Nuc}_3_10x_*.bed |  awk -v DIR=$DIR -v TFNUC=${Ref}_${Nuc} '{
            max_val = $8;
            max_col = 8;
            max_count = 1;  # Initially assume max_val occurs once (in column 8)
            for (i = 9; i <= 17; i++) {
                if ($i > max_val) {
                    max_val = $i;
                    max_col = i;
                    max_count = 1;  # Reset count to 1 for the new maximum
                } else if ($i == max_val) {
                    max_count++;  # Increment count if we find another occurrence of the same max value
                }
            }
            if (max_count == 1) {
                if (max_col == 8) {
                    printf "%s\n", $0 > (DIR "/" TFNUC "_4_10x_+1.bed");
                } 
            else if (max_col == 9) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_4_10x_+2.bed");
                }
            else if (max_col == 10) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_4_10x_+3.bed");
                }
            else if (max_col == 11) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_4_10x_+4.bed");
                }
            else if (max_col == 12) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_4_10x_+5.bed");
                }
            else if (max_col == 13) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_4_10x_+6.bed");
                }
            else if (max_col == 14) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_4_10x_+7.bed");
                }
            else if (max_col == 15) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_4_10x_+8.bed");
                }
            else if (max_col == 16) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_4_10x_+9.bed");
                }
            else if (max_col == 17) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_4_10x_+0.bed");
                } 
            }
        }' 

cat $DIR/${Ref}_${Nuc}_*_10x_+1.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_Phase0",$17}' | sort -k8,8nr > $DIR/${Ref}_${Nuc}_${Sort}_phase_0.bed
cat $DIR/${Ref}_${Nuc}_*_10x_+2.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_Phase1",$8}' | sort -k8,8nr > $DIR/${Ref}_${Nuc}_${Sort}_phase_1.bed
cat $DIR/${Ref}_${Nuc}_*_10x_+3.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_Phase2",$9}' | sort -k8,8nr > $DIR/${Ref}_${Nuc}_${Sort}_phase_2.bed
cat $DIR/${Ref}_${Nuc}_*_10x_+4.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_Phase3",$10}' | sort -k8,8nr > $DIR/${Ref}_${Nuc}_${Sort}_phase_3.bed
cat $DIR/${Ref}_${Nuc}_*_10x_+5.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_Phase4",$11}' | sort -k8,8nr > $DIR/${Ref}_${Nuc}_${Sort}_phase_4.bed
cat $DIR/${Ref}_${Nuc}_*_10x_+6.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_Phase5",$12}' | sort -k8,8nr > $DIR/${Ref}_${Nuc}_${Sort}_phase_5.bed
cat $DIR/${Ref}_${Nuc}_*_10x_+7.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_Phase6",$13}' | sort -k8,8nr > $DIR/${Ref}_${Nuc}_${Sort}_phase_6.bed
cat $DIR/${Ref}_${Nuc}_*_10x_+8.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_Phase7",$14}' | sort -k8,8nr > $DIR/${Ref}_${Nuc}_${Sort}_phase_7.bed
cat $DIR/${Ref}_${Nuc}_*_10x_+9.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_Phase8",$15}' | sort -k8,8nr > $DIR/${Ref}_${Nuc}_${Sort}_phase_8.bed
cat $DIR/${Ref}_${Nuc}_*_10x_+0.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_Phase9",$16}' | sort -k8,8nr > $DIR/${Ref}_${Nuc}_${Sort}_phase_9.bed

mv $DIR/${Ref}_${Nuc}_*_10x_*.bed $DIR/temp
rm $DIR/${BASE}_original_senseanti_*_peak.bed
bedtools intersect -v -a $DIR/${BED}_1bp.bed -b $DIR/${Ref}_${Nuc}_${Sort}_phase_*.bed | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_no10x","0"}' > $DIR/${Ref}_${Nuc}_${Sort}_phase_no10x.bed
cat $DIR/${Ref}_${Nuc}_${Sort}_phase_1.bed $DIR/${Ref}_${Nuc}_${Sort}_phase_2.bed $DIR/${Ref}_${Nuc}_${Sort}_phase_3.bed $DIR/${Ref}_${Nuc}_${Sort}_phase_4.bed $DIR/${Ref}_${Nuc}_${Sort}_phase_5.bed $DIR/${Ref}_${Nuc}_${Sort}_phase_6.bed $DIR/${Ref}_${Nuc}_${Sort}_phase_7.bed $DIR/${Ref}_${Nuc}_${Sort}_phase_8.bed $DIR/${Ref}_${Nuc}_${Sort}_phase_9.bed $DIR/${Ref}_${Nuc}_${Sort}_phase_0.bed $DIR/${Ref}_${Nuc}_${Sort}_phase_no10x.bed | bedtools sort -i | uniq > $DIR/${Ref}_${Nuc}_${Sort}_all.bed
## make heatmap figure of original cdt, sorted by phase and phase score.

head -n 1  $DIR/SCORES/${BASE}_original_anti.cdt > $DIR/${BASE}_original_anti_head.cdt
head -n 1  $DIR/SCORES/${BASE}_original_sense.cdt > $DIR/${BASE}_original_sense_head.cdt
tail -n +2  $DIR/SCORES/${BASE}_original_anti.cdt > $DIR/${BASE}_original_anti.cdt
tail -n +2  $DIR/SCORES/${BASE}_original_sense.cdt > $DIR/${BASE}_original_sense.cdt
awk '{OFS="\t"; print $7,$8}' $DIR/${Ref}_${Nuc}_${Sort}_all.bed | paste $DIR/${BASE}_original_anti.cdt -  | sort -k1004,1004nr  > $DIR/${BASE}_original_anti_scoresort.cdt 
awk '{OFS="\t"; print $7,$8}' $DIR/${Ref}_${Nuc}_${Sort}_all.bed | paste $DIR/${BASE}_original_sense.cdt - | sort -k1004,1004nr  > $DIR/${BASE}_original_sense_scoresort.cdt 
rm $DIR/${BASE}_original_anti.cdt $DIR/${BASE}_original_sense.cdt 

awk -v DIR=$DIR -v BASE=${BASE} '{
    if ($1003 ~ /Phase0/) {
        print $0 > (DIR "/" BASE "_original_phase_0_anti.cdt");
    }
    else if ($1003 ~ /Phase1/) {
        print $0 > (DIR "/" BASE "_original_phase_1_anti.cdt");
    }
    else if ($1003 ~ /Phase2/) {
        print $0 > (DIR "/" BASE "_original_phase_2_anti.cdt");
    }
    else if ($1003 ~ /Phase3/) {
        print $0 > (DIR "/" BASE "_original_phase_3_anti.cdt");
    }
    else if ($1003 ~ /Phase4/) {
        print $0 > (DIR "/" BASE "_original_phase_4_anti.cdt");
    }
    else if ($1003 ~ /Phase5/) {
        print $0 > (DIR "/" BASE "_original_phase_5_anti.cdt");
    }
    else if ($1003 ~ /Phase6/) {
        print $0 > (DIR "/" BASE "_original_phase_6_anti.cdt");
    }
    else if ($1003 ~ /Phase7/) {
        print $0 > (DIR "/" BASE "_original_phase_7_anti.cdt");
    }
    else if ($1003 ~ /Phase8/) {
        print $0 > (DIR "/" BASE "_original_phase_8_anti.cdt");
    }
    else if ($1003 ~ /Phase9/) {
        print $0 > (DIR "/" BASE "_original_phase_9_anti.cdt");
    }
    else {
        print $0 > (DIR "/" BASE "_original_phase_no10x_anti.cdt");
    }
}' $DIR/${BASE}_original_anti_scoresort.cdt 

awk -v DIR=$DIR -v BASE=${BASE} '{
    if ($1003 ~ /Phase0/) {
        print $0 > (DIR "/" BASE "_original_phase_0_sense.cdt");
    }
    else if ($1003 ~ /Phase1/) {
        print $0 > (DIR "/" BASE "_original_phase_1_sense.cdt");
    }
    else if ($1003 ~ /Phase2/) {
        print $0 > (DIR "/" BASE "_original_phase_2_sense.cdt");
    }
    else if ($1003 ~ /Phase3/) {
        print $0 > (DIR "/" BASE "_original_phase_3_sense.cdt");
    }
    else if ($1003 ~ /Phase4/) {
        print $0 > (DIR "/" BASE "_original_phase_4_sense.cdt");
    }
    else if ($1003 ~ /Phase5/) {
        print $0 > (DIR "/" BASE "_original_phase_5_sense.cdt");
    }
    else if ($1003 ~ /Phase6/) {
        print $0 > (DIR "/" BASE "_original_phase_6_sense.cdt");
    }
    else if ($1003 ~ /Phase7/) {
        print $0 > (DIR "/" BASE "_original_phase_7_sense.cdt");
    }
    else if ($1003 ~ /Phase8/) {
        print $0 > (DIR "/" BASE "_original_phase_8_sense.cdt");
    }
    else if ($1003 ~ /Phase9/) {
        print $0 > (DIR "/" BASE "_original_phase_9_sense.cdt");
    }
    else {
        print $0 > (DIR "/" BASE "_original_phase_no10x_sense.cdt");
    }
}' $DIR/${BASE}_original_sense_scoresort.cdt 

rm $DIR/${BASE}_original_sense_scoresort.cdt  $DIR/${BASE}_original_anti_scoresort.cdt    

cut -f 1-253  $DIR/${BASE}_original_sense_head.cdt > $DIR/${BASE}_original_sense_head_1-253.cdt
cut -f 1-253  $DIR/${BASE}_original_anti_head.cdt > $DIR/${BASE}_original_anti_head_1-253.cdt
if  [[ "$Nuc" == "+1Nuc" ]]; then  
    cat  $DIR/${BASE}_original_phase_0_sense.cdt $DIR/${BASE}_original_phase_1_sense.cdt $DIR/${BASE}_original_phase_2_sense.cdt $DIR/${BASE}_original_phase_3_sense.cdt $DIR/${BASE}_original_phase_4_sense.cdt $DIR/${BASE}_original_phase_5_sense.cdt $DIR/${BASE}_original_phase_6_sense.cdt $DIR/${BASE}_original_phase_7_sense.cdt $DIR/${BASE}_original_phase_8_sense.cdt $DIR/${BASE}_original_phase_9_sense.cdt $DIR/${BASE}_original_phase_no10x_sense.cdt | \
    cut -f 1-2,502-751 | cat $DIR/${BASE}_original_sense_head_1-253.cdt - >  $DIR/${BASE}_original_250bp_phase_sort_sense.cdt

    cat $DIR/${BASE}_original_phase_0_anti.cdt $DIR/${BASE}_original_phase_1_anti.cdt $DIR/${BASE}_original_phase_2_anti.cdt $DIR/${BASE}_original_phase_3_anti.cdt $DIR/${BASE}_original_phase_4_anti.cdt $DIR/${BASE}_original_phase_5_anti.cdt $DIR/${BASE}_original_phase_6_anti.cdt $DIR/${BASE}_original_phase_7_anti.cdt $DIR/${BASE}_original_phase_8_anti.cdt $DIR/${BASE}_original_phase_9_anti.cdt $DIR/${BASE}_original_phase_no10x_anti.cdt | \
    cut -f 1-2,502-751 | cat $DIR/${BASE}_original_anti_head_1-253.cdt - >  $DIR/${BASE}_original_250bp_phase_sort_anti.cdt

    rm   $DIR/${BASE}_original_sense_head_1-253.cdt  $DIR/${BASE}_original_anti_head_1-253.cdt
    for file in $DIR/${BASE}_original_phase_*_sense.cdt ; do
        filename=$(basename "$file" "_sense.cdt")  
        cat "$DIR/${BASE}_original_sense_head.cdt" "$file" | cut -f 1-1002 > "$DIR/SCORES/${filename}_sense.cdt"
        rm $file
    done
    for file in $DIR/${BASE}_original_phase_*_anti.cdt ; do
        filename=$(basename "$file" "_anti.cdt")  
        cat "$DIR/${BASE}_original_anti_head.cdt" "$file" | cut -f 1-1002 >  "$DIR/SCORES/${filename}_anti.cdt"
        rm $file
    done
    
    rm $DIR/${BASE}_original_anti_head.cdt $DIR/${BASE}_original_sense_head.cdt 

## make heatmap figure

    java -jar $SCRIPTMANAGER figure-generation heatmap -p .95 --blue $DIR/${BASE}_original_250bp_phase_sort_sense.cdt -o  $DIR/${BASE}_original_250bp_phase_sort_sense.png
    java -jar $SCRIPTMANAGER figure-generation heatmap -p .95 --red $DIR/${BASE}_original_250bp_phase_sort_anti.cdt  -o $DIR/${BASE}_original_250bp_phase_sort_anti.png
    java -jar $SCRIPTMANAGER figure-generation merge-heatmap $DIR/${BASE}_original_250bp_phase_sort_sense.png $DIR/${BASE}_original_250bp_phase_sort_anti.png -o $DIR/${BASE}_original_250bp_phase_sort_merge.png
    java -jar $SCRIPTMANAGER figure-generation label-heatmap $DIR/${BASE}_original_250bp_phase_sort_merge.png -f 20 -l 0 -m +125 -r +250  -o $DIR/${BASE}_${Sort}_original_250bp_phase_sort_merge.svg
    java -jar $SCRIPTMANAGER figure-generation label-heatmap $DIR/${BASE}_original_250bp_phase_sort_sense.png -f 20 -l 0 -m +125 -r +250  -o $DIR/${BASE}_original_250bp_phase_sort_sense.svg
    java -jar $SCRIPTMANAGER figure-generation label-heatmap $DIR/${BASE}_original_250bp_phase_sort_anti.png -f 20 -l 0 -m +125 -r +250  -o $DIR/${BASE}_original_250bp_phase_sort_anti.svg
    rm $DIR/${BASE}_original_250bp_phase_sort_sense.png $DIR/${BASE}_original_250bp_phase_sort_anti.png $DIR/${BASE}_original_250bp_phase_sort_merge.png
    rm $DIR/${BASE}_original_250bp_phase_sort_sense.cdt $DIR/${BASE}_original_250bp_phase_sort_anti.cdt
elif [[ "$Nuc" == "-1Nuc" ]]; then 
   cat  $DIR/${BASE}_original_phase_0_sense.cdt $DIR/${BASE}_original_phase_1_sense.cdt $DIR/${BASE}_original_phase_2_sense.cdt $DIR/${BASE}_original_phase_3_sense.cdt $DIR/${BASE}_original_phase_4_sense.cdt $DIR/${BASE}_original_phase_5_sense.cdt $DIR/${BASE}_original_phase_6_sense.cdt $DIR/${BASE}_original_phase_7_sense.cdt $DIR/${BASE}_original_phase_8_sense.cdt $DIR/${BASE}_original_phase_9_sense.cdt $DIR/${BASE}_original_phase_no10x_sense.cdt | \
   cut -f 1-2,253-502 | cat $DIR/${BASE}_original_sense_head_1-253.cdt - >  $DIR/${BASE}_original_250bp_phase_sort_sense.cdt
   cat $DIR/${BASE}_original_phase_0_anti.cdt $DIR/${BASE}_original_phase_1_anti.cdt $DIR/${BASE}_original_phase_2_anti.cdt $DIR/${BASE}_original_phase_3_anti.cdt $DIR/${BASE}_original_phase_4_anti.cdt $DIR/${BASE}_original_phase_5_anti.cdt $DIR/${BASE}_original_phase_6_anti.cdt $DIR/${BASE}_original_phase_7_anti.cdt $DIR/${BASE}_original_phase_8_anti.cdt $DIR/${BASE}_original_phase_9_anti.cdt $DIR/${BASE}_original_phase_no10x_anti.cdt | \
   cut -f 1-2,253-502 | cat $DIR/${BASE}_original_anti_head_1-253.cdt - >  $DIR/${BASE}_original_250bp_phase_sort_anti.cdt

   rm $DIR/${BASE}_original_anti_head_1-253.cdt $DIR/${BASE}_original_sense_head_1-253.cdt

   for file in $DIR/${BASE}_original_phase_*_sense.cdt ; do
      filename=$(basename "$file" "_sense.cdt")  
      cat "$DIR/${BASE}_original_sense_head.cdt" "$file" | cut -f 1-1002 > "$DIR/SCORES/${filename}_sense.cdt"
      rm $file
    done

    for file in $DIR/${BASE}_original_phase_*_anti.cdt ; do
      filename=$(basename "$file" "_anti.cdt")  
       cat "$DIR/${BASE}_original_anti_head.cdt" "$file" | cut -f 1-1002 > "$DIR/SCORES/${filename}_anti.cdt"
      rm $file
    done
    rm   $DIR/${BASE}_original_sense_head.cdt  $DIR/${BASE}_original_anti_head.cdt

## make heatmap figure 

    java -jar $SCRIPTMANAGER figure-generation heatmap -p .95 --blue $DIR/${BASE}_original_250bp_phase_sort_sense.cdt -o  $DIR/${BASE}_original_250bp_phase_sort_sense.png
    java -jar $SCRIPTMANAGER figure-generation heatmap -p .95 --red $DIR/${BASE}_original_250bp_phase_sort_anti.cdt  -o $DIR/${BASE}_original_250bp_phase_sort_anti.png
    java -jar $SCRIPTMANAGER figure-generation merge-heatmap $DIR/${BASE}_original_250bp_phase_sort_sense.png $DIR/${BASE}_original_250bp_phase_sort_anti.png -o $DIR/${BASE}_original_250bp_phase_sort_merge.png
    java -jar $SCRIPTMANAGER figure-generation label-heatmap $DIR/${BASE}_original_250bp_phase_sort_merge.png -f 20 -l -250 -m -125 -r 0  -o $DIR/${BASE}_${Sort}_original_250bp_phase_sort_merge.svg
    java -jar $SCRIPTMANAGER figure-generation label-heatmap $DIR/${BASE}_original_250bp_phase_sort_anti.png -f 20 -l -250 -m -125 -r 0  -o $DIR/${BASE}_original_250bp_phase_sort_anti.svg
    java -jar $SCRIPTMANAGER figure-generation label-heatmap $DIR/${BASE}_original_250bp_phase_sort_sense.png -f 20 -l -250 -m -125 -r 0  -o $DIR/${BASE}_original_250bp_phase_sort_sense.svg

    rm $DIR/${BASE}_original_250bp_phase_sort_sense.png $DIR/${BASE}_original_250bp_phase_sort_anti.png $DIR/${BASE}_original_250bp_phase_sort_merge.png
    rm $DIR/${BASE}_original_250bp_phase_sort_sense.cdt $DIR/${BASE}_original_250bp_phase_sort_anti.cdt
fi


if [[ "$Nuc" == "+1Nuc" ]]; then  
    cat $DIR/${BASE}_shuffle_senseanti_1_peak.bed  | awk -v DIR=$DIR -v TFNUC=${Ref}_${Nuc} '{
            max_val = $8;
            max_col = 8;
            max_count = 1;  # Initially assume max_val occurs once (in column 8)
            for (i = 9; i <= 17; i++) {
                if ($i > max_val) {
                    max_val = $i;
                    max_col = i;
                    max_count = 1;  # Reset count to 1 for the new maximum
                } else if ($i == max_val) {
                    max_count++;  # Increment count if we find another occurrence of the same max value
                }
            }
            if (max_count == 1) {
                if (max_col == 8) {
                    printf "%s\n", $0 > (DIR "/" TFNUC "_1_10x_+1.bed");
                } 
            else if (max_col == 9) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_1_10x_+2.bed");
                }
            else if (max_col == 10) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_1_10x_+3.bed");
                }
            else if (max_col == 11) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_1_10x_+4.bed");
                }
            else if (max_col == 12) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_1_10x_+5.bed");
                }
            else if (max_col == 13) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_1_10x_+6.bed");
                }
            else if (max_col == 14) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_1_10x_+7.bed");
                }
            else if (max_col == 15) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_1_10x_+8.bed");
                }
            else if (max_col == 16) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_1_10x_+9.bed");
                }
            else if (max_col == 17) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_1_10x_+0.bed");
                } 
            }
        }' 
        cat $DIR/${BASE}_shuffle_senseanti_2_peak.bed  |  bedtools intersect -v -a - -b $DIR/${Ref}_${Nuc}_1_10x_*.bed |  awk -v DIR=$DIR -v TFNUC=${Ref}_${Nuc} '{
            max_val = $8;
            max_col = 8;
            max_count = 1;  # Initially assume max_val occurs once (in column 8)
            for (i = 9; i <= 17; i++) {
                if ($i > max_val) {
                    max_val = $i;
                    max_col = i;
                    max_count = 1;  # Reset count to 1 for the new maximum
                } else if ($i == max_val) {
                    max_count++;  # Increment count if we find another occurrence of the same max value
                }
            }
            if (max_count == 1) {
                if (max_col == 8) {
                    printf "%s\n", $0 > (DIR "/" TFNUC "_2_10x_+1.bed");
                } 
            else if (max_col == 9) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_2_10x_+2.bed");
                }
            else if (max_col == 10) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_2_10x_+3.bed");
                }
            else if (max_col == 11) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_2_10x_+4.bed");
                }
            else if (max_col == 12) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_2_10x_+5.bed");
                }
            else if (max_col == 13) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_2_10x_+6.bed");
                }
            else if (max_col == 14) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_2_10x_+7.bed");
                }
            else if (max_col == 15) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_2_10x_+8.bed");
                }
            else if (max_col == 16) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_2_10x_+9.bed");
                }
            else if (max_col == 17) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_2_10x_+0.bed");
                } 
            }
        }' 
        cat $DIR/${BASE}_shuffle_senseanti_3_peak.bed  |  bedtools intersect -v -a - -b $DIR/${Ref}_${Nuc}_1_10x_*.bed $DIR/${Ref}_${Nuc}_2_10x_*.bed |  awk -v DIR=$DIR -v TFNUC=${Ref}_${Nuc} '{
            max_val = $8;
            max_col = 8;
            max_count = 1;  # Initially assume max_val occurs once (in column 8)
            for (i = 9; i <= 17; i++) {
                if ($i > max_val) {
                    max_val = $i;
                    max_col = i;
                    max_count = 1;  # Reset count to 1 for the new maximum
                } else if ($i == max_val) {
                    max_count++;  # Increment count if we find another occurrence of the same max value
                }
            }
            if (max_count == 1) {
                if (max_col == 8) {
                    printf "%s\n", $0 > (DIR "/" TFNUC "_3_10x_+1.bed");
                } 
            else if (max_col == 9) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_3_10x_+2.bed");
                }
            else if (max_col == 10) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_3_10x_+3.bed");
                }
            else if (max_col == 11) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_3_10x_+4.bed");
                }
            else if (max_col == 12) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_3_10x_+5.bed");
                }
            else if (max_col == 13) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_3_10x_+6.bed");
                }
            else if (max_col == 14) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_3_10x_+7.bed");
                }
            else if (max_col == 15) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_3_10x_+8.bed");
                }
            else if (max_col == 16) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_3_10x_+9.bed");
                }
            else if (max_col == 17) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_3_10x_+0.bed");
                } 
            }
        }' 
        cat $DIR/${BASE}_shuffle_senseanti_4_peak.bed  |  bedtools intersect -v -a - -b $DIR/${Ref}_${Nuc}_1_10x_*.bed $DIR/${Ref}_${Nuc}_2_10x_*.bed $DIR/${Ref}_${Nuc}_3_10x_*.bed  |  awk -v DIR=$DIR -v TFNUC=${Ref}_${Nuc} '{
            max_val = $8;
            max_col = 8;
            max_count = 1;  # Initially assume max_val occurs once (in column 8)
            for (i = 9; i <= 17; i++) {
                if ($i > max_val) {
                    max_val = $i;
                    max_col = i;
                    max_count = 1;  # Reset count to 1 for the new maximum
                } else if ($i == max_val) {
                    max_count++;  # Increment count if we find another occurrence of the same max value
                }
            }
            if (max_count == 1) {
                if (max_col == 8) {
                    printf "%s\n", $0 > (DIR "/" TFNUC "_4_10x_+1.bed");
                } 
            else if (max_col == 9) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_4_10x_+2.bed");
                }
            else if (max_col == 10) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_4_10x_+3.bed");
                }
            else if (max_col == 11) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_4_10x_+4.bed");
                }
            else if (max_col == 12) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_4_10x_+5.bed");
                }
            else if (max_col == 13) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_4_10x_+6.bed");
                }
            else if (max_col == 14) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_4_10x_+7.bed");
                }
            else if (max_col == 15) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_4_10x_+8.bed");
                }
            else if (max_col == 16) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_4_10x_+9.bed");
                }
            else if (max_col == 17) {
                printf "%s\n", $0 > (DIR "/" TFNUC "_4_10x_+0.bed");
                } 
            }
        }' 

cat $DIR/${Ref}_${Nuc}_*_10x_+1.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_Phase0",$17}' > $DIR/${Ref}_${Nuc}_shuffle_phase_0.bed
cat $DIR/${Ref}_${Nuc}_*_10x_+2.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_Phase1",$8}'  > $DIR/${Ref}_${Nuc}_shuffle_phase_1.bed
cat $DIR/${Ref}_${Nuc}_*_10x_+3.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_Phase2",$9}'  > $DIR/${Ref}_${Nuc}_shuffle_phase_2.bed
cat $DIR/${Ref}_${Nuc}_*_10x_+4.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_Phase3",$10}'  > $DIR/${Ref}_${Nuc}_shuffle_phase_3.bed
cat $DIR/${Ref}_${Nuc}_*_10x_+5.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_Phase4",$11}' > $DIR/${Ref}_${Nuc}_shuffle_phase_4.bed
cat $DIR/${Ref}_${Nuc}_*_10x_+6.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_Phase5",$12}' > $DIR/${Ref}_${Nuc}_shuffle_phase_5.bed
cat $DIR/${Ref}_${Nuc}_*_10x_+7.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_Phase6",$13}'  > $DIR/${Ref}_${Nuc}_shuffle_phase_6.bed
cat $DIR/${Ref}_${Nuc}_*_10x_+8.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_Phase7",$14}'  > $DIR/${Ref}_${Nuc}_shuffle_phase_7.bed
cat $DIR/${Ref}_${Nuc}_*_10x_+9.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_Phase8",$15}'  > $DIR/${Ref}_${Nuc}_shuffle_phase_8.bed
cat $DIR/${Ref}_${Nuc}_*_10x_+0.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_Phase9",$16}'  > $DIR/${Ref}_${Nuc}_shuffle_phase_9.bed


rm $DIR/${Ref}_${Nuc}_*_10x_*.bed 

bedtools intersect -v -a $DIR/${BED}_1bp.bed -b $DIR/${Ref}_${Nuc}_shuffle_phase_*.bed | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_no10x","0"}' > $DIR/${Ref}_${Nuc}_shuffle_phase_no10x.bed

rm $DIR/${BASE}_shuffle_senseanti_*_peak.bed

cat $DIR/${Ref}_${Nuc}_shuffle_phase_*.bed | bedtools sort -i | uniq > $DIR/${Ref}_${Nuc}_${Sort}_shuffle_all.bed
rm $DIR/${Ref}_${Nuc}_shuffle_phase_*.bed
## make heatmap figure of original cdt, sorted by phase and phase score.
head -n 1  $DIR/SCORES/${BASE}_shuffle_anti.cdt > $DIR/${BASE}_shuffle_anti_head.cdt
head -n 1  $DIR/SCORES/${BASE}_shuffle_sense.cdt > $DIR/${BASE}_shuffle_sense_head.cdt
tail -n +2  $DIR/SCORES/${BASE}_shuffle_anti.cdt > $DIR/${BASE}_shuffle_anti.cdt
tail -n +2  $DIR/SCORES/${BASE}_shuffle_sense.cdt > $DIR/${BASE}_shuffle_sense.cdt
awk '{OFS="\t"; print $7,$8}' $DIR/${Ref}_${Nuc}_${Sort}_shuffle_all.bed | paste $DIR/${BASE}_shuffle_anti.cdt -  | sort -k1004,1004nr  > $DIR/${BASE}_shuffle_anti_scoresort.cdt 
awk '{OFS="\t"; print $7,$8}' $DIR/${Ref}_${Nuc}_${Sort}_shuffle_all.bed | paste $DIR/${BASE}_shuffle_sense.cdt - | sort -k1004,1004nr  > $DIR/${BASE}_shuffle_sense_scoresort.cdt 
rm $DIR/${BASE}_shuffle_anti.cdt $DIR/${BASE}_shuffle_sense.cdt 

awk -v DIR=$DIR -v BASE=${BASE} '{
    if ($1003 ~ /Phase0/) {
        print $0 > (DIR "/" BASE "_shuffle_phase_0_anti.cdt");
    }
    else if ($1003 ~ /Phase1/) {
        print $0 > (DIR "/" BASE "_shuffle_phase_1_anti.cdt");
    }
    else if ($1003 ~ /Phase2/) {
        print $0 > (DIR "/" BASE "_shuffle_phase_2_anti.cdt");
    }
    else if ($1003 ~ /Phase3/) {
        print $0 > (DIR "/" BASE "_shuffle_phase_3_anti.cdt");
    }
    else if ($1003 ~ /Phase4/) {
        print $0 > (DIR "/" BASE "_shuffle_phase_4_anti.cdt");
    }
    else if ($1003 ~ /Phase5/) {
        print $0 > (DIR "/" BASE "_shuffle_phase_5_anti.cdt");
    }
    else if ($1003 ~ /Phase6/) {
        print $0 > (DIR "/" BASE "_shuffle_phase_6_anti.cdt");
    }
    else if ($1003 ~ /Phase7/) {
        print $0 > (DIR "/" BASE "_shuffle_phase_7_anti.cdt");
    }
    else if ($1003 ~ /Phase8/) {
        print $0 > (DIR "/" BASE "_shuffle_phase_8_anti.cdt");
    }
    else if ($1003 ~ /Phase9/) {
        print $0 > (DIR "/" BASE "_shuffle_phase_9_anti.cdt");
    }
    else {
        print $0 > (DIR "/" BASE "_shuffle_phase_no10x_anti.cdt");
    }
}' $DIR/${BASE}_shuffle_anti_scoresort.cdt 


awk -v DIR=$DIR -v BASE=${BASE} '{
    if ($1003 ~ /Phase0/) {
        print $0 > (DIR "/" BASE "_shuffle_phase_0_sense.cdt");
    }
    else if ($1003 ~ /Phase1/) {
        print $0 > (DIR "/" BASE "_shuffle_phase_1_sense.cdt");
    }
    else if ($1003 ~ /Phase2/) {
        print $0 > (DIR "/" BASE "_shuffle_phase_2_sense.cdt");
    }
    else if ($1003 ~ /Phase3/) {
        print $0 > (DIR "/" BASE "_shuffle_phase_3_sense.cdt");
    }
    else if ($1003 ~ /Phase4/) {
        print $0 > (DIR "/" BASE "_shuffle_phase_4_sense.cdt");
    }
    else if ($1003 ~ /Phase5/) {
        print $0 > (DIR "/" BASE "_shuffle_phase_5_sense.cdt");
    }
    else if ($1003 ~ /Phase6/) {
        print $0 > (DIR "/" BASE "_shuffle_phase_6_sense.cdt");
    }
    else if ($1003 ~ /Phase7/) {
        print $0 > (DIR "/" BASE "_shuffle_phase_7_sense.cdt");
    }
    else if ($1003 ~ /Phase8/) {
        print $0 > (DIR "/" BASE "_shuffle_phase_8_sense.cdt");
    }
    else if ($1003 ~ /Phase9/) {
        print $0 > (DIR "/" BASE "_shuffle_phase_9_sense.cdt");
    }
    else {
        print $0 > (DIR "/" BASE "_shuffle_phase_no10x_sense.cdt");
    }
}' $DIR/${BASE}_shuffle_sense_scoresort.cdt 

rm $DIR/${BASE}_shuffle_sense_scoresort.cdt  $DIR/${BASE}_shuffle_anti_scoresort.cdt 

cut -f 1-252  $DIR/${BASE}_shuffle_sense_head.cdt > $DIR/${BASE}_shuffle_sense_head_1-253.cdt
cut -f 1-252  $DIR/${BASE}_shuffle_anti_head.cdt > $DIR/${BASE}_shuffle_anti_head_1-253.cdt


cat  $DIR/${BASE}_shuffle_phase_0_sense.cdt $DIR/${BASE}_shuffle_phase_1_sense.cdt $DIR/${BASE}_shuffle_phase_2_sense.cdt $DIR/${BASE}_shuffle_phase_3_sense.cdt $DIR/${BASE}_shuffle_phase_4_sense.cdt $DIR/${BASE}_shuffle_phase_5_sense.cdt $DIR/${BASE}_shuffle_phase_6_sense.cdt $DIR/${BASE}_shuffle_phase_7_sense.cdt $DIR/${BASE}_shuffle_phase_8_sense.cdt $DIR/${BASE}_shuffle_phase_9_sense.cdt $DIR/${BASE}_shuffle_phase_no10x_sense.cdt | \
cut -f 1-2,502-751 | cat $DIR/${BASE}_shuffle_sense_head_1-253.cdt - >  $DIR/${BASE}_shuffle_250bp_phase_sort_sense.cdt

cat  $DIR/${BASE}_shuffle_phase_0_anti.cdt $DIR/${BASE}_shuffle_phase_1_anti.cdt $DIR/${BASE}_shuffle_phase_2_anti.cdt $DIR/${BASE}_shuffle_phase_3_anti.cdt $DIR/${BASE}_shuffle_phase_4_anti.cdt $DIR/${BASE}_shuffle_phase_5_anti.cdt $DIR/${BASE}_shuffle_phase_6_anti.cdt $DIR/${BASE}_shuffle_phase_7_anti.cdt $DIR/${BASE}_shuffle_phase_8_anti.cdt $DIR/${BASE}_shuffle_phase_9_anti.cdt $DIR/${BASE}_shuffle_phase_no10x_anti.cdt | \
cut -f 1-2,502-751 | cat $DIR/${BASE}_shuffle_anti_head_1-253.cdt - >  $DIR/${BASE}_shuffle_250bp_phase_sort_anti.cdt
rm $DIR/${BASE}_shuffle_sense_head_1-253.cdt  $DIR/${BASE}_shuffle_anti_head_1-253.cdt 

for file in $DIR/${BASE}_shuffle_phase_*_sense.cdt ; do
    filename=$(basename "$file" "_sense.cdt")  
    cat "$DIR/${BASE}_shuffle_sense_head.cdt" "$file" | cut -f 1-1002 > "$DIR/SCORES/${filename}_sense.cdt"
    rm $file
done

for file in $DIR/${BASE}_shuffle_phase_*_anti.cdt ; do
    filename=$(basename "$file" "_anti.cdt")  
    cat "$DIR/${BASE}_shuffle_sense_head.cdt" "$file" | cut -f 1-1002 > "$DIR/SCORES/${filename}_anti.cdt"
    rm $file
done
rm  $DIR/${BASE}_shuffle_sense_head.cdt  $DIR/${BASE}_shuffle_anti_head.cdt

## make heatmap figure of shuffled cdt

java -jar $SCRIPTMANAGER figure-generation heatmap -p .95 --blue $DIR/${BASE}_shuffle_250bp_phase_sort_sense.cdt -o  $DIR/${BASE}_shuffle_250bp_phase_sort_sense.png
java -jar $SCRIPTMANAGER figure-generation heatmap -p .95 --red $DIR/${BASE}_shuffle_250bp_phase_sort_anti.cdt  -o $DIR/${BASE}_shuffle_250bp_phase_sort_anti.png
java -jar $SCRIPTMANAGER figure-generation merge-heatmap $DIR/${BASE}_shuffle_250bp_phase_sort_sense.png $DIR/${BASE}_shuffle_250bp_phase_sort_anti.png -o $DIR/${BASE}_shuffle_250bp_phase_sort_merge.png
java -jar $SCRIPTMANAGER figure-generation label-heatmap $DIR/${BASE}_shuffle_250bp_phase_sort_merge.png -f 20 -l 0 -m +125 -r +250  -o $DIR/${BASE}_${Sort}_shuffle_250bp_phase_sort_merge.svg
rm $DIR/${BASE}_shuffle_250bp_phase_sort_sense.png $DIR/${BASE}_shuffle_250bp_phase_sort_anti.png  $DIR/${BASE}_shuffle_250bp_phase_sort_merge.png
rm $DIR/${BASE}_shuffle_250bp_phase_sort_sense.cdt $DIR/${BASE}_shuffle_250bp_phase_sort_anti.cdt


elif [[ "$Nuc" == "-1Nuc" ]]; then 
    :
fi


mkdir -p $DIR/10plot
mkdir -p $DIR/10Dinucleotide
mkdir -p $DIR/10CDT

for file in $DIR/SCORES/${BASE}_original_phase_*_sense.cdt ; do
    filename=`basename $file "_sense.cdt"`
    perl $COMPOSITE $DIR/SCORES/${filename}_sense.cdt $DIR/${filename}_sense
    perl $COMPOSITE $DIR/SCORES/${filename}_anti.cdt $DIR/${filename}_anti
    tail -1 $DIR/${filename}_sense | cat $DIR/${filename}_anti - > $DIR/10plot/${filename}.out
    rm $DIR/${filename}_sense $DIR/${filename}_anti
done


for file in  $DIR/${Ref}_${Nuc}_${Sort}_phase_*.bed  ; do
  filename=$(basename "$file" ".bed")
  java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $file -o $DIR/10Dinucleotide/${filename}_1000bp.bed
  java -jar $SCRIPTMANAGER sequence-analysis fasta-extract $GENOME   $DIR/10Dinucleotide/${filename}_1000bp.bed -o $DIR/10Dinucleotide/${filename}_1000bp.fa
  rm $DIR/10Dinucleotide/${filename}_1000bp.bed
done


for file in $DIR/${Ref}_${Nuc}_${Sort}_phase_*.bed  ; do
  filename=$(basename "$file" ".bed")
  python $MOTIFSCAN -i  $DIR/10Dinucleotide/${filename}_1000bp.fa -m RR -o $DIR/10CDT/RR_${filename}
  python $MOTIFSCAN -i  $DIR/10Dinucleotide/${filename}_1000bp.fa -m YY -o $DIR/10CDT/YY_${filename}
  python $MOTIFSCAN -i  $DIR/10Dinucleotide/${filename}_1000bp.fa -m WW -o $DIR/10CDT/WW_${filename}
  python $MOTIFSCAN -i  $DIR/10Dinucleotide/${filename}_1000bp.fa -m SS -o $DIR/10CDT/SS_${filename}
done


for file in  $DIR/${Ref}_${Nuc}_${Sort}_phase_*.bed ; do
  filename=$(basename "$file" ".bed")
   perl $COMPOSITE $DIR/10CDT/RR_${filename}_sense.cdt $DIR/RR_${filename}_sense
   perl $COMPOSITE $DIR/10CDT/RR_${filename}_anti.cdt $DIR/RR_${filename}_anti
    tail -1 $DIR/RR_${filename}_sense | cat $DIR/RR_${filename}_anti - > $DIR/10plot/RR_${filename}.out
   rm $DIR/RR_${filename}_sense $DIR/RR_${filename}_anti

   perl $COMPOSITE $DIR/10CDT/YY_${filename}_sense.cdt $DIR/YY_${filename}_sense
   perl $COMPOSITE $DIR/10CDT/YY_${filename}_anti.cdt $DIR/YY_${filename}_anti
    tail -1 $DIR/YY_${filename}_sense | cat $DIR/YY_${filename}_anti - > $DIR/10plot/YY_${filename}.out
   rm $DIR/YY_${filename}_sense $DIR/YY_${filename}_anti

   perl $COMPOSITE $DIR/10CDT/SS_${filename}_sense.cdt $DIR/SS_${filename}_sense
   perl $COMPOSITE $DIR/10CDT/SS_${filename}_anti.cdt $DIR/SS_${filename}_anti
    tail -1 $DIR/SS_${filename}_sense | cat $DIR/SS_${filename}_anti - > $DIR/10plot/SS_${filename}.out
   rm $DIR/SS_${filename}_sense $DIR/SS_${filename}_anti

   perl $COMPOSITE $DIR/10CDT/WW_${filename}_sense.cdt $DIR/WW_${filename}_sense
   perl $COMPOSITE $DIR/10CDT/WW_${filename}_anti.cdt $DIR/WW_${filename}_anti
    tail -1 $DIR/WW_${filename}_sense | cat $DIR/WW_${filename}_anti - > $DIR/10plot/WW_${filename}.out
   rm $DIR/WW_${filename}_sense $DIR/WW_${filename}_anti
done

#Statiscal
for file in  $DIR/${Ref}_${Nuc}_${Sort}_phase_*.bed  ; do
  filename=$(basename "$file" ".bed")
  cat $DIR/10CDT/RR_${filename}_sense.cdt | cut -f  1-2  > $DIR/RR_${filename}_sense_ref.cdt
  cat $DIR/10CDT/RR_${filename}_sense.cdt | cut -f  283-482  | paste $DIR/RR_${filename}_sense_ref.cdt - > $DIR/RR_${filename}_283-482_sense.cdt
  cat $DIR/10CDT/RR_${filename}_sense.cdt | cut -f  523-722  | paste $DIR/RR_${filename}_sense_ref.cdt - > $DIR/RR_${filename}_523-722_sense.cdt
  rm $DIR/RR_${filename}_sense_ref.cdt
  cat $DIR/10CDT/YY_${filename}_sense.cdt | cut -f  1-2  > $DIR/YY_${filename}_sense_ref.cdt
  cat $DIR/10CDT/YY_${filename}_sense.cdt | cut -f  283-482  | paste $DIR/YY_${filename}_sense_ref.cdt - > $DIR/YY_${filename}_283-482_sense.cdt
  cat $DIR/10CDT/YY_${filename}_sense.cdt | cut -f  523-722  | paste $DIR/YY_${filename}_sense_ref.cdt - > $DIR/YY_${filename}_523-722_sense.cdt
  rm $DIR/YY_${filename}_sense_ref.cdt
  cat $DIR/10CDT/SS_${filename}_sense.cdt | cut -f  1-2  > $DIR/SS_${filename}_sense_ref.cdt
  cat $DIR/10CDT/SS_${filename}_sense.cdt | cut -f  283-482  | paste $DIR/SS_${filename}_sense_ref.cdt - > $DIR/SS_${filename}_283-482_sense.cdt
  cat $DIR/10CDT/SS_${filename}_sense.cdt | cut -f  523-722  | paste $DIR/SS_${filename}_sense_ref.cdt - > $DIR/SS_${filename}_523-722_sense.cdt
  rm $DIR/SS_${filename}_sense_ref.cdt
  cat $DIR/10CDT/WW_${filename}_sense.cdt | cut -f  1-2  > $DIR/WW_${filename}_sense_ref.cdt
  cat $DIR/10CDT/WW_${filename}_sense.cdt | cut -f  283-482  | paste $DIR/WW_${filename}_sense_ref.cdt - > $DIR/WW_${filename}_283-482_sense.cdt
  cat $DIR/10CDT/WW_${filename}_sense.cdt | cut -f  523-722  | paste $DIR/WW_${filename}_sense_ref.cdt - > $DIR/WW_${filename}_523-722_sense.cdt
  rm $DIR/WW_${filename}_sense_ref.cdt
done

for file in $DIR/SCORES/${BASE}_original_phase_*_anti.cdt ; do
    filename=$(basename "$file" "_anti.cdt")
    cat $DIR/SCORES/${filename}_sense.cdt | cut -f  1-2  > $DIR/${filename}_sense_ref.cdt
    cat $DIR/SCORES/${filename}_sense.cdt | cut -f  283-482  | paste $DIR/${filename}_sense_ref.cdt - > $DIR/${filename}_283-482_sense.cdt
    cat $DIR/SCORES/${filename}_sense.cdt | cut -f  523-722  | paste $DIR/${filename}_sense_ref.cdt - > $DIR/${filename}_523-722_sense.cdt
    rm $DIR/${filename}_sense_ref.cdt
    cat $DIR/SCORES/${filename}_anti.cdt | cut -f  1-2  > $DIR/${filename}_anti_ref.cdt
    cat $DIR/SCORES/${filename}_anti.cdt | cut -f  283-482  | paste $DIR/${filename}_anti_ref.cdt - > $DIR/${filename}_283-482_anti.cdt
    cat $DIR/SCORES/${filename}_anti.cdt | cut -f  523-722  | paste $DIR/${filename}_anti_ref.cdt - > $DIR/${filename}_523-722_anti.cdt
    rm $DIR/${filename}_anti_ref.cdt
done


mkdir -p $DIR/${Ref}_${Nuc}_${Sort}_BZ
mkdir -p $DIR/${Ref}_${Nuc}_${Sort}_Dinucleotide


mv $DIR/${BAM}_*.cdt $DIR/${Ref}_${Nuc}_${Sort}_BZ
mv $DIR/RR_*.cdt  $DIR/WW_*.cdt $DIR/YY_*.cdt  $DIR/SS_*.cdt $DIR/${Ref}_${Nuc}_${Sort}_Dinucleotide

for file in $DIR/${Ref}_${Nuc}_${Sort}_Dinucleotide/*_sense.cdt ; do
      filename=$(basename "$file" ".cdt")
      java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum $file -o $DIR/${Ref}_${Nuc}_${Sort}_Dinucleotide/${filename}_SCORES.out 
      tail -n +2 $DIR/${Ref}_${Nuc}_${Sort}_Dinucleotide/${filename}_SCORES.out  > $DIR/${Ref}_${Nuc}_${Sort}_Dinucleotide/${filename}_score.out
      
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
        }'  | paste $DIR/${Ref}_${Nuc}_${Sort}_Dinucleotide/${filename}_score.out  - |  awk '{
    OFS="\t";
    print $1,$2,$3/($2+1),$4/($2+1),$5/($2+1),$6/($2+1),$7/($2+1),$8/($2+1),$9/($2+1),$10/($2+1),$11/($2+1),$12/($2+1)}' > $DIR/${Ref}_${Nuc}_${Sort}_Dinucleotide/${filename}_score_peak.out
    rm $DIR/${Ref}_${Nuc}_${Sort}_Dinucleotide/${filename}_score.out
    rm $DIR/${Ref}_${Nuc}_${Sort}_Dinucleotide/${filename}_SCORES.out
  done


for file in $DIR/${Ref}_${Nuc}_${Sort}_BZ/*.cdt  ; do
      filename=$(basename "$file" ".cdt")
      java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum $file -o $DIR/${Ref}_${Nuc}_${Sort}_BZ/${filename}_SCORES.out 
      tail -n +2 $DIR/${Ref}_${Nuc}_${Sort}_BZ/${filename}_SCORES.out  > $DIR/${Ref}_${Nuc}_${Sort}_BZ/${filename}_score.out
      
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
        }'  | paste $DIR/${Ref}_${Nuc}_${Sort}_BZ/${filename}_score.out  - |  awk '{
    OFS="\t";
    print $1,$2,$3/($2+1),$4/($2+1),$5/($2+1),$6/($2+1),$7/($2+1),$8/($2+1),$9/($2+1),$10/($2+1),$11/($2+1),$12/($2+1)}' > $DIR/${Ref}_${Nuc}_${Sort}_BZ/${filename}_score_peak.out
    rm $DIR/${Ref}_${Nuc}_${Sort}_BZ/${filename}_score.out
    rm $DIR/${Ref}_${Nuc}_${Sort}_BZ/${filename}_SCORES.out
  done



  for folder in $DIR/${Ref}_${Nuc}_${Sort}_BZ $DIR/${Ref}_${Nuc}_${Sort}_Dinucleotide ; do
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

for folder in $DIR/${Ref}_${Nuc}_${Sort}_BZ $DIR/${Ref}_${Nuc}_${Sort}_Dinucleotide ; do
    
    if [ -d "$folder" ]; then
        X=$folder
        output_file="${X}.out"
        python $chisquare "$X" "$output_file" 
    fi
done

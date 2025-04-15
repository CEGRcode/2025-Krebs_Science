#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 00:10:00
#SBATCH -A open
#SBATCH -o logs/3_10phase_of_CTCF_Q1_Q4.log.out-%a
#SBATCH -e logs/3_10phase_of_CTCF_Q1_Q4.log.err-%a
#SBATCH --array 1-1

### CHANGE ME
WRK=/storage/group/bfp2/default/hxc585_HainingChen/2025_Chen_TF-Nuc/X_Bulk_Processing
#WRK=/ocean/projects/see180003p/owlang/2024-Krebs_Science/0X_Bulk_Processing
#WRK=/scratch/owl5022/2024-Krebs_Science/0X_Bulk_Processing
METADATA=$WRK/10phase_of_FOXA_uHepG2.txt
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
MOTIF=$WRK/../data/RefPT-Motif/1bp/
OUTDIR=$WRK/Library/10phase

# Create output directories if they don't exist
[ -d logs ] || mkdir logs
[ -d $OUTDIR ] || mkdir $OUTDIR
# cd to data file
cd $WRK/
# Script shortcuts
ORIGINAL_SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15-$SLURM_ARRAY_TASK_ID.jar
cp $ORIGINAL_SCRIPTMANAGER $SCRIPTMANAGER
COMPOSITE=$WRK/../bin/sum_Col_CDT.pl
SHUFFLE=$WRK/../bin/shuffle_script.py
chisquare=$WRK/../bin/chisquare.py
originalshuffleTtest=$WRK/../bin/original-shuffle_t-test.py
MOTIFSCAN=$WRK/../bin/scan_FASTA_for_motif_as_binary_string.py
Nucleosomes=$WRK/../04_Call_Motifs/temp-4_NFIA/Nucleosomes_SORT-Genomic.bed

# Determine BED file for the current job array index
BEDFILE=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $1}'`
BED=`basename $BEDFILE "_1bp.bed"`
Sort=`basename $BEDFILE".bed" | cut -d "_" -f 3 | cut -d "-" -f 2`
Ref=`basename $BEDFILE".bed" | cut -d "_" -f 1`

# Determine BAM file for the current job array index
BAMFILE=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $3}'`
BAM=`basename $BAMFILE ".bam"`

# extract
DIR=$OUTDIR/$BED
## Create output directories with ENCODE file name
[ -d $DIR ] || mkdir $DIR
[[ -d $DIR/SCORES ]] || mkdir $DIR/SCORES
[[ -d $DIR/temp ]] || mkdir $DIR/temp

## add Q group to Ref
awk -v Ref=${Ref} -v Sort=${Sort} '{ OFS="\t"; print $1, $2, $3, $4, $5, $6, Ref"_"Sort }' $BEDFILE | bedtools sort -i | uniq > $DIR/temp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1 $DIR/temp.bed -o $DIR/${BED}_1bp.bed
rm $DIR/temp.bed
## sort by distance to closest nucleosome and pile up with BI
bedtools closest -a $DIR/${BED}_1bp.bed -b $Nucleosomes -d -D a -t first | awk  '{ OFS="\t"; print $1, $2, $3, $4, $14, $6, $7 }'| bedtools sort -i | uniq > $DIR/${BED}_nearNuc.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $DIR/${BED}_nearNuc.bed -o $DIR/${BED}_nearNuc_1000bp.bed
    awk -v DIR=$DIR -v BED=$BED '{
        if (($5 >= 0) && ($5 < 20))      
            print $0 > (DIR"/" BED "_nearNuc_1000bp_523-592.bed");
        else if (($5 >= 20) && ($5 < 40))      
            print $0 > (DIR"/" BED "_nearNuc_1000bp_523-612.bed");
        else if (($5 >= 40) && ($5 < 60))      
            print $0 > (DIR"/" BED "_nearNuc_1000bp_523-632.bed");
        else if (($5 >= 60) && ($5 < 80))      
            print $0 > (DIR"/" BED "_nearNuc_1000bp_523-652.bed");
        else if (($5 >= 80) && ($5 < 100))     
            print $0 > (DIR"/" BED "_nearNuc_1000bp_523-672.bed");
        else if (($5 >= 100) && ($5 < 120))     
            print $0 > (DIR"/" BED "_nearNuc_1000bp_533-692.bed");
        else if (($5 >= 120) && ($5 < 140))     
            print $0 > (DIR"/" BED "_nearNuc_1000bp_553-712.bed");
        else if (($5 <= 0) && ($5 > -20))      
            print $0 > (DIR"/" BED "_nearNuc_1000bp_413-482.bed");
        else if (($5 <= -20) && ($5 > -40))      
            print $0 > (DIR"/" BED "_nearNuc_1000bp_393-482.bed");
        else if (($5 <= -40) && ($5 > -60))     
            print $0 > (DIR"/" BED "_nearNuc_1000bp_373-482.bed");
        else if (($5 <= -60) && ($5 > -80))     
            print $0 > (DIR"/" BED "_nearNuc_1000bp_353-482.bed");
        else if (($5 <= -80) && ($5 > -100))     
            print $0 > (DIR"/" BED "_nearNuc_1000bp_333-482.bed");
        else if (($5 <= -100) && ($5 > -120))     
            print $0 > (DIR"/" BED "_nearNuc_1000bp_313-472.bed");
        else if (($5 <= -120) && ($5 > -140))     
            print $0 > (DIR"/" BED "_nearNuc_1000bp_293-452.bed");
    }' $DIR/${BED}_nearNuc_1000bp.bed

rm $DIR/${BED}_nearNuc_1000bp.bed

## negative control by shuffle
for file in $DIR/${BED}_*_1000bp_*.bed ; do
    filename=`basename $file ".bed"`
    Extract=`basename $filename | cut -d "_" -f 6`
    Nuc=`basename $filename | cut -d "_" -f 4`
    Ref=`basename $filename | cut -d "_" -f 1`
    Sort=`basename $filename | cut -d "_" -f 3 | cut -d "-" -f 2`
    java -jar $SCRIPTMANAGER read-analysis tag-pileup $file $BAMFILE --cpu 4 -5 -1 -M $DIR/SCORES/${BAM}_${filename}_read1_original
    python $SHUFFLE  $DIR/SCORES/${BAM}_${filename}_read1_original_sense.cdt $DIR/SCORES/${BAM}_${filename}_read1_shuffle_sense.cdt
    python $SHUFFLE  $DIR/SCORES/${BAM}_${filename}_read1_original_anti.cdt $DIR/SCORES/${BAM}_${filename}_read1_shuffle_anti.cdt
done
## extract cdt in specific window and determine the phase, and sort by distance to closest nuclosome

for file in $DIR/${BED}_*_1000bp_*.bed ; do
    filename=`basename $file ".bed"`
    Extract=`basename $filename | cut -d "_" -f 6`
    Nuc=`basename $filename | cut -d "_" -f 4`
    Ref=`basename $filename | cut -d "_" -f 1`
    Sort=`basename $filename | cut -d "_" -f 3 | cut -d "-" -f 2`
    tail -n +2 $DIR/SCORES/${BAM}_${filename}_read1_original_sense.cdt | cut -f ${Extract} | \
    awk  -v max_col=160 '{
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
        }'  | paste $file - >  $DIR/${filename}_sense_peak.bed
    tail -n +2 $DIR/SCORES/${BAM}_${filename}_read1_original_anti.cdt | cut -f ${Extract} | \
    awk  -v max_col=160 '{
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
        }'  |  paste $DIR/${filename}_sense_peak.bed - | awk '{ OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,+$21,$22,$23,$24,$25,$26,$27,$8+$26+$9+$27+$10+$18+$11+$19+$12+$20+$13+$21+$14+$22+$15+$23+$16+$24+$17+$25 }' | \
        awk '{ OFS="\t"; print $1, $2, $3, $4, $5, $6, $7, ($8+$26)/($28+1) , ($9+$27)/($28+1) , ($10+$18)/($28+1), ($11+$19)/($28+1), ($12+$20)/($28+1), ($13+$21)/($28+1), ($14+$22)/($28+1), ($15+$23)/($28+1) ,($16+$24)/($28+1),($17+$25)/($28+1) }' > $DIR/${filename}_senseanti_original_peak.bed
        rm $DIR/${filename}_sense_peak.bed


        tail -n +2 $DIR/SCORES/${BAM}_${filename}_read1_shuffle_sense.cdt | cut -f ${Extract} | \
        awk  -v max_col=160 '{
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
        }'  | paste $file - >  $DIR/${filename}_sense_peak.bed
        tail -n +2 $DIR/SCORES/${BAM}_${filename}_read1_shuffle_anti.cdt | cut -f ${Extract} | \
        awk  -v max_col=160 '{
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
        }'  |  paste $DIR/${filename}_sense_peak.bed  - | awk '{ OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,+$21,$22,$23,$24,$25,$26,$27,$8+$26+$9+$27+$10+$18+$11+$19+$12+$20+$13+$21+$14+$22+$15+$23+$16+$24+$17+$25 }' | \
        awk '{ OFS="\t"; print $1, $2, $3, $4, $5, $6, $7, ($8+$26)/($28+1) , ($9+$27)/($28+1) , ($10+$18)/($28+1), ($11+$19)/($28+1), ($12+$20)/($28+1), ($13+$21)/($28+1), ($14+$22)/($28+1), ($15+$23)/($28+1) ,($16+$24)/($28+1),($17+$25)/($28+1) }' > $DIR/${filename}_senseanti_shuffle_peak.bed
        rm  $DIR/${filename}_sense_peak.bed

        cat $DIR/${filename}_senseanti_original_peak.bed  | awk -v DIR=$DIR -v filename=${filename}_original '{
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
                    printf "%s\n", $0 > (DIR "/" filename "_10x_+1.bed");
                } 
            else if (max_col == 9) {
                printf "%s\n", $0 > (DIR "/" filename "_10x_+2.bed");
                }
            else if (max_col == 10) {
                printf "%s\n", $0 > (DIR "/" filename "_10x_+3.bed");
                }
            else if (max_col == 11) {
                printf "%s\n", $0 > (DIR "/" filename "_10x_+4.bed");
                }
            else if (max_col == 12) {
                printf "%s\n", $0 > (DIR "/" filename "_10x_+5.bed");
                }
            else if (max_col == 13) {
                printf "%s\n", $0 > (DIR "/" filename "_10x_+6.bed");
                }
            else if (max_col == 14) {
                printf "%s\n", $0 > (DIR "/" filename "_10x_+7.bed");
                }
            else if (max_col == 15) {
                printf "%s\n", $0 > (DIR "/" filename "_10x_+8.bed");
                }
            else if (max_col == 16) {
                printf "%s\n", $0 > (DIR "/" filename "_10x_+9.bed");
                }
            else if (max_col == 17) {
                printf "%s\n", $0 > (DIR "/" filename "_10x_+0.bed");
                } 
            }
        }' 
        rm $DIR/${filename}_senseanti_original_peak.bed 

        cat $DIR/${filename}_senseanti_shuffle_peak.bed  | awk -v DIR=$DIR -v filename=${filename}_shuffle '{
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
                    printf "%s\n", $0 > (DIR "/" filename "_10x_+1.bed");
                } 
            else if (max_col == 9) {
                printf "%s\n", $0 > (DIR "/" filename "_10x_+2.bed");
                }
            else if (max_col == 10) {
                printf "%s\n", $0 > (DIR "/" filename "_10x_+3.bed");
                }
            else if (max_col == 11) {
                printf "%s\n", $0 > (DIR "/" filename "_10x_+4.bed");
                }
            else if (max_col == 12) {
                printf "%s\n", $0 > (DIR "/" filename "_10x_+5.bed");
                }
            else if (max_col == 13) {
                printf "%s\n", $0 > (DIR "/" filename "_10x_+6.bed");
                }
            else if (max_col == 14) {
                printf "%s\n", $0 > (DIR "/" filename "_10x_+7.bed");
                }
            else if (max_col == 15) {
                printf "%s\n", $0 > (DIR "/" filename "_10x_+8.bed");
                }
            else if (max_col == 16) {
                printf "%s\n", $0 > (DIR "/" filename "_10x_+9.bed");
                }
            else if (max_col == 17) {
                printf "%s\n", $0 > (DIR "/" filename "_10x_+0.bed");
                } 
            }
        }' 
        rm $DIR/${filename}_senseanti_shuffle_peak.bed 
done

cat $DIR/${BED}_nearNuc_1000bp_*_original_10x_+1.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_original_Phase0",$17}'  > $DIR/${Ref}_nearNuc_original_phase_0.bed
cat $DIR/${BED}_nearNuc_1000bp_*_original_10x_+2.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_original_Phase1",$8}'  > $DIR/${Ref}_nearNuc_original_phase_1.bed
cat $DIR/${BED}_nearNuc_1000bp_*_original_10x_+3.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_original_Phase2",$9}' > $DIR/${Ref}_nearNuc_original_phase_2.bed
cat $DIR/${BED}_nearNuc_1000bp_*_original_10x_+4.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_original_Phase3",$10}'  > $DIR/${Ref}_nearNuc_original_phase_3.bed
cat $DIR/${BED}_nearNuc_1000bp_*_original_10x_+5.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_original_Phase4",$11}'  > $DIR/${Ref}_nearNuc_original_phase_4.bed
cat $DIR/${BED}_nearNuc_1000bp_*_original_10x_+6.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_original_Phase5",$12}'  > $DIR/${Ref}_nearNuc_original_phase_5.bed
cat $DIR/${BED}_nearNuc_1000bp_*_original_10x_+7.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_original_Phase6",$13}'  > $DIR/${Ref}_nearNuc_original_phase_6.bed
cat $DIR/${BED}_nearNuc_1000bp_*_original_10x_+8.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_original_Phase7",$14}'  > $DIR/${Ref}_nearNuc_original_phase_7.bed
cat $DIR/${BED}_nearNuc_1000bp_*_original_10x_+9.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_original_Phase8",$15}'  > $DIR/${Ref}_nearNuc_original_phase_8.bed
cat $DIR/${BED}_nearNuc_1000bp_*_original_10x_+0.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_original_Phase9",$16}'  > $DIR/${Ref}_nearNuc_original_phase_9.bed
rm $DIR/${BED}_nearNuc_1000bp_*_original_10x_*.bed 

cat $DIR/${BED}_nearNuc_1000bp_*_shuffle_10x_+1.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_shuffle_Phase0",$17}'  > $DIR/${Ref}_nearNuc_shuffle_phase_0.bed
cat $DIR/${BED}_nearNuc_1000bp_*_shuffle_10x_+2.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_shuffle_Phase1",$8}'  > $DIR/${Ref}_nearNuc_shuffle_phase_1.bed
cat $DIR/${BED}_nearNuc_1000bp_*_shuffle_10x_+3.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_shuffle_Phase2",$9}'  > $DIR/${Ref}_nearNuc_shuffle_phase_2.bed
cat $DIR/${BED}_nearNuc_1000bp_*_shuffle_10x_+4.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_shuffle_Phase3",$10}'  > $DIR/${Ref}_nearNuc_shuffle_phase_3.bed
cat $DIR/${BED}_nearNuc_1000bp_*_shuffle_10x_+5.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_shuffle_Phase4",$11}'  > $DIR/${Ref}_nearNuc_shuffle_phase_4.bed
cat $DIR/${BED}_nearNuc_1000bp_*_shuffle_10x_+6.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_shuffle_Phase5",$12}'  > $DIR/${Ref}_nearNuc_shuffle_phase_5.bed
cat $DIR/${BED}_nearNuc_1000bp_*_shuffle_10x_+7.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_shuffle_Phase6",$13}' > $DIR/${Ref}_nearNuc_shuffle_phase_6.bed
cat $DIR/${BED}_nearNuc_1000bp_*_shuffle_10x_+8.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_shuffle_Phase7",$14}'  > $DIR/${Ref}_nearNuc_shuffle_phase_7.bed
cat $DIR/${BED}_nearNuc_1000bp_*_shuffle_10x_+9.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_shuffle_Phase8",$15}'  > $DIR/${Ref}_nearNuc_shuffle_phase_8.bed
cat $DIR/${BED}_nearNuc_1000bp_*_shuffle_10x_+0.bed | bedtools sort -i | uniq | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_shuffle_Phase9",$16}'  > $DIR/${Ref}_nearNuc_shuffle_phase_9.bed
rm $DIR/${BED}_nearNuc_1000bp_*_shuffle_10x_*.bed 

cat $DIR/${BED}_nearNuc_1000bp_*.bed | bedtools intersect -v -a - -b $DIR/${Ref}_nearNuc_original_phase_*.bed | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_original_no10x","0"}' > $DIR/${Ref}_nearNuc_original_phase_no10x.bed
cat $DIR/${BED}_nearNuc_1000bp_*.bed | bedtools intersect -v -a - -b $DIR/${Ref}_nearNuc_shuffle_phase_*.bed | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7"_shuffle_no10x","0"}' > $DIR/${Ref}_nearNuc_shuffle_phase_no10x.bed
rm $DIR/${BED}_nearNuc_1000bp_*.bed
cat $DIR/${Ref}_nearNuc_original_phase_0.bed  $DIR/${Ref}_nearNuc_original_phase_1.bed  $DIR/${Ref}_nearNuc_original_phase_2.bed  $DIR/${Ref}_nearNuc_original_phase_3.bed  $DIR/${Ref}_nearNuc_original_phase_4.bed  $DIR/${Ref}_nearNuc_original_phase_5.bed $DIR/${Ref}_nearNuc_original_phase_6.bed  $DIR/${Ref}_nearNuc_original_phase_7.bed  $DIR/${Ref}_nearNuc_original_phase_8.bed  $DIR/${Ref}_nearNuc_original_phase_9.bed $DIR/${Ref}_nearNuc_original_phase_no10x.bed | bedtools sort -i | uniq > $DIR/${Ref}_nearNuc_original_all.bed
mv $DIR/${Ref}_nearNuc_original_phase_0.bed  $DIR/${Ref}_nearNuc_original_phase_1.bed  $DIR/${Ref}_nearNuc_original_phase_2.bed  $DIR/${Ref}_nearNuc_original_phase_3.bed  $DIR/${Ref}_nearNuc_original_phase_4.bed  $DIR/${Ref}_nearNuc_original_phase_5.bed $DIR/${Ref}_nearNuc_original_phase_6.bed  $DIR/${Ref}_nearNuc_original_phase_7.bed  $DIR/${Ref}_nearNuc_original_phase_8.bed  $DIR/${Ref}_nearNuc_original_phase_9.bed $DIR/${Ref}_nearNuc_original_phase_no10x.bed $DIR/temp
cat $DIR/${Ref}_nearNuc_shuffle_phase_0.bed  $DIR/${Ref}_nearNuc_shuffle_phase_1.bed  $DIR/${Ref}_nearNuc_shuffle_phase_2.bed  $DIR/${Ref}_nearNuc_shuffle_phase_3.bed  $DIR/${Ref}_nearNuc_shuffle_phase_4.bed  $DIR/${Ref}_nearNuc_shuffle_phase_5.bed $DIR/${Ref}_nearNuc_shuffle_phase_6.bed  $DIR/${Ref}_nearNuc_shuffle_phase_7.bed  $DIR/${Ref}_nearNuc_shuffle_phase_8.bed  $DIR/${Ref}_nearNuc_shuffle_phase_9.bed $DIR/${Ref}_nearNuc_shuffle_phase_no10x.bed | bedtools sort -i | uniq > $DIR/${Ref}_nearNuc_shuffle_all.bed
rm $DIR/${Ref}_nearNuc_shuffle_phase_0.bed  $DIR/${Ref}_nearNuc_shuffle_phase_1.bed  $DIR/${Ref}_nearNuc_shuffle_phase_2.bed  $DIR/${Ref}_nearNuc_shuffle_phase_3.bed  $DIR/${Ref}_nearNuc_shuffle_phase_4.bed  $DIR/${Ref}_nearNuc_shuffle_phase_5.bed $DIR/${Ref}_nearNuc_shuffle_phase_6.bed  $DIR/${Ref}_nearNuc_shuffle_phase_7.bed  $DIR/${Ref}_nearNuc_shuffle_phase_8.bed  $DIR/${Ref}_nearNuc_shuffle_phase_9.bed $DIR/${Ref}_nearNuc_shuffle_phase_no10x.bed 

## make heatmap , attach phase, distance, score to each cdt

 awk -v DIR=$DIR -v BED=$BED '{
        if (($5 >= 0) && ($5 < 20))      
            print $0 > (DIR"/" BED "_nearNuc_1000bp_523-592.bed");
        else if (($5 >= 20) && ($5 < 40))      
            print $0 > (DIR"/" BED "_nearNuc_1000bp_523-612.bed");
        else if (($5 >= 40) && ($5 < 60))      
            print $0 > (DIR"/" BED "_nearNuc_1000bp_523-632.bed");
        else if (($5 >= 60) && ($5 < 80))      
            print $0 > (DIR"/" BED "_nearNuc_1000bp_523-652.bed");
        else if (($5 >= 80) && ($5 < 100))     
            print $0 > (DIR"/" BED "_nearNuc_1000bp_523-672.bed");
        else if (($5 >= 100) && ($5 < 120))     
            print $0 > (DIR"/" BED "_nearNuc_1000bp_533-692.bed");
        else if (($5 >= 120) && ($5 < 140))     
            print $0 > (DIR"/" BED "_nearNuc_1000bp_553-712.bed");
        else if (($5 <= 0) && ($5 > -20))      
            print $0 > (DIR"/" BED "_nearNuc_1000bp_413-482.bed");
        else if (($5 <= -20) && ($5 > -40))      
            print $0 > (DIR"/" BED "_nearNuc_1000bp_393-482.bed");
        else if (($5 <= -40) && ($5 > -60))     
            print $0 > (DIR"/" BED "_nearNuc_1000bp_373-482.bed");
        else if (($5 <= -60) && ($5 > -80))     
            print $0 > (DIR"/" BED "_nearNuc_1000bp_353-482.bed");
        else if (($5 <= -80) && ($5 > -100))     
            print $0 > (DIR"/" BED "_nearNuc_1000bp_333-482.bed");
        else if (($5 <= -100) && ($5 > -120))     
            print $0 > (DIR"/" BED "_nearNuc_1000bp_313-472.bed");
        else if (($5 <= -120) && ($5 > -140))     
            print $0 > (DIR"/" BED "_nearNuc_1000bp_293-452.bed");
    }' $DIR/${Ref}_nearNuc_original_all.bed

    for file in $DIR/${BED}_nearNuc_1000bp_*.bed  ; do
        filename=`basename $file ".bed"`
        tail -n +2  $DIR/SCORES/${BAM}_${filename}_read1_original_anti.cdt > $DIR/${BAM}_${filename}_read1_original_anti.cdt
        tail -n +2  $DIR/SCORES/${BAM}_${filename}_read1_original_sense.cdt > $DIR/${BAM}_${filename}_read1_original_sense.cdt
        awk '{OFS="\t"; print $7,$5,$8}' $file  | paste $DIR/${BAM}_${filename}_read1_original_anti.cdt - | sort -k1004,1004n  >  $DIR/${BAM}_${filename}_read1_original_anti_scoresort.cdt 
        awk '{OFS="\t"; print $7,$5,$8}' $file  | paste $DIR/${BAM}_${filename}_read1_original_sense.cdt - | sort -k1004,1004n  >  $DIR/${BAM}_${filename}_read1_original_sense_scoresort.cdt
        rm $DIR/${BAM}_${filename}_read1_original_anti.cdt $DIR/${BAM}_${filename}_read1_original_sense.cdt
        rm $file
    done

for Nuc in  nearNuc ; do 
    cat $DIR/${BAM}_${BED}_${Nuc}_1000bp_*_read1_original_sense_scoresort.cdt | sort -k1004,1004n | \
    awk -v DIR=$DIR -v BASE=${BAM}_${BED}_${Nuc}_read1 '{
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
    }' 
    rm $DIR/${BAM}_${BED}_${Nuc}_1000bp_*_read1_original_sense_scoresort.cdt

    cat $DIR/${BAM}_${BED}_${Nuc}_1000bp_*_read1_original_anti_scoresort.cdt | sort -k1004,1004n | \
    awk -v DIR=$DIR -v BASE=${BAM}_${BED}_${Nuc}_read1 '{
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
    }'  
    rm $DIR/${BAM}_${BED}_${Nuc}_1000bp_*_read1_original_anti_scoresort.cdt
done

    head -n 1  $DIR/SCORES/${BAM}_${BED}_nearNuc_1000bp_523-592_read1_original_anti.cdt > $DIR/${BAM}_${BED}_read1_original_head.cdt
    cut -f 1-503  $DIR/${BAM}_${BED}_read1_original_head.cdt > $DIR/${BAM}_${BED}_read1_original_head_1-503.cdt
    cat  $DIR/${BAM}_${BED}_nearNuc_read1_original_phase_0_sense.cdt $DIR/${BAM}_${BED}_nearNuc_read1_original_phase_1_sense.cdt $DIR/${BAM}_${BED}_nearNuc_read1_original_phase_2_sense.cdt $DIR/${BAM}_${BED}_nearNuc_read1_original_phase_3_sense.cdt $DIR/${BAM}_${BED}_nearNuc_read1_original_phase_4_sense.cdt $DIR/${BAM}_${BED}_nearNuc_read1_original_phase_5_sense.cdt $DIR/${BAM}_${BED}_nearNuc_read1_original_phase_6_sense.cdt $DIR/${BAM}_${BED}_nearNuc_read1_original_phase_7_sense.cdt $DIR/${BAM}_${BED}_nearNuc_read1_original_phase_8_sense.cdt $DIR/${BAM}_${BED}_nearNuc_read1_original_phase_9_sense.cdt $DIR/${BAM}_${BED}_nearNuc_read1_original_phase_no10x_sense.cdt | \
    cut -f 1-2,253-751 | paste $DIR/${BAM}_${BED}_read1_original_head_1-503.cdt - > $DIR/${BAM}_${BED}_nearNuc_read1_original_500bp_phase_sort_sense.cdt
    cat  $DIR/${BAM}_${BED}_nearNuc_read1_original_phase_0_anti.cdt $DIR/${BAM}_${BED}_nearNuc_read1_original_phase_1_anti.cdt $DIR/${BAM}_${BED}_nearNuc_read1_original_phase_2_anti.cdt $DIR/${BAM}_${BED}_nearNuc_read1_original_phase_3_anti.cdt $DIR/${BAM}_${BED}_nearNuc_read1_original_phase_4_anti.cdt $DIR/${BAM}_${BED}_nearNuc_read1_original_phase_5_anti.cdt $DIR/${BAM}_${BED}_nearNuc_read1_original_phase_6_anti.cdt $DIR/${BAM}_${BED}_nearNuc_read1_original_phase_7_anti.cdt $DIR/${BAM}_${BED}_nearNuc_read1_original_phase_8_anti.cdt $DIR/${BAM}_${BED}_nearNuc_read1_original_phase_9_anti.cdt $DIR/${BAM}_${BED}_nearNuc_read1_original_phase_no10x_anti.cdt | \
    cut -f 1-2,253-751 | paste $DIR/${BAM}_${BED}_read1_original_head_1-503.cdt - > $DIR/${BAM}_${BED}_nearNuc_read1_original_500bp_phase_sort_anti.cdt

    java -jar $SCRIPTMANAGER figure-generation heatmap -p .95 --blue $DIR/${BAM}_${BED}_nearNuc_read1_original_500bp_phase_sort_sense.cdt -o  $DIR/${BAM}_${BED}_nearNuc_read1_original_500bp_phase_sort_sense.png
    java -jar $SCRIPTMANAGER figure-generation heatmap -p .95 --red $DIR/${BAM}_${BED}_nearNuc_read1_original_500bp_phase_sort_anti.cdt  -o $DIR/${BAM}_${BED}_nearNuc_read1_original_500bp_phase_sort_anti.png
    java -jar $SCRIPTMANAGER figure-generation merge-heatmap $DIR/${BAM}_${BED}_nearNuc_read1_original_500bp_phase_sort_sense.png $DIR/${BAM}_${BED}_nearNuc_read1_original_500bp_phase_sort_anti.png -o $DIR/${BAM}_${BED}_nearNuc_read1_original_250bp_phase_sort_merge.png
    java -jar $SCRIPTMANAGER figure-generation label-heatmap $DIR/${BAM}_${BED}_nearNuc_read1_original_250bp_phase_sort_merge.png -f 20 -l -250 -m 0 -r +250  -o $DIR/${BAM}_${BED}_nearNuc_read1_original_250bp_phase_sort_merge.svg
    java -jar $SCRIPTMANAGER figure-generation label-heatmap $DIR/${BAM}_${BED}_nearNuc_read1_original_500bp_phase_sort_sense.png -f 20 -l -250 -m 0 -r +250  -o $DIR/${BAM}_${BED}_nearNuc_read1_original_500bp_phase_sort_sense.svg
    java -jar $SCRIPTMANAGER figure-generation label-heatmap $DIR/${BAM}_${BED}_nearNuc_read1_original_500bp_phase_sort_anti.png -f 20 -l -250 -m 0 -r +250  -o $DIR/${BAM}_${BED}_nearNuc_read1_original_500bp_phase_sort_anti.svg
    rm $DIR/${BAM}_${BED}_nearNuc_read1_original_500bp_phase_sort_sense.png $DIR/${BAM}_${BED}_nearNuc_read1_original_500bp_phase_sort_anti.png $DIR/${BAM}_${BED}_nearNuc_read1_original_250bp_phase_sort_merge.png
    rm $DIR/${BAM}_${BED}_nearNuc_read1_original_500bp_phase_sort_sense.cdt $DIR/${BAM}_${BED}_nearNuc_read1_original_500bp_phase_sort_anti.cdt

    mkdir -p $DIR/10xCDT/
    for file in $DIR/${BAM}_${BED}_nearNuc_read1_original_phase_*_sense.cdt $DIR/${BAM}_${BED}_nearNuc_read1_original_phase_*_anti.cdt ; do
        filename=$(basename "$file" ".cdt")  
        cat "$DIR/${BAM}_${BED}_read1_original_head.cdt" "$file" | cut -f 1-1002 > "$DIR/10xCDT/${filename}.cdt"
        rm $file
    done



## make heatmap for shuffle cdt , attach phase, distance, score to each cdt

awk -v DIR=$DIR -v BED=$BED '{
        if (($5 >= 0) && ($5 < 20))      
            print $0 > (DIR"/" BED "_nearNuc_1000bp_523-592.bed");
        else if (($5 >= 20) && ($5 < 40))      
            print $0 > (DIR"/" BED "_nearNuc_1000bp_523-612.bed");
        else if (($5 >= 40) && ($5 < 60))      
            print $0 > (DIR"/" BED "_nearNuc_1000bp_523-632.bed");
        else if (($5 >= 60) && ($5 < 80))      
            print $0 > (DIR"/" BED "_nearNuc_1000bp_523-652.bed");
        else if (($5 >= 80) && ($5 < 100))     
            print $0 > (DIR"/" BED "_nearNuc_1000bp_523-672.bed");
        else if (($5 >= 100) && ($5 < 120))     
            print $0 > (DIR"/" BED "_nearNuc_1000bp_533-692.bed");
        else if (($5 >= 120) && ($5 < 140))     
            print $0 > (DIR"/" BED "_nearNuc_1000bp_553-712.bed");
        else if (($5 <= 0) && ($5 > -20))      
            print $0 > (DIR"/" BED "_nearNuc_1000bp_413-482.bed");
        else if (($5 <= -20) && ($5 > -40))      
            print $0 > (DIR"/" BED "_nearNuc_1000bp_393-482.bed");
        else if (($5 <= -40) && ($5 > -60))     
            print $0 > (DIR"/" BED "_nearNuc_1000bp_373-482.bed");
        else if (($5 <= -60) && ($5 > -80))     
            print $0 > (DIR"/" BED "_nearNuc_1000bp_353-482.bed");
        else if (($5 <= -80) && ($5 > -100))     
            print $0 > (DIR"/" BED "_nearNuc_1000bp_333-482.bed");
        else if (($5 <= -100) && ($5 > -120))     
            print $0 > (DIR"/" BED "_nearNuc_1000bp_313-472.bed");
        else if (($5 <= -120) && ($5 > -140))     
            print $0 > (DIR"/" BED "_nearNuc_1000bp_293-452.bed");
    }' $DIR/${Ref}_nearNuc_shuffle_all.bed


for file in $DIR/${BED}_nearNuc_1000bp_*.bed  ; do
        filename=`basename $file ".bed"`
        tail -n +2  $DIR/SCORES/${BAM}_${filename}_read1_shuffle_anti.cdt > $DIR/${BAM}_${filename}_read1_shuffle_anti.cdt
        tail -n +2  $DIR/SCORES/${BAM}_${filename}_read1_shuffle_sense.cdt > $DIR/${BAM}_${filename}_read1_shuffle_sense.cdt
        awk '{OFS="\t"; print $7,$5,$8}' $file  | paste $DIR/${BAM}_${filename}_read1_shuffle_anti.cdt - | sort -k1004,1004n  >  $DIR/${BAM}_${filename}_read1_shuffle_anti_scoresort.cdt 
        awk '{OFS="\t"; print $7,$5,$8}' $file  | paste $DIR/${BAM}_${filename}_read1_shuffle_sense.cdt - | sort -k1004,1004n  >  $DIR/${BAM}_${filename}_read1_shuffle_sense_scoresort.cdt
        rm $DIR/${BAM}_${filename}_read1_shuffle_anti.cdt $DIR/${BAM}_${filename}_read1_shuffle_sense.cdt
        rm $file
done

for Nuc in  nearNuc ; do 
    cat $DIR/${BAM}_${BED}_${Nuc}_1000bp_*_read1_shuffle_sense_scoresort.cdt | sort -k1004,1004n | \
    awk -v DIR=$DIR -v BASE=${BAM}_${BED}_${Nuc}_read1 '{
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
    }' 
    rm $DIR/${BAM}_${BED}_${Nuc}_1000bp_*_read1_shuffle_sense_scoresort.cdt

    cat $DIR/${BAM}_${BED}_${Nuc}_1000bp_*_read1_shuffle_anti_scoresort.cdt | sort -k1004,1004n | \
    awk -v DIR=$DIR -v BASE=${BAM}_${BED}_${Nuc}_read1 '{
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
    }'  
    rm $DIR/${BAM}_${BED}_${Nuc}_1000bp_*_read1_shuffle_anti_scoresort.cdt
done

    head -n 1  $DIR/SCORES/${BAM}_${BED}_nearNuc_1000bp_523-592_read1_shuffle_anti.cdt > $DIR/${BAM}_${BED}_read1_shuffle_head.cdt
    cut -f 1-503  $DIR/${BAM}_${BED}_read1_shuffle_head.cdt > $DIR/${BAM}_${BED}_read1_shuffle_head_1-503.cdt
    cat  $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_phase_0_sense.cdt $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_phase_1_sense.cdt $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_phase_2_sense.cdt $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_phase_3_sense.cdt $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_phase_4_sense.cdt $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_phase_5_sense.cdt $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_phase_6_sense.cdt $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_phase_7_sense.cdt $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_phase_8_sense.cdt $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_phase_9_sense.cdt $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_phase_no10x_sense.cdt | \
    cut -f 1-2,253-751 | paste $DIR/${BAM}_${BED}_read1_shuffle_head_1-503.cdt - > $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_500bp_phase_sort_sense.cdt
    cat  $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_phase_0_anti.cdt $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_phase_1_anti.cdt $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_phase_2_anti.cdt $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_phase_3_anti.cdt $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_phase_4_anti.cdt $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_phase_5_anti.cdt $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_phase_6_anti.cdt $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_phase_7_anti.cdt $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_phase_8_anti.cdt $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_phase_9_anti.cdt $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_phase_no10x_anti.cdt | \
    cut -f 1-2,253-751 | paste $DIR/${BAM}_${BED}_read1_shuffle_head_1-503.cdt - > $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_500bp_phase_sort_anti.cdt

    java -jar $SCRIPTMANAGER figure-generation heatmap -p .95 --blue $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_500bp_phase_sort_sense.cdt -o  $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_500bp_phase_sort_sense.png
    java -jar $SCRIPTMANAGER figure-generation heatmap -p .95 --red $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_500bp_phase_sort_anti.cdt  -o $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_500bp_phase_sort_anti.png
    java -jar $SCRIPTMANAGER figure-generation merge-heatmap $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_500bp_phase_sort_sense.png $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_500bp_phase_sort_anti.png -o $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_250bp_phase_sort_merge.png
    java -jar $SCRIPTMANAGER figure-generation label-heatmap $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_250bp_phase_sort_merge.png -f 20 -l -250 -m 0 -r +250  -o $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_250bp_phase_sort_merge.svg
    java -jar $SCRIPTMANAGER figure-generation label-heatmap $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_500bp_phase_sort_sense.png -f 20 -l -250 -m 0 -r +250  -o $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_500bp_phase_sort_sense.svg
    java -jar $SCRIPTMANAGER figure-generation label-heatmap $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_500bp_phase_sort_anti.png -f 20 -l -250 -m 0 -r +250  -o $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_500bp_phase_sort_anti.svg
    rm $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_250bp_phase_sort_merge.png $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_500bp_phase_sort_sense.png $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_500bp_phase_sort_anti.png
    rm $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_500bp_phase_sort_sense.cdt  $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_500bp_phase_sort_anti.cdt

    mkdir -p $DIR/10xCDT/
    for file in $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_phase_*_anti.cdt $DIR/${BAM}_${BED}_nearNuc_read1_shuffle_phase_*_sense.cdt ; do
        filename=$(basename "$file" ".cdt")  
        cat "$DIR/${BAM}_${BED}_read1_shuffle_head.cdt" "$file" | cut -f 1-1002 > "$DIR/10xCDT/${filename}.cdt"
        rm $file
    done

## make heatmap of adjust nucleosome dyad.

java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1 $DIR/${Ref}_nearNuc_original_all.bed -o $DIR/${Ref}_nearNuc_original_all_1bp.bed
    bedtools closest -D a -t first -a $DIR/${Ref}_nearNuc_original_all_1bp.bed -b ../04_Call_Motifs/temp-4_NFIA/Nucleosomes_SORT-Genomic.bed | awk '{ OFS="\t"; print $1, $2, $3, $4, $5, $6, $7, $10, $8 }' | \
    awk -v DIR=$DIR -v BASE=${Ref}_nearNuc_original '{
    if ($5 > 0) {
        print $0 > (DIR "/" BASE "_down.bed");
    }
    else if ($5 < 0) {
        print $0 > (DIR "/" BASE "_up.bed");
    }
    }'  

    awk '{ OFS="\t"; print $1, $2, $3, $4, $5+1, $6, $7, $8 , $9}' $DIR/${Ref}_nearNuc_original_down.bed  > $DIR/${Ref}_nearNuc_original_down_1.bed
    awk '{ OFS="\t"; print $1, $2, $3, $4, $5-1, $6, $7, $8 , $9}' $DIR/${Ref}_nearNuc_original_up.bed  > $DIR/${Ref}_nearNuc_original_up_1.bed

    rm $DIR/${Ref}_nearNuc_original_down.bed $DIR/${Ref}_nearNuc_original_up.bed $DIR/${Ref}_nearNuc_original_all_1bp.bed

    cat $DIR/${Ref}_nearNuc_original_down_1.bed $DIR/${Ref}_nearNuc_original_up_1.bed | bedtools sort -i | uniq | \
    awk -v DIR=$DIR -v BASE=${Ref}_nearNuc_original '{
    if ($7 ~ /Phase0/) {
        print $0 > (DIR "/" BASE "_phase0.bed");
    }
    else if ($7 ~ /Phase1/) {
        print $0 > (DIR "/" BASE "_phase1.bed");
    }
    else if ($7 ~ /Phase2/) {
        print $0 > (DIR "/" BASE "_phase2.bed");
    }
    else if ($7 ~ /Phase3/) {
        print $0 > (DIR "/" BASE "_phase3.bed");
    }
    else if ($7 ~ /Phase4/) {
        print $0 > (DIR "/" BASE "_phase4.bed");
    }
    else if ($7 ~ /Phase5/) {
        print $0 > (DIR "/" BASE "_phase5.bed");
    }
    else if ($7 ~ /Phase6/) {
        print $0 > (DIR "/" BASE "_phase6.bed");
    }
    else if ($7 ~ /Phase7/) {
        print $0 > (DIR "/" BASE "_phase7.bed");
    }
    else if ($7 ~ /Phase8/) {
        print $0 > (DIR "/" BASE "_phase8.bed");
    }
    else if ($7 ~ /Phase9/) {
        print $0 > (DIR "/" BASE "_phase9.bed");
    }
    else {
        print $0 > (DIR "/" BASE "_phaseno10x.bed");
    }
    }'  
    rm $DIR/${Ref}_nearNuc_original_down_1.bed $DIR/${Ref}_nearNuc_original_up_1.bed
for file in $DIR/${Ref}_nearNuc_original_phase*.bed ; do
    filename=$(basename "$file" .bed)
    awk -v DIR=$DIR -v filename="$filename" '{
        abs_value = ($5 < 0) ? -$5 : $5
        last_digit = abs_value % 10
        if ($5 >= 0) {
            print $0 >> (DIR "/" filename "_" last_digit ".bed")
        } else {
            new_digit = (last_digit == 0) ? 0 : 10 - last_digit
            print $0 >> (DIR "/" filename "_" new_digit ".bed")
        }
    }' "$file"
done
    rm $DIR/${Ref}_nearNuc_original_phaseno10x_*.bed

    cat ${DIR}/${Ref}_nearNuc_original_phase1_1.bed ${DIR}/${Ref}_nearNuc_original_phase2_2.bed ${DIR}/${Ref}_nearNuc_original_phase3_3.bed ${DIR}/${Ref}_nearNuc_original_phase4_4.bed ${DIR}/${Ref}_nearNuc_original_phase5_5.bed ${DIR}/${Ref}_nearNuc_original_phase6_6.bed ${DIR}/${Ref}_nearNuc_original_phase7_7.bed ${DIR}/${Ref}_nearNuc_original_phase8_8.bed ${DIR}/${Ref}_nearNuc_original_phase9_9.bed ${DIR}/${Ref}_nearNuc_original_phase0_0.bed ${DIR}/${Ref}_nearNuc_original_phaseno10x.bed | awk '{OFS="\t"} {print $1,$8,$8,$1"_"$8"_"$8,$9,$6,$1,$2,$3,$4,$5,$6,$7}' > $DIR/nearNuc_${Ref}_Nucshift_0.bed
    rm ${DIR}/${Ref}_nearNuc_original_phase1_1.bed ${DIR}/${Ref}_nearNuc_original_phase2_2.bed ${DIR}/${Ref}_nearNuc_original_phase3_3.bed ${DIR}/${Ref}_nearNuc_original_phase4_4.bed ${DIR}/${Ref}_nearNuc_original_phase5_5.bed ${DIR}/${Ref}_nearNuc_original_phase6_6.bed ${DIR}/${Ref}_nearNuc_original_phase7_7.bed ${DIR}/${Ref}_nearNuc_original_phase8_8.bed ${DIR}/${Ref}_nearNuc_original_phase9_9.bed ${DIR}/${Ref}_nearNuc_original_phase0_0.bed ${DIR}/${Ref}_nearNuc_original_phaseno10x.bed

    cat ${DIR}/${Ref}_nearNuc_original_phase1_2.bed ${DIR}/${Ref}_nearNuc_original_phase2_3.bed ${DIR}/${Ref}_nearNuc_original_phase3_4.bed ${DIR}/${Ref}_nearNuc_original_phase4_5.bed ${DIR}/${Ref}_nearNuc_original_phase5_6.bed ${DIR}/${Ref}_nearNuc_original_phase6_7.bed ${DIR}/${Ref}_nearNuc_original_phase7_8.bed ${DIR}/${Ref}_nearNuc_original_phase8_9.bed ${DIR}/${Ref}_nearNuc_original_phase9_0.bed ${DIR}/${Ref}_nearNuc_original_phase0_1.bed | \
    awk '{OFS="\t"} {print $1,$8,$8,$1"_"$8"_"$8,$9,$6,$1,$2,$3,$4,$5,$6,$7}' |  bedtools shift -i - -g $Genome -p -1 -m 1 | awk '{OFS="\t"} {print $1,$2,$3,$1"_"$2"_"$3,$5,$6,$7,$8,$9,$10,$11,$12,$13}'  > $DIR/nearNuc_${Ref}_Nucshift_up1.bed
    rm ${DIR}/${Ref}_nearNuc_original_phase1_2.bed ${DIR}/${Ref}_nearNuc_original_phase2_3.bed ${DIR}/${Ref}_nearNuc_original_phase3_4.bed ${DIR}/${Ref}_nearNuc_original_phase4_5.bed ${DIR}/${Ref}_nearNuc_original_phase5_6.bed ${DIR}/${Ref}_nearNuc_original_phase6_7.bed ${DIR}/${Ref}_nearNuc_original_phase7_8.bed ${DIR}/${Ref}_nearNuc_original_phase8_9.bed ${DIR}/${Ref}_nearNuc_original_phase9_0.bed ${DIR}/${Ref}_nearNuc_original_phase0_1.bed

    cat ${DIR}/${Ref}_nearNuc_original_phase1_3.bed ${DIR}/${Ref}_nearNuc_original_phase2_4.bed ${DIR}/${Ref}_nearNuc_original_phase3_5.bed ${DIR}/${Ref}_nearNuc_original_phase4_6.bed ${DIR}/${Ref}_nearNuc_original_phase5_7.bed ${DIR}/${Ref}_nearNuc_original_phase6_8.bed ${DIR}/${Ref}_nearNuc_original_phase7_9.bed ${DIR}/${Ref}_nearNuc_original_phase8_0.bed ${DIR}/${Ref}_nearNuc_original_phase9_1.bed ${DIR}/${Ref}_nearNuc_original_phase0_2.bed | \
    awk '{OFS="\t"} {print $1,$8,$8,$1"_"$8"_"$8,$9,$6,$1,$2,$3,$4,$5,$6,$7}' |  bedtools shift -i - -g $Genome -p -2 -m 2 | awk '{OFS="\t"} {print $1,$2,$3,$1"_"$2"_"$3,$5,$6,$7,$8,$9,$10,$11,$12,$13}'  > $DIR/nearNuc_${Ref}_Nucshift_up2.bed
    rm ${DIR}/${Ref}_nearNuc_original_phase1_3.bed ${DIR}/${Ref}_nearNuc_original_phase2_4.bed ${DIR}/${Ref}_nearNuc_original_phase3_5.bed ${DIR}/${Ref}_nearNuc_original_phase4_6.bed ${DIR}/${Ref}_nearNuc_original_phase5_7.bed ${DIR}/${Ref}_nearNuc_original_phase6_8.bed ${DIR}/${Ref}_nearNuc_original_phase7_9.bed ${DIR}/${Ref}_nearNuc_original_phase8_0.bed ${DIR}/${Ref}_nearNuc_original_phase9_1.bed ${DIR}/${Ref}_nearNuc_original_phase0_2.bed

    cat ${DIR}/${Ref}_nearNuc_original_phase1_4.bed ${DIR}/${Ref}_nearNuc_original_phase2_5.bed ${DIR}/${Ref}_nearNuc_original_phase3_6.bed ${DIR}/${Ref}_nearNuc_original_phase4_7.bed ${DIR}/${Ref}_nearNuc_original_phase5_8.bed ${DIR}/${Ref}_nearNuc_original_phase6_9.bed ${DIR}/${Ref}_nearNuc_original_phase7_0.bed ${DIR}/${Ref}_nearNuc_original_phase8_1.bed ${DIR}/${Ref}_nearNuc_original_phase9_2.bed ${DIR}/${Ref}_nearNuc_original_phase0_3.bed | \
    awk '{OFS="\t"} {print $1,$8,$8,$1"_"$8"_"$8,$9,$6,$1,$2,$3,$4,$5,$6,$7}' |  bedtools shift -i - -g $Genome -p -3 -m 3 | awk '{OFS="\t"} {print $1,$2,$3,$1"_"$2"_"$3,$5,$6,$7,$8,$9,$10,$11,$12,$13}'  > $DIR/nearNuc_${Ref}_Nucshift_up3.bed
    rm ${DIR}/${Ref}_nearNuc_original_phase1_4.bed ${DIR}/${Ref}_nearNuc_original_phase2_5.bed ${DIR}/${Ref}_nearNuc_original_phase3_6.bed ${DIR}/${Ref}_nearNuc_original_phase4_7.bed ${DIR}/${Ref}_nearNuc_original_phase5_8.bed ${DIR}/${Ref}_nearNuc_original_phase6_9.bed ${DIR}/${Ref}_nearNuc_original_phase7_0.bed ${DIR}/${Ref}_nearNuc_original_phase8_1.bed ${DIR}/${Ref}_nearNuc_original_phase9_2.bed ${DIR}/${Ref}_nearNuc_original_phase0_3.bed 

    cat ${DIR}/${Ref}_nearNuc_original_phase1_5.bed ${DIR}/${Ref}_nearNuc_original_phase2_6.bed ${DIR}/${Ref}_nearNuc_original_phase3_7.bed ${DIR}/${Ref}_nearNuc_original_phase4_8.bed ${DIR}/${Ref}_nearNuc_original_phase5_9.bed ${DIR}/${Ref}_nearNuc_original_phase6_0.bed ${DIR}/${Ref}_nearNuc_original_phase7_1.bed ${DIR}/${Ref}_nearNuc_original_phase8_2.bed ${DIR}/${Ref}_nearNuc_original_phase9_3.bed ${DIR}/${Ref}_nearNuc_original_phase0_4.bed | \
    awk '{OFS="\t"} {print $1,$8,$8,$1"_"$8"_"$8,$9,$6,$1,$2,$3,$4,$5,$6,$7}' |  bedtools shift -i - -g $Genome -p -4 -m 4 | awk '{OFS="\t"} {print $1,$2,$3,$1"_"$2"_"$3,$5,$6,$7,$8,$9,$10,$11,$12,$13}'  > $DIR/nearNuc_${Ref}_Nucshift_up4.bed
    rm ${DIR}/${Ref}_nearNuc_original_phase1_5.bed ${DIR}/${Ref}_nearNuc_original_phase2_6.bed ${DIR}/${Ref}_nearNuc_original_phase3_7.bed ${DIR}/${Ref}_nearNuc_original_phase4_8.bed ${DIR}/${Ref}_nearNuc_original_phase5_9.bed ${DIR}/${Ref}_nearNuc_original_phase6_0.bed ${DIR}/${Ref}_nearNuc_original_phase7_1.bed ${DIR}/${Ref}_nearNuc_original_phase8_2.bed ${DIR}/${Ref}_nearNuc_original_phase9_3.bed ${DIR}/${Ref}_nearNuc_original_phase0_4.bed

    cat ${DIR}/${Ref}_nearNuc_original_phase1_6.bed ${DIR}/${Ref}_nearNuc_original_phase2_7.bed ${DIR}/${Ref}_nearNuc_original_phase3_8.bed ${DIR}/${Ref}_nearNuc_original_phase4_9.bed ${DIR}/${Ref}_nearNuc_original_phase5_0.bed ${DIR}/${Ref}_nearNuc_original_phase6_1.bed ${DIR}/${Ref}_nearNuc_original_phase7_2.bed ${DIR}/${Ref}_nearNuc_original_phase8_3.bed ${DIR}/${Ref}_nearNuc_original_phase9_4.bed ${DIR}/${Ref}_nearNuc_original_phase0_5.bed | \
    awk '{OFS="\t"} {print $1,$8,$8,$1"_"$8"_"$8,$9,$6,$1,$2,$3,$4,$5,$6,$7}' |  bedtools shift -i - -g $Genome -p -5 -m 5 | awk '{OFS="\t"} {print $1,$2,$3,$1"_"$2"_"$3,$5,$6,$7,$8,$9,$10,$11,$12,$13}' > $DIR/nearNuc_${Ref}_Nucshift_up5.bed
    rm ${DIR}/${Ref}_nearNuc_original_phase1_6.bed ${DIR}/${Ref}_nearNuc_original_phase2_7.bed ${DIR}/${Ref}_nearNuc_original_phase3_8.bed ${DIR}/${Ref}_nearNuc_original_phase4_9.bed ${DIR}/${Ref}_nearNuc_original_phase5_0.bed ${DIR}/${Ref}_nearNuc_original_phase6_1.bed ${DIR}/${Ref}_nearNuc_original_phase7_2.bed ${DIR}/${Ref}_nearNuc_original_phase8_3.bed ${DIR}/${Ref}_nearNuc_original_phase9_4.bed ${DIR}/${Ref}_nearNuc_original_phase0_5.bed

    cat ${DIR}/${Ref}_nearNuc_original_phase1_7.bed ${DIR}/${Ref}_nearNuc_original_phase2_8.bed ${DIR}/${Ref}_nearNuc_original_phase3_9.bed ${DIR}/${Ref}_nearNuc_original_phase4_0.bed ${DIR}/${Ref}_nearNuc_original_phase5_1.bed ${DIR}/${Ref}_nearNuc_original_phase6_2.bed ${DIR}/${Ref}_nearNuc_original_phase7_3.bed ${DIR}/${Ref}_nearNuc_original_phase8_4.bed ${DIR}/${Ref}_nearNuc_original_phase9_5.bed ${DIR}/${Ref}_nearNuc_original_phase0_6.bed | \
    awk '{OFS="\t"} {print $1,$8,$8,$1"_"$8"_"$8,$9,$6,$1,$2,$3,$4,$5,$6,$7}' |  bedtools shift -i - -g $Genome -p 4 -m -4 | awk '{OFS="\t"} {print $1,$2,$3,$1"_"$2"_"$3,$5,$6,$7,$8,$9,$10,$11,$12,$13}' > $DIR/nearNuc_${Ref}_Nucshift_down4.bed
    rm ${DIR}/${Ref}_nearNuc_original_phase1_7.bed ${DIR}/${Ref}_nearNuc_original_phase2_8.bed ${DIR}/${Ref}_nearNuc_original_phase3_9.bed ${DIR}/${Ref}_nearNuc_original_phase4_0.bed ${DIR}/${Ref}_nearNuc_original_phase5_1.bed ${DIR}/${Ref}_nearNuc_original_phase6_2.bed ${DIR}/${Ref}_nearNuc_original_phase7_3.bed ${DIR}/${Ref}_nearNuc_original_phase8_4.bed ${DIR}/${Ref}_nearNuc_original_phase9_5.bed ${DIR}/${Ref}_nearNuc_original_phase0_6.bed

    cat ${DIR}/${Ref}_nearNuc_original_phase1_8.bed ${DIR}/${Ref}_nearNuc_original_phase2_9.bed ${DIR}/${Ref}_nearNuc_original_phase3_0.bed ${DIR}/${Ref}_nearNuc_original_phase4_1.bed ${DIR}/${Ref}_nearNuc_original_phase5_2.bed ${DIR}/${Ref}_nearNuc_original_phase6_3.bed ${DIR}/${Ref}_nearNuc_original_phase7_4.bed ${DIR}/${Ref}_nearNuc_original_phase8_5.bed ${DIR}/${Ref}_nearNuc_original_phase9_6.bed ${DIR}/${Ref}_nearNuc_original_phase0_7.bed | \
    awk '{OFS="\t"} {print $1,$8,$8,$1"_"$8"_"$8,$9,$6,$1,$2,$3,$4,$5,$6,$7}' |  bedtools shift -i - -g $Genome -p 3 -m -3 | awk '{OFS="\t"} {print $1,$2,$3,$1"_"$2"_"$3,$5,$6,$7,$8,$9,$10,$11,$12,$13}'  > $DIR/nearNuc_${Ref}_Nucshift_down3.bed
    rm ${DIR}/${Ref}_nearNuc_original_phase1_8.bed ${DIR}/${Ref}_nearNuc_original_phase2_9.bed ${DIR}/${Ref}_nearNuc_original_phase3_0.bed ${DIR}/${Ref}_nearNuc_original_phase4_1.bed ${DIR}/${Ref}_nearNuc_original_phase5_2.bed ${DIR}/${Ref}_nearNuc_original_phase6_3.bed ${DIR}/${Ref}_nearNuc_original_phase7_4.bed ${DIR}/${Ref}_nearNuc_original_phase8_5.bed ${DIR}/${Ref}_nearNuc_original_phase9_6.bed ${DIR}/${Ref}_nearNuc_original_phase0_7.bed

    cat ${DIR}/${Ref}_nearNuc_original_phase1_9.bed ${DIR}/${Ref}_nearNuc_original_phase2_0.bed ${DIR}/${Ref}_nearNuc_original_phase3_1.bed ${DIR}/${Ref}_nearNuc_original_phase4_2.bed ${DIR}/${Ref}_nearNuc_original_phase5_3.bed ${DIR}/${Ref}_nearNuc_original_phase6_4.bed ${DIR}/${Ref}_nearNuc_original_phase7_5.bed ${DIR}/${Ref}_nearNuc_original_phase8_6.bed ${DIR}/${Ref}_nearNuc_original_phase9_7.bed ${DIR}/${Ref}_nearNuc_original_phase0_8.bed | \
    awk '{OFS="\t"} {print $1,$8,$8,$1"_"$8"_"$8,$9,$6,$1,$2,$3,$4,$5,$6,$7}' |  bedtools shift -i - -g $Genome -p 2 -m -2 | awk '{OFS="\t"} {print $1,$2,$3,$1"_"$2"_"$3,$5,$6,$7,$8,$9,$10,$11,$12,$13}' > $DIR/nearNuc_${Ref}_Nucshift_down2.bed
    rm ${DIR}/${Ref}_nearNuc_original_phase1_9.bed ${DIR}/${Ref}_nearNuc_original_phase2_0.bed ${DIR}/${Ref}_nearNuc_original_phase3_1.bed ${DIR}/${Ref}_nearNuc_original_phase4_2.bed ${DIR}/${Ref}_nearNuc_original_phase5_3.bed ${DIR}/${Ref}_nearNuc_original_phase6_4.bed ${DIR}/${Ref}_nearNuc_original_phase7_5.bed ${DIR}/${Ref}_nearNuc_original_phase8_6.bed ${DIR}/${Ref}_nearNuc_original_phase9_7.bed ${DIR}/${Ref}_nearNuc_original_phase0_8.bed

    cat ${DIR}/${Ref}_nearNuc_original_phase1_0.bed ${DIR}/${Ref}_nearNuc_original_phase2_1.bed ${DIR}/${Ref}_nearNuc_original_phase3_2.bed ${DIR}/${Ref}_nearNuc_original_phase4_3.bed ${DIR}/${Ref}_nearNuc_original_phase5_4.bed ${DIR}/${Ref}_nearNuc_original_phase6_5.bed ${DIR}/${Ref}_nearNuc_original_phase7_6.bed ${DIR}/${Ref}_nearNuc_original_phase8_7.bed ${DIR}/${Ref}_nearNuc_original_phase9_8.bed ${DIR}/${Ref}_nearNuc_original_phase0_9.bed | \
    awk '{OFS="\t"} {print $1,$8,$8,$1"_"$8"_"$8,$9,$6,$1,$2,$3,$4,$5,$6,$7}' |  bedtools shift -i - -g $Genome -p 1 -m -1 | awk '{OFS="\t"} {print $1,$2,$3,$1"_"$2"_"$3,$5,$6,$7,$8,$9,$10,$11,$12,$13}' > $DIR/nearNuc_${Ref}_Nucshift_down1.bed
    rm ${DIR}/${Ref}_nearNuc_original_phase1_0.bed ${DIR}/${Ref}_nearNuc_original_phase2_1.bed ${DIR}/${Ref}_nearNuc_original_phase3_2.bed ${DIR}/${Ref}_nearNuc_original_phase4_3.bed ${DIR}/${Ref}_nearNuc_original_phase5_4.bed ${DIR}/${Ref}_nearNuc_original_phase6_5.bed ${DIR}/${Ref}_nearNuc_original_phase7_6.bed ${DIR}/${Ref}_nearNuc_original_phase8_7.bed ${DIR}/${Ref}_nearNuc_original_phase9_8.bed ${DIR}/${Ref}_nearNuc_original_phase0_9.bed

    rm $DIR/${Ref}_nearNuc_original_phase*.bed

    cat $DIR/nearNuc_${Ref}_Nucshift_*.bed | bedtools sort -i | uniq  | sort -k5,5nr | \
    awk -v DIR=$DIR -v BASE="AdjustnearNuc_${Ref}" '{
    if ($13 ~ /Phase0/) {
        print $0 > (DIR "/" BASE "_phase0.bed");
    }
    else if ($13 ~ /Phase1/) {
        print $0 > (DIR "/" BASE "_phase1.bed");
    }
    else if ($13 ~ /Phase2/) {
        print $0 > (DIR "/" BASE "_phase2.bed");
    }
    else if ($13 ~ /Phase3/) {
        print $0 > (DIR "/" BASE "_phase3.bed");
    }
    else if ($13 ~ /Phase4/) {
        print $0 > (DIR "/" BASE "_phase4.bed");
    }
    else if ($13 ~ /Phase5/) {
        print $0 > (DIR "/" BASE "_phase5.bed");
    }
    else if ($13 ~ /Phase6/) {
        print $0 > (DIR "/" BASE "_phase6.bed");
    }
    else if ($13 ~ /Phase7/) {
        print $0 > (DIR "/" BASE "_phase7.bed");
    }
    else if ($13 ~ /Phase8/) {
        print $0 > (DIR "/" BASE "_phase8.bed");
    }
    else if ($13 ~ /Phase9/) {
        print $0 > (DIR "/" BASE "_phase9.bed");
    }
    else {
        print $0 > (DIR "/" BASE "_phaseno10x.bed");
    }
    }'  

    cat $DIR/AdjustnearNuc_${Ref}_phase0.bed $DIR/AdjustnearNuc_${Ref}_phase1.bed $DIR/AdjustnearNuc_${Ref}_phase2.bed $DIR/AdjustnearNuc_${Ref}_phase3.bed $DIR/AdjustnearNuc_${Ref}_phase4.bed $DIR/AdjustnearNuc_${Ref}_phase5.bed $DIR/AdjustnearNuc_${Ref}_phase6.bed $DIR/AdjustnearNuc_${Ref}_phase7.bed $DIR/AdjustnearNuc_${Ref}_phase8.bed $DIR/AdjustnearNuc_${Ref}_phase9.bed $DIR/AdjustnearNuc_${Ref}_phaseno10x.bed > $DIR/AdjustnearNuc_${Ref}_all.bed
    mv $DIR/AdjustnearNuc_${Ref}_phase0.bed $DIR/AdjustnearNuc_${Ref}_phase1.bed $DIR/AdjustnearNuc_${Ref}_phase2.bed $DIR/AdjustnearNuc_${Ref}_phase3.bed $DIR/AdjustnearNuc_${Ref}_phase4.bed $DIR/AdjustnearNuc_${Ref}_phase5.bed $DIR/AdjustnearNuc_${Ref}_phase6.bed $DIR/AdjustnearNuc_${Ref}_phase7.bed $DIR/AdjustnearNuc_${Ref}_phase8.bed $DIR/AdjustnearNuc_${Ref}_phase9.bed $DIR/AdjustnearNuc_${Ref}_phaseno10x.bed $DIR/temp
    rm $DIR/nearNuc_${Ref}_Nucshift_*.bed
    for file in $DIR/Adjust*Nuc_${Ref}_all.bed ; do
        filename=`basename $file ".bed"`
        java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 250 $file -o $DIR/${filename}_250bp.bed
        java -jar $SCRIPTMANAGER read-analysis tag-pileup $DIR/${filename}_250bp.bed $BAMFILE --cpu 4 -5 -1 -M $DIR/SCORES/${BAM}_${filename}_read1
        java -jar $SCRIPTMANAGER figure-generation heatmap -p .95 --blue $DIR/SCORES/${BAM}_${filename}_read1_sense.cdt -o  $DIR/${BAM}_${filename}_read1_sense.png
        java -jar $SCRIPTMANAGER figure-generation heatmap -p .95 --red $DIR/SCORES/${BAM}_${filename}_read1_anti.cdt  -o $DIR/${BAM}_${filename}_read1_anti.png
        java -jar $SCRIPTMANAGER figure-generation merge-heatmap $DIR/${BAM}_${filename}_read1_sense.png $DIR/${BAM}_${filename}_read1_anti.png -o $DIR/${BAM}_${filename}_read1_merge.png
        java -jar $SCRIPTMANAGER figure-generation label-heatmap $DIR/${BAM}_${filename}_read1_merge.png -f 20 -l -125 -m 0 -r +125  -o $DIR/${BAM}_${filename}_read1_merge.svg
        java -jar $SCRIPTMANAGER figure-generation label-heatmap $DIR/${BAM}_${filename}_read1_sense.png -f 20 -l -125 -m 0 -r +125  -o $DIR/${BAM}_${filename}_read1_sense.svg
        java -jar $SCRIPTMANAGER figure-generation label-heatmap $DIR/${BAM}_${filename}_read1_anti.png -f 20 -l -125 -m 0 -r +125  -o $DIR/${BAM}_${filename}_read1_anti.svg
        rm  $DIR/${BAM}_${filename}_read1_sense.png $DIR/${BAM}_${filename}_read1_anti.png  $DIR/${BAM}_${filename}_read1_merge.png
        rm $DIR/${filename}_250bp.bed
    done


mkdir -p $DIR/10xplot

for file in $DIR/10xCDT/*_sense.cdt ; do
    filename=`basename $file "_sense.cdt"`
    perl $COMPOSITE $DIR/10xCDT/${filename}_sense.cdt $DIR/${filename}_sense
    perl $COMPOSITE $DIR/10xCDT/${filename}_anti.cdt $DIR/${filename}_anti
    tail -1 $DIR/${filename}_sense | cat $DIR/${filename}_anti - > $DIR/10xplot/${filename}.out
    rm $DIR/${filename}_sense $DIR/${filename}_anti
done

for file in $DIR/temp/Adjust*_${Ref}_phase*.bed  ; do
    filename=`basename $file ".bed"`
    java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $file -o $DIR/${filenam}_1000bp.bed
    java -jar $SCRIPTMANAGER read-analysis tag-pileup  $DIR/${filenam}_1000bp.bed $BAMFILE --cpu 4 -5 -1 -M $DIR/10xCDT/${BAM}_${filename}_read1 -o $DIR/10xplot/${BAM}_${filename}_read1.out
    rm $DIR/${filenam}_1000bp.bed
done

mkdir -p $DIR/10Dinucleotide

for file in  $DIR/temp/${Ref}*.bed  ; do
  filename=$(basename "$file" ".bed")
  java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $file -o $DIR/10Dinucleotide/${filename}_1000bp.bed
  java -jar $SCRIPTMANAGER sequence-analysis fasta-extract $GENOME   $DIR/10Dinucleotide/${filename}_1000bp.bed -o $DIR/10Dinucleotide/${filename}_1000bp.fa
  rm $DIR/10Dinucleotide/${filename}_1000bp.bed
done

for file in $DIR/temp/${Ref}*.bed  ; do
  filename=$(basename "$file" ".bed")
  python $MOTIFSCAN -i  $DIR/10Dinucleotide/${filename}_1000bp.fa -m RR -o $DIR/10xCDT/RR_${filename}
  python $MOTIFSCAN -i  $DIR/10Dinucleotide/${filename}_1000bp.fa -m YY -o $DIR/10xCDT/YY_${filename}
  python $MOTIFSCAN -i  $DIR/10Dinucleotide/${filename}_1000bp.fa -m WW -o $DIR/10xCDT/WW_${filename}
  python $MOTIFSCAN -i  $DIR/10Dinucleotide/${filename}_1000bp.fa -m SS -o $DIR/10xCDT/SS_${filename}
done

for file in  $DIR/temp/${Ref}*.bed ; do
  filename=$(basename "$file" ".bed")
   perl $COMPOSITE $DIR/10xCDT/RR_${filename}_sense.cdt $DIR/10xplot/RR_${filename}.out
   perl $COMPOSITE $DIR/10xCDT/YY_${filename}_sense.cdt $DIR/10xplot/YY_${filename}.out
   perl $COMPOSITE $DIR/10xCDT/SS_${filename}_sense.cdt $DIR/10xplot/SS_${filename}.out
   perl $COMPOSITE $DIR/10xCDT/WW_${filename}_sense.cdt $DIR/10xplot/WW_${filename}.out
done


#Statiscal
mkdir -p $DIR/${Ref}_BZ
for file in $DIR/${Ref}_*_original_all.bed ; do
    Nuc=`basename $file".bed" | cut -d "_" -f 2`
    cat $DIR/${Ref}_${Nuc}_original_all.bed | awk '{OFS="\t"} {print $4,$7,$8}' |  awk -v DIR=$DIR -v BASE="${Ref}_${Nuc}_original" '{
    if ($2 ~ /Phase0/) {
        print $0 > (DIR "/" BASE "_phase0.csv");
    }
    else if ($2 ~ /Phase1/) {
        print $0 > (DIR "/" BASE "_phase1.csv");
    }
    else if ($2 ~ /Phase2/) {
        print $0 > (DIR "/" BASE "_phase2.csv");
    }
    else if ($2 ~ /Phase3/) {
        print $0 > (DIR "/" BASE "_phase3.csv");
    }
    else if ($2 ~ /Phase4/) {
        print $0 > (DIR "/" BASE "_phase4.csv");
    }
    else if ($2 ~ /Phase5/) {
        print $0 > (DIR "/" BASE "_phase5.csv");
    }
    else if ($2 ~ /Phase6/) {
        print $0 > (DIR "/" BASE "_phase6.csv");
    }
    else if ($2 ~ /Phase7/) {
        print $0 > (DIR "/" BASE "_phase7.csv");
    }
    else if ($2 ~ /Phase8/) {
        print $0 > (DIR "/" BASE "_phase8.csv");
    }
    else if ($2 ~ /Phase9/) {
        print $0 > (DIR "/" BASE "_phase9.csv");
    }
    else {
        print $0 > (DIR "/" BASE "_phaseno10x.csv");
    }
    }'  

    cat $DIR/${Ref}_${Nuc}_shuffle_all.bed | awk '{OFS="\t"} {print $4,$7,$8}' |  awk -v DIR=$DIR -v BASE="${Ref}_${Nuc}_shuffle" '{
    if ($2 ~ /Phase0/) {
        print $0 > (DIR "/" BASE "_phase0.csv");
    }
    else if ($2 ~ /Phase1/) {
        print $0 > (DIR "/" BASE "_phase1.csv");
    }
    else if ($2 ~ /Phase2/) {
        print $0 > (DIR "/" BASE "_phase2.csv");
    }
    else if ($2 ~ /Phase3/) {
        print $0 > (DIR "/" BASE "_phase3.csv");
    }
    else if ($2 ~ /Phase4/) {
        print $0 > (DIR "/" BASE "_phase4.csv");
    }
    else if ($2 ~ /Phase5/) {
        print $0 > (DIR "/" BASE "_phase5.csv");
    }
    else if ($2 ~ /Phase6/) {
        print $0 > (DIR "/" BASE "_phase6.csv");
    }
    else if ($2 ~ /Phase7/) {
        print $0 > (DIR "/" BASE "_phase7.csv");
    }
    else if ($2 ~ /Phase8/) {
        print $0 > (DIR "/" BASE "_phase8.csv");
    }
    else if ($2 ~ /Phase9/) {
        print $0 > (DIR "/" BASE "_phase9.csv");
    }
    else {
        print $0 > (DIR "/" BASE "_phaseno10x.csv");
    }
    }'  

    for X in 1 2 3 4 5 6 7 8 9 0 no10x ; do
        echo -e "Region,Nucleosomephase,enrichment" > $DIR/${Ref}_BZ/${Ref}_${Nuc}_phase${X}_score_peak.csv
        cat $DIR/${Ref}_${Nuc}_*_phase${X}.csv | awk 'BEGIN {OFS=","} { print $1, $2, $3 }' >> $DIR/${Ref}_BZ/${Ref}_${Nuc}_phase${X}_score_peak.csv
        rm $DIR/${Ref}_${Nuc}_*_phase${X}.csv
    done
done

mkdir -p $WRK/Library/10phase_10bp
for folder in $DIR/${Ref}_BZ ; do
    if [ -d "$folder" ]; then
        X="$folder"
        output_file="${X}.out"
        python "$originalshuffleTtest" "$X" "$output_file"
        mv "$output_file" $WRK/Library/10phase_10bp
    fi
done


for file in  $DIR/temp/${Ref}_*_original_phase_*.bed  ; do
        filename=$(basename "$file" ".bed")
        cat $DIR/10xCDT/RR_${filename}_sense.cdt | cut -f  1-2  > $DIR/RR_${filename}_sense_ref.cdt
        cat $DIR/10xCDT/RR_${filename}_sense.cdt | cut -f  383-482,523-622  | paste $DIR/RR_${filename}_sense_ref.cdt - > $DIR/RR_${filename}_383-482-523-622_sense.cdt
        rm $DIR/RR_${filename}_sense_ref.cdt
        cat $DIR/10xCDT/YY_${filename}_sense.cdt | cut -f  1-2  > $DIR/YY_${filename}_sense_ref.cdt
        cat $DIR/10xCDT/YY_${filename}_sense.cdt | cut -f  383-482,523-622  | paste $DIR/YY_${filename}_sense_ref.cdt - > $DIR/YY_${filename}_383-482-523-622_sense.cdt
        rm $DIR/YY_${filename}_sense_ref.cdt
        cat $DIR/10xCDT/SS_${filename}_sense.cdt | cut -f  1-2  > $DIR/SS_${filename}_sense_ref.cdt
        cat $DIR/10xCDT/SS_${filename}_sense.cdt | cut -f  383-482,523-622  | paste $DIR/SS_${filename}_sense_ref.cdt - > $DIR/SS_${filename}_383-482-523-622_sense.cdt
        rm $DIR/SS_${filename}_sense_ref.cdt
        cat $DIR/10xCDT/WW_${filename}_sense.cdt | cut -f  1-2  > $DIR/WW_${filename}_sense_ref.cdt
        cat $DIR/10xCDT/WW_${filename}_sense.cdt | cut -f  383-482,523-622  | paste $DIR/WW_${filename}_sense_ref.cdt - > $DIR/WW_${filename}_383-482-523-622_sense.cdt
        rm $DIR/WW_${filename}_sense_ref.cdt
done

mkdir -p $DIR/${Ref}_Dinucleotide

mv $DIR/RR_*.cdt  $DIR/WW_*.cdt $DIR/YY_*.cdt  $DIR/SS_*.cdt $DIR/${Ref}_Dinucleotide

for file in $DIR/${Ref}_Dinucleotide/*_sense.cdt ; do
        filename=$(basename "$file" ".cdt")
        java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum $file -o $DIR/${Ref}_Dinucleotide/${filename}_SCORES.out 
        tail -n +2 $DIR/${Ref}_Dinucleotide/${filename}_SCORES.out  > $DIR/${Ref}_Dinucleotide/${filename}_score.out
      
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
            }'  | paste $DIR/${Ref}_Dinucleotide/${filename}_score.out  - |  awk '{
            OFS="\t";
            print $1,$2,$3/($2+1),$4/($2+1),$5/($2+1),$6/($2+1),$7/($2+1),$8/($2+1),$9/($2+1),$10/($2+1),$11/($2+1),$12/($2+1)}' > $DIR/${Ref}_Dinucleotide/${filename}_score_peak.out
        rm $DIR/${Ref}_Dinucleotide/${filename}_score.out
        rm $DIR/${Ref}_Dinucleotide/${filename}_SCORES.out
done

for folder in $DIR/${Ref}_Dinucleotide ; do
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


for folder in  $DIR/${Ref}_Dinucleotide ; do
    if [ -d "$folder" ]; then
        X="$folder"
        output_file="${X}.out"
        python "$chisquare" "$X" "$output_file"
        mv "$output_file" $WRK/Library/10phase_10bp
    fi
done

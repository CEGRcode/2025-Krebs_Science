#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=48gb
#SBATCH -t 6:00:00
#SBATCH -A open
#SBATCH -o logs/Benzonase-seq_download-merge-dedup-merge.log.out
#SBATCH -e logs/Benzonase-seq_download-merge-dedup-merge.log.err

# Script to hardcode the merging/renaming of PEGR BAM files into a standard file naming system

# data/BAM
#   |--BNase-seq_50U-3min_merge_hg38.bam
#   |--BNase-seq_50U-30min_merge_hg38.bam
#   |--BNase-seq_50U-10min_merge_hg38.bam

# Dependencies
# - java
# - samtools

set -exo
source activate bx

[ -d Benzonase ] || mkdir Benzonase

cd Benzonase

# Script shortcuts
PICARD=../../bin/picard.jar
BAMDIR=../../data/BAM

# ===== Standardize bio-rep filenames (merge & markdup as needed) ======

# Naked DNA - 50U titration
mv 40726_Input_-_K562_-_IMDM_-_BI.bam $BAMDIR/NakedDNA_BNase-seq_0.04U-10min-BI_hg38.bam
mv 40727_Input_-_K562_-_IMDM_-_BI.bam $BAMDIR/NakedDNA_BNase-seq_0.125U-10min-BI_hg38.bam
mv 40728_Input_-_K562_-_IMDM_-_BI.bam $BAMDIR/NakedDNA_BNase-seq_0.375U-10min-BI_hg38.bam

# SEM cells (+IAA/+DMSO)
mv 41859_Input_-_SEM_-_RPMI_-_BI_hg38_hg38.bam $BAMDIR/SEM-CTCFAID2-DMSO_-_BI_rep1_hg38.bam
mv 41860_Input_-_SEM_-_RPMI_1uM5PhIAA_BI_hg38_hg38.bam $BAMDIR/SEM-CTCFAID2-1uM5PhIAA_-_BI_rep1_hg38.bam
mv 42076_Input_-_SEM_-_RPMI_-_BI_hg38_hg38.bam $BAMDIR/SEM-CTCFAID2-DMSO_-_BI_rep2_hg38.bam
mv 42077_Input_-_SEM_-_RPMI_1uM5PhIAA_BI_hg38_hg38.bam $BAMDIR/SEM-CTCFAID2-1uM5PhIAA_-_BI_rep2_hg38.bam

# 50U - 3min
# Rename non-resequenced replicates to standard naming
mv 28385_Input_-_K562_-_-_50Unuclease3min-0cycSonic_BI.bam $BAMDIR/BNase-seq_50U-3min_1_hg38.bam
mv 28386_Input_-_K562_-_-_50Unuclease3min-0cycSonic_BI.bam $BAMDIR/BNase-seq_50U-3min_2_hg38.bam

# 50U - 10min
# Rename non-resequenced replicates to standard naming
mv 28393_Input_-_K562_-_-_50Unuclease30min-0cycSonic_BI.bam $BAMDIR/BNase-seq_50U-30min_1_hg38.bam
mv 28394_Input_-_K562_-_-_50Unuclease30min-0cycSonic_BI.bam $BAMDIR/BNase-seq_50U-30min_2_hg38.bam

# 50U - 30min (rep 1-11, 18-20)
# Rename non-resequenced replicates to standard naming
mv 25769_Input_-_K562_-_-_50Unuclease-0cycSonic_BI.bam $BAMDIR/BNase-seq_50U-10min_1_hg38.bam
mv 25976_Input_-_K562_-_-_50Unuclease10min-0cycSonic_BI.bam $BAMDIR/BNase-seq_50U-10min_2_hg38.bam
mv 28326_Input_-_K562_-_-_50Unuclease10min-0cycSonic-newQuenchingbuffer_BI.bam BNase-seq_50U-10min_3_hg38.bam
mv 28389_Input_-_K562_-_-_50Unuclease10min-0cycSonic_BI.bam BNase-seq_50U-10min_4_hg38.bam
mv 28390_Input_-_K562_-_-_50Unuclease10min-0cycSonic_BI.bam BNase-seq_50U-10min_5_hg38.bam
mv 28468_Input_-_K562_-_-_50Unuclease10min-0cycSonic-Total_BI.bam BNase-seq_50U-10min_6_hg38.bam
mv 28469_Input_-_K562_-_-_50Unuclease10min-0cycSonic-Pellet_BI.bam BNase-seq_50U-10min_7_hg38.bam
mv 28471_Input_-_K562_-_-_50Unuclease10min-0cycSonic-Super_BI.bam BNase-seq_50U-10min_8_hg38.bam
mv 28472_Input_-_K562_-_-_50Unuclease10min-0cycSonic-Super_BI.bam BNase-seq_50U-10min_9_hg38.bam
mv 28473_Input_-_K562_-_-_50Unuclease10min-0cycSonic-Super_BI.bam BNase-seq_50U-10min_10_hg38.bam
mv 28820_Input_-_K562_-_-_50Unuclease10min-0cycSonic-Super_BI.bam BNase-seq_50U-10min_11_hg38.bam
mv 33646_Input_-_K562_-_-_XL-50Unuclease-0cycSonic-super-EndRepair-Atail-ATligation_BI.bam BNase-seq_50U-10min_18_hg38.bam
mv 33647_Input_-_K562_-_-_XL-50Unuclease-0cycSonic-super-noEndRepair-Atail-ATligation_BI.bam BNase-seq_50U-10min_19_hg38.bam
mv 33648_Input_-_K562_-_-_XL-50Unuclease-0cycSonic-super-EndRepair-Atail-ATligation_BI.bam BNase-seq_50U-10min_20_hg38.bam

# 50U - 30min (biorep 12-17)
# Merge and MarkDuplicates on resequencing replicates
# Mark and remove duplicates for technical replicates pairs and single (non-technical sequencing) replicates
# Optional arguments: VALIDATION_STRINGENCY='LENIENT' QUIET=true VERBOSITY=ERROR
java -jar $PICARD MarkDuplicates -O $BAMDIR/BNase-seq_50U-10min_12_hg38.bam -M BNase-seq_50U-10min_12_hg38.metrics \
	--REMOVE_DUPLICATES 'true' --ASSUME_SORTED 'true' --DUPLICATE_SCORING_STRATEGY 'SUM_OF_BASE_QUALITIES' \
	--READ_NAME_REGEX '[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*.' --OPTICAL_DUPLICATE_PIXEL_DISTANCE '100' \
	-I 33023_Input_-_K562_-_-_XL-50Unuclease-0cycSonic-total-input_BI.bam \
	-I 33159_Input_-_K562_-_-_XL-50Unuclease-0cycSonic-total-input_BI.bam
java -jar $PICARD MarkDuplicates -O $BAMDIR/BNase-seq_50U-10min_13_hg38.bam -M BNase-seq_50U-10min_13_hg38.metrics \
	--REMOVE_DUPLICATES 'true' --ASSUME_SORTED 'true' --DUPLICATE_SCORING_STRATEGY 'SUM_OF_BASE_QUALITIES' \
	--READ_NAME_REGEX '[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*.' --OPTICAL_DUPLICATE_PIXEL_DISTANCE '100' \
	-I 33024_Input_-_K562_-_-_XL-50Unuclease-0cycSonic-total-input_BI.bam \
	-I 33160_Input_-_K562_-_-_XL-50Unuclease-0cycSonic-total-input_BI.bam
java -jar $PICARD MarkDuplicates -O $BAMDIR/BNase-seq_50U-10min_14_hg38.bam -M BNase-seq_50U-10min_14_hg38.metrics \
	--REMOVE_DUPLICATES 'true' --ASSUME_SORTED 'true' --DUPLICATE_SCORING_STRATEGY 'SUM_OF_BASE_QUALITIES' \
	--READ_NAME_REGEX '[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*.' --OPTICAL_DUPLICATE_PIXEL_DISTANCE '100' \
	-I 33025_Input_-_K562_-_-_XL-50Unuclease-0cycSonic-super-input_BI.bam \
	-I 33161_Input_-_K562_-_-_XL-50Unuclease-0cycSonic-super-input_BI.bam
java -jar $PICARD MarkDuplicates -O $BAMDIR/BNase-seq_50U-10min_15_hg38.bam -M BNase-seq_50U-10min_15_hg38.metrics \
	--REMOVE_DUPLICATES 'true' --ASSUME_SORTED 'true' --DUPLICATE_SCORING_STRATEGY 'SUM_OF_BASE_QUALITIES' \
	--READ_NAME_REGEX '[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*.' --OPTICAL_DUPLICATE_PIXEL_DISTANCE '100' \
	-I 33026_Input_-_K562_-_-_XL-50Unuclease-0cycSonic-super-input_BI.bam \
	-I 33162_Input_-_K562_-_-_XL-50Unuclease-0cycSonic-super-input_BI.bam
java -jar $PICARD MarkDuplicates -O $BAMDIR/BNase-seq_50U-10min_16_hg38.bam -M BNase-seq_50U-10min_16_hg38.metrics \
	--REMOVE_DUPLICATES 'true' --ASSUME_SORTED 'true' --DUPLICATE_SCORING_STRATEGY 'SUM_OF_BASE_QUALITIES' \
	--READ_NAME_REGEX '[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*.' --OPTICAL_DUPLICATE_PIXEL_DISTANCE '100' \
	-I 33027_Input_-_K562_-_-_XL-50Unuclease-0cycSonic-pellet-input_BI.bam \
	-I 33163_Input_-_K562_-_-_XL-50Unuclease-0cycSonic-pellet-input_BI.bam
java -jar $PICARD MarkDuplicates -O $BAMDIR/BNase-seq_50U-10min_17_hg38.bam -M BNase-seq_50U-10min_17_hg38.metrics \
	--REMOVE_DUPLICATES 'true' --ASSUME_SORTED 'true' --DUPLICATE_SCORING_STRATEGY 'SUM_OF_BASE_QUALITIES' \
	--READ_NAME_REGEX '[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*.' --OPTICAL_DUPLICATE_PIXEL_DISTANCE '100' \
	-I 33028_Input_-_K562_-_-_XL-50Unuclease-0cycSonic-pellet-input_BI.bam \
	-I 33164_Input_-_K562_-_-_XL-50Unuclease-0cycSonic-pellet-input_BI.bam


# ===== Create merged versions ======

# 50U - 3min
java -jar $PICARD MergeSamFiles -O $BAMDIR/BNase-seq_50U-3min_merge_hg38.bam \
	--USE_THREADING 'true' \
	-I $BAMDIR/BNase-seq_50U-3min_1_hg38.bam -I $BAMDIR/BNase-seq_50U-3min_2_hg38.bam
# 50U - 30min
java -jar $PICARD MergeSamFiles -O $BAMDIR/BNase-seq_50U-30min_merge_hg38.bam \
	--USE_THREADING 'true' \
	-I $BAMDIR/BNase-seq_50U-30min_1_hg38.bam -I $BAMDIR/BNase-seq_50U-30min_2_hg38.bam
# 50U - 10min
java -jar $PICARD MergeSamFiles -O $BAMDIR/BNase-seq_50U-10min_merge_hg38.bam \
	--USE_THREADING 'true' \
	-I $BAMDIR/BNase-seq_50U-10min_1_hg38.bam -I $BAMDIR/BNase-seq_50U-10min_2_hg38.bam \
	-I BNase-seq_50U-10min_3_hg38.bam -I BNase-seq_50U-10min_4_hg38.bam \
	-I BNase-seq_50U-10min_5_hg38.bam -I BNase-seq_50U-10min_6_hg38.bam \
	-I BNase-seq_50U-10min_7_hg38.bam -I BNase-seq_50U-10min_8_hg38.bam \
	-I BNase-seq_50U-10min_9_hg38.bam -I BNase-seq_50U-10min_10_hg38.bam \
	-I BNase-seq_50U-10min_11_hg38.bam -I BNase-seq_50U-10min_12_hg38.bam \
	-I BNase-seq_50U-10min_13_hg38.bam -I BNase-seq_50U-10min_14_hg38.bam \
	-I BNase-seq_50U-10min_15_hg38.bam -I BNase-seq_50U-10min_16_hg38.bam \
	-I BNase-seq_50U-10min_17_hg38.bam -I BNase-seq_50U-10min_18_hg38.bam \
	-I BNase-seq_50U-10min_19_hg38.bam -I BNase-seq_50U-10min_20_hg38.bam

for FILE in $BAMDIR/BNase-seq_*.bam;
do
	samtools index $FILE
done

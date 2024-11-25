#!/bin/bash

# Copy over heatmap and composite data for S7 from Library and generate custom figures.

### CHANGE ME
WRK=/path/to/2024-Krebs_Science/Z_Figures
WRK=/ocean/projects/see180003p/owlang/2024-Krebs_Science/Z_Figures
WRK=/storage/home/owl5022/scratch/2024-Krebs_Science/Z_Figures
THREADS=4
###

# Dependencies
# - java
# - python
# - pandas
# - seaborn

set -exo
module load anaconda3
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Inputs and outputs
LIBRARY=$WRK/../X_Bulk_Processing/Library

# Script shortcuts
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar
VIOLIN=$WRK/../bin/make_violin-split_plot.py
SUMCDT=$WRK/../bin/sum_each_CDT.py

[ -d S7 ] || mkdir S7

# ===============================================================================================================================

[ -d S7/b ] || mkdir S7/b

# Heatmaps
BED=PlusOneDyad_SORT-Expression_2000bp
cp $LIBRARY/$BED/SVG/BNase-ChIP_H2A_merge_hg38_$BED\_midpoint-MIN128-MAX164_TotalTag_combined.svg S7/b
cp $LIBRARY/$BED/SVG/BNase-ChIP_H2AZ_merge_hg38_$BED\_midpoint-MIN128-MAX164_TotalTag_combined.svg S7/b
cp $LIBRARY/$BED/SVG/BNase-ChIP_H2B_merge_hg38_$BED\_midpoint-MIN128-MAX164_TotalTag_combined.svg S7/b
cp $LIBRARY/$BED/SVG/BNase-ChIP_H3_merge_hg38_$BED\_midpoint-MIN128-MAX164_TotalTag_combined.svg S7/b
cp $LIBRARY/$BED/SVG/BNase-ChIP_H3K4me1_merge_hg38_$BED\_midpoint-MIN128-MAX164_TotalTag_combined.svg S7/b
cp $LIBRARY/$BED/SVG/BNase-ChIP_H3K4me3_merge_hg38_$BED\_midpoint-MIN128-MAX164_TotalTag_combined.svg S7/b
cp $LIBRARY/$BED/SVG/BNase-ChIP_H3K9ac_merge_hg38_$BED\_midpoint-MIN128-MAX164_TotalTag_combined.svg S7/b
cp $LIBRARY/$BED/SVG/BNase-ChIP_H3K27ac_merge_hg38_$BED\_midpoint-MIN128-MAX164_TotalTag_combined.svg S7/b
cp $LIBRARY/$BED/SVG/BNase-ChIP_H3K27me3_merge_hg38_$BED\_midpoint-MIN128-MAX164_TotalTag_combined.svg S7/b
cp $LIBRARY/$BED/SVG/BNase-ChIP_H3K36me3_merge_hg38_$BED\_midpoint-MIN128-MAX164_TotalTag_combined.svg S7/b
cp $LIBRARY/$BED/SVG/BNase-ChIP_H4_merge_hg38_$BED\_midpoint-MIN128-MAX164_TotalTag_combined.svg S7/b

# Composites
cp $LIBRARY/$BED/Composites/BNase-ChIP_H2A_merge_hg38_$BED\_midpoint-MIN128-MAX164_TotalTag_combined.out S7/b/H2A.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H2AZ_merge_hg38_$BED\_midpoint-MIN128-MAX164_TotalTag_combined.out S7/b/H2AZ.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H2B_merge_hg38_$BED\_midpoint-MIN128-MAX164_TotalTag_combined.out S7/b/H2B.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3_merge_hg38_$BED\_midpoint-MIN128-MAX164_TotalTag_combined.out S7/b/H3.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3K4me1_merge_hg38_$BED\_midpoint-MIN128-MAX164_TotalTag_combined.out S7/b/H3K4me1.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3K4me3_merge_hg38_$BED\_midpoint-MIN128-MAX164_TotalTag_combined.out S7/b/H3K4me3.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3K9ac_merge_hg38_$BED\_midpoint-MIN128-MAX164_TotalTag_combined.out S7/b/H3K9ac.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3K27ac_merge_hg38_$BED\_midpoint-MIN128-MAX164_TotalTag_combined.out S7/b/H3K27ac.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3K27me3_merge_hg38_$BED\_midpoint-MIN128-MAX164_TotalTag_combined.out S7/b/H3K27me3.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3K36me3_merge_hg38_$BED\_midpoint-MIN128-MAX164_TotalTag_combined.out S7/b/H3K36me3.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H4_merge_hg38_$BED\_midpoint-MIN128-MAX164_TotalTag_combined.out S7/b/H4.out


# ===============================================================================================================================

[ -d S7/c ] || mkdir S7/c

# Heatmaps
BED=PlusOneDyad_SORT-Expression_2000bp
cp $LIBRARY/$BED/SVG/CUTRUN_H2AZ_merge_hg38_$BED\_midpoint-MIN128-MAX164_TotalTag_combined.svg S7/c
cp $LIBRARY/$BED/SVG/CUTRUN_H3K4me1_merge_hg38_$BED\_midpoint-MIN128-MAX164_TotalTag_combined.svg S7/c
cp $LIBRARY/$BED/SVG/CUTRUN_H3K4me3_merge_hg38_$BED\_midpoint-MIN128-MAX164_TotalTag_combined.svg S7/c
cp $LIBRARY/$BED/SVG/CUTRUN_H3K27ac_merge_hg38_$BED\_midpoint-MIN128-MAX164_TotalTag_combined.svg S7/c
cp $LIBRARY/$BED/SVG/CUTRUN_H3K27me3_merge_hg38_$BED\_midpoint-MIN128-MAX164_TotalTag_combined.svg S7/c
cp $LIBRARY/$BED/SVG/CUTRUN_IgG_merge_hg38_$BED\_midpoint-MIN128-MAX164_TotalTag_combined.svg S7/c

# Composites
cp $LIBRARY/$BED/Composites/CUTRUN_H2AZ_merge_hg38_$BED\_midpoint-MIN128-MAX164_TotalTag_combined.out S7/c
cp $LIBRARY/$BED/Composites/CUTRUN_H3K4me1_merge_hg38_$BED\_midpoint-MIN128-MAX164_TotalTag_combined.out S7/c
cp $LIBRARY/$BED/Composites/CUTRUN_H3K4me3_merge_hg38_$BED\_midpoint-MIN128-MAX164_TotalTag_combined.out S7/c
cp $LIBRARY/$BED/Composites/CUTRUN_H3K27ac_merge_hg38_$BED\_midpoint-MIN128-MAX164_TotalTag_combined.out S7/c
cp $LIBRARY/$BED/Composites/CUTRUN_H3K27me3_merge_hg38_$BED\_midpoint-MIN128-MAX164_TotalTag_combined.out S7/c
cp $LIBRARY/$BED/Composites/CUTRUN_IgG_merge_hg38_$BED\_midpoint-MIN128-MAX164_TotalTag_combined.out S7/c


# ===============================================================================================================================

[ -d S7/d ] || mkdir S7/d

# Composites
BED=PlusOneDyad_SORT-Expression_2000bp
cp $LIBRARY/$BED/Composites/BNase-ChIP_H2AZ_merge_hg38_$BED\_5read1-MIN128-MAX164.out S7/d/H2AZ.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H2B_merge_hg38_$BED\_5read1-MIN128-MAX164.out S7/d/H2B.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3_merge_hg38_$BED\_5read1-MIN128-MAX164.out S7/d/H3.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H4_merge_hg38_$BED\_5read1-MIN128-MAX164.out S7/d/H4.out
cp $LIBRARY/$BED/Composites/BNase-seq_50U-10min_merge_hg38_$BED\_midpoint-MIN128-MAX164_TotalTag_combined.out S7/d/BNase.out


# ===============================================================================================================================

[ -d S7/e ] || mkdir S7/e

# (F4d) Calculate density violin plots for active histone mod(s) w/ respepct to their base histone.
# see 04_Figures/Fig_4d.sh

BED=PlusOneDyad_SORT-Expression_2000bp
CDIR=../X_Bulk_Processing/Library/$BED/CDT/
FDIR=../data/BAM/NormalizationFactors

TEMP=temp-S7e
[ -d $TEMP ] || mkdir $TEMP


for TARGET in "H2AZ" "H2B" "H3" "H3K4me3" "H3K9ac" "H3K27ac";
do
	BAM=BNase-ChIP_${TARGET}_merge_hg38
	FACTOR=`grep '^Scaling factor' $FDIR/$BAM\_TotalTag_ScalingFactors.out | awk '{print $3}'`

	# Slice Proximal/Distal tag counts
	cut -f1,2,930-999   $CDIR/$BAM\_$BED\_5read1-MIN128-MAX164_sense.cdt > $TEMP/$TARGET\_Proximal_NOSCALE.cdt
	cut -f1,2,1003-1072 $CDIR/$BAM\_$BED\_5read1-MIN128-MAX164_anti.cdt  > $TEMP/$TARGET\_Distal_NOSCALE.cdt

	# Scale by scaling factor
	java -jar $SCRIPTMANAGER read-analysis scale-matrix -s $FACTOR $TEMP/$TARGET\_Proximal_NOSCALE.cdt -o $TEMP/$TARGET\_Proximal.cdt
	java -jar $SCRIPTMANAGER read-analysis scale-matrix -s $FACTOR $TEMP/$TARGET\_Distal_NOSCALE.cdt -o $TEMP/$TARGET\_Distal.cdt
done

# Merge each (Proximal/Distal) occupancy pair into one tab-delimited file
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -o $TEMP/H2AZ-H2B_Proximal.tab   $TEMP/H2AZ_Proximal.cdt    $TEMP/H2B_Proximal.cdt
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -o $TEMP/H2AZ-H2B_Distal.tab     $TEMP/H2AZ_Distal.cdt      $TEMP/H2B_Distal.cdt
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -o $TEMP/H3K4me3-H3_Proximal.tab $TEMP/H3K4me3_Proximal.cdt $TEMP/H3_Proximal.cdt
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -o $TEMP/H3K4me3-H3_Distal.tab   $TEMP/H3K4me3_Distal.cdt   $TEMP/H3_Distal.cdt
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -o $TEMP/H3K9ac-H3_Proximal.tab  $TEMP/H3K9ac_Proximal.cdt  $TEMP/H3_Proximal.cdt
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -o $TEMP/H3K9ac-H3_Distal.tab    $TEMP/H3K9ac_Distal.cdt    $TEMP/H3_Distal.cdt
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -o $TEMP/H3K27ac-H3_Proximal.tab $TEMP/H3K27ac_Proximal.cdt $TEMP/H3_Proximal.cdt
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -o $TEMP/H3K27ac-H3_Distal.tab   $TEMP/H3K27ac_Distal.cdt   $TEMP/H3_Distal.cdt

# Calculate log2 density at each coordinate given the pair of values
sed '1d' $TEMP/H2AZ-H2B_Proximal.tab   | awk 'BEGIN {OFS="\t"}{z = (log(($2+1)/($3+1))/log(2)); print $1,z,"H2AZ-H2B","Proximal"}'   > $TEMP/H2AZ-H2B_Proximal.density
sed '1d' $TEMP/H2AZ-H2B_Distal.tab     | awk 'BEGIN {OFS="\t"}{z = (log(($2+1)/($3+1))/log(2)); print $1,z,"H2AZ-H2B","Distal"}'     > $TEMP/H2AZ-H2B_Distal.density
sed '1d' $TEMP/H3K4me3-H3_Proximal.tab | awk 'BEGIN {OFS="\t"}{z = (log(($2+1)/($3+1))/log(2)); print $1,z,"H3K4me3-H3","Proximal"}' > $TEMP/H3K4me3-H3_Proximal.density
sed '1d' $TEMP/H3K4me3-H3_Distal.tab   | awk 'BEGIN {OFS="\t"}{z = (log(($2+1)/($3+1))/log(2)); print $1,z,"H3K4me3-H3","Distal"}'   > $TEMP/H3K4me3-H3_Distal.density
sed '1d' $TEMP/H3K9ac-H3_Proximal.tab  | awk 'BEGIN {OFS="\t"}{z = (log(($2+1)/($3+1))/log(2)); print $1,z,"H3K9ac-H3","Proximal"}'  > $TEMP/H3K9ac-H3_Proximal.density
sed '1d' $TEMP/H3K9ac-H3_Distal.tab    | awk 'BEGIN {OFS="\t"}{z = (log(($2+1)/($3+1))/log(2)); print $1,z,"H3K9ac-H3","Distal"}'    > $TEMP/H3K9ac-H3_Distal.density
sed '1d' $TEMP/H3K27ac-H3_Proximal.tab | awk 'BEGIN {OFS="\t"}{z = (log(($2+1)/($3+1))/log(2)); print $1,z,"H3K27ac-H3","Proximal"}' > $TEMP/H3K27ac-H3_Proximal.density
sed '1d' $TEMP/H3K27ac-H3_Distal.tab   | awk 'BEGIN {OFS="\t"}{z = (log(($2+1)/($3+1))/log(2)); print $1,z,"H3K27ac-H3","Distal"}'   > $TEMP/H3K27ac-H3_Distal.density

# Compile density info
cat $TEMP/*.density | gzip -dc > S7/e/DensityInfo.tab

# Generate violin plot
python $VIOLIN -i <(gzip -dc S7/h/DensityInfo.tab.gz | cut -f2,3,4) -o S7/e/DensityInfo.svg \
	--width 4 --height 4 --preset1 \
	--title "Density at +1 nucleosome" \
	--xlabel "modification" --ylabel "Density (log2)"

# Clean-up
rm -r $TEMP


# ===============================================================================================================================

[ -d S7/f ] || mkdir S7/f

# Heatmaps
BED=PlusOneDyad_SORT-Expression_WithUnexpressed_2000bp
cp $LIBRARY/$BED/SVG/CoPRO_Capped_merge_hg38_$BED\_5read2_merge_treeview_label.svg S7/f
cp $LIBRARY/$BED/SVG/CoPRO_Capped_merge_hg38_$BED\_5read1_merge_treeview_label.svg S7/f
cp $LIBRARY/$BED/SVG/BNase-seq_50U-10min_merge_hg38_$BED\_midpoint_TotalTag_combined.svg S7/f

# Composites
cp $LIBRARY/$BED/Composites/CoPRO_Capped_merge_hg38_$BED\_5read2.out S7/f
cp $LIBRARY/$BED/Composites/CoPRO_Capped_merge_hg38_$BED\_5read1.out S7/f
cp $LIBRARY/$BED/Composites/ChIP-exo_Pol2_merge_hg38_$BED\_5read1_TotalTag.out S7/f
cp $LIBRARY/$BED/Composites/BNase-seq_50U-10min_merge_hg38_$BED\_midpoint_TotalTag_combined.out S7/f

# Custom combined matrix Pol2 heatmap
python $SUMCDT -o S7/f/ChIP-exo_Pol2_merge_hg38_$BED\_5read1_TotalTag_combined.cdt \
				-1 $LIBRARY/$BED/CDT/ChIP-exo_Pol2_merge_hg38_$BED\_5read1_TotalTag_sense.cdt \
				-2 $LIBRARY/$BED/CDT/ChIP-exo_Pol2_merge_hg38_$BED\_5read1_TotalTag_anti.cdt
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation heatmap --color 833C0C -p 0.95 \
		S7/f/ChIP-exo_Pol2_merge_hg38_$BED\_5read1_TotalTag_combined.cdt \
		-o S7/f/ChIP-exo_Pol2_merge_hg38_$BED\_5read1_TotalTag_combined.png
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation label-heatmap \
		S7/f/ChIP-exo_Pol2_merge_hg38_$BED\_5read1_TotalTag_combined.png \
		-l "-1" -m "0" -r "+1" -w 1 -f 20 \
		-x "Distance from TSS (kb)" -y "${BED}" \
		-o F7/a/$BASE\_combined.svg

# ===============================================================================================================================

[ -d S7/g ] || mkdir S7/g

# Composites
BED=PlusOneDyad_SORT-pHN-dHN_GROUP-TOP-2500_1000bp
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3K4me3_merge_hg38_$BED\_5read1-MAX80_TotalTag.out S7/g/TOP-H3K4me3.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3K9ac_merge_hg38_$BED\_5read1-MAX80_TotalTag.out S7/g/TOP-H3K9ac.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3K27ac_merge_hg38_$BED\_5read1-MAX80_TotalTag.out S7/g/TOP-H3K27ac.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3K4me3_merge_hg38_$BED\_5read1-MAX80_TotalTag.out S7/g/TOP-H3K4me3.out
cp $LIBRARY/$BED/Composites/BNase-seq_50U-10min_merge_hg38_$BED\_midpoint-MAX80_TotalTag_combined.out S7/g/TOP-BNase.out

BED=PlusOneDyad_SORT-pHN-dHN_GROUP-BOTTOM-2500_1000bp
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3K9ac_merge_hg38_$BED\_5read1-MAX80_TotalTag.out S7/g/BOTTOM-H3K9ac.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3K27ac_merge_hg38_$BED\_5read1-MAX80_TotalTag.out S7/g/BOTTOM-H3K27ac.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3K4me3_merge_hg38_$BED\_5read1-MAX80_TotalTag.out S7/g/BOTTOM-H3K4me3.out
cp $LIBRARY/$BED/Composites/BNase-seq_50U-10min_merge_hg38_$BED\_midpoint-MAX80_TotalTag_combined.out S7/g/BOTTOM-BNase.out


# ===============================================================================================================================

[ -d S7/h ] || mkdir S7/h


# (F6d) Calculate density violin plots for active histone mod(s) w/ respepct to their base histone.
# see 04_Figures/Fig_6d.sh

BED=PlusOneDyad_SORT-Expression_2000bp
CDIR=../X_Bulk_Processing/Library/$BED/CDT

TEMP=temp-S7h
[ -d $TEMP ] || mkdir $TEMP

for TARGET in "H3" "H3K4me3" "H3K9ac" "H3K27ac";
do
	BAM=BNase-ChIP_${TARGET}_merge_hg38

	# Slice Proximal/Distal tag counts
	cut -f1,2,920-999   $CDIR/$BAM\_$BED\_5read1-MAX80_TotalTag_sense.cdt > $TEMP/$TARGET\_Proximal.cdt
	cut -f1,2,1003-1082 $CDIR/$BAM\_$BED\_5read1-MAX80_TotalTag_anti.cdt  > $TEMP/$TARGET\_Distal.cdt

done

# Merge each (Proximal/Distal) occupancy pair into one tab-delimited file
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -o $TEMP/H3K4me3-H3_Proximal.tab $TEMP/H3K4me3_Proximal.cdt $TEMP/H3_Proximal.cdt
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -o $TEMP/H3K4me3-H3_Distal.tab   $TEMP/H3K4me3_Distal.cdt   $TEMP/H3_Distal.cdt
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -o $TEMP/H3K9ac-H3_Proximal.tab  $TEMP/H3K9ac_Proximal.cdt  $TEMP/H3_Proximal.cdt
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -o $TEMP/H3K9ac-H3_Distal.tab    $TEMP/H3K9ac_Distal.cdt    $TEMP/H3_Distal.cdt
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -o $TEMP/H3K27ac-H3_Proximal.tab $TEMP/H3K27ac_Proximal.cdt $TEMP/H3_Proximal.cdt
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -o $TEMP/H3K27ac-H3_Distal.tab   $TEMP/H3K27ac_Distal.cdt   $TEMP/H3_Distal.cdt

# Calculate log2 density at each coordinate given the pair of values
sed '1d' $TEMP/H3K4me3-H3_Proximal.tab | awk 'BEGIN {OFS="\t"}{z = (log(($2+1)/($3+1))/log(2)); print $1,z,"H3K4me3-H3","Proximal"}' > $TEMP/H3K4me3-H3_Proximal.density
sed '1d' $TEMP/H3K4me3-H3_Distal.tab   | awk 'BEGIN {OFS="\t"}{z = (log(($2+1)/($3+1))/log(2)); print $1,z,"H3K4me3-H3","Distal"}'   > $TEMP/H3K4me3-H3_Distal.density
sed '1d' $TEMP/H3K9ac-H3_Proximal.tab  | awk 'BEGIN {OFS="\t"}{z = (log(($2+1)/($3+1))/log(2)); print $1,z,"H3K9ac-H3","Proximal"}'  > $TEMP/H3K9ac-H3_Proximal.density
sed '1d' $TEMP/H3K9ac-H3_Distal.tab    | awk 'BEGIN {OFS="\t"}{z = (log(($2+1)/($3+1))/log(2)); print $1,z,"H3K9ac-H3","Distal"}'    > $TEMP/H3K9ac-H3_Distal.density
sed '1d' $TEMP/H3K27ac-H3_Proximal.tab | awk 'BEGIN {OFS="\t"}{z = (log(($2+1)/($3+1))/log(2)); print $1,z,"H3K27ac-H3","Proximal"}' > $TEMP/H3K27ac-H3_Proximal.density
sed '1d' $TEMP/H3K27ac-H3_Distal.tab   | awk 'BEGIN {OFS="\t"}{z = (log(($2+1)/($3+1))/log(2)); print $1,z,"H3K27ac-H3","Distal"}'   > $TEMP/H3K27ac-H3_Distal.density

# Compile density info
cat $TEMP/*.density | gzip -c > S7/h/DensityInfo.tab.gz

# Generate violin plot
python $VIOLIN -i <(gzip -dc S7/h/DensityInfo.tab.gz | cut -f2,3,4) -o S7/h/DensityInfo.svg \
	--width 3 --height 4 --preset2 \
	--title "Density at +1 nucleosome" \
	--xlabel "modification" --ylabel "Density (log2)"

# Clean-up
rm -r $TEMP

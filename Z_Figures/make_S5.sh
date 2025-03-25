#!/bin/bash

# Organize data from X_Bulk_Processing into Z_Figures for S5

### CHANGE ME
WRK=/storage/group/bfp2/default/hxc585_HainingChen/2025_Chen_TF-Nuc/Z_Figures
#WRK=/storage/home/owl5022/scratch/2024-Krebs_Science/Z_Figures
###

set -exo
module load anaconda3
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

LIBRARY=$WRK/../X_Bulk_Processing/Library
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar

[ -d S5 ] || mkdir S5

# ===============================================================================================================================

[ -d S5/a ] || mkdir S5/a
##heatmap
cp $LIBRARY/10phase/FOXA_LABEL-uHepG2_SORT-ClosestDyad/BNase-seq_50U-10min_merge_hg38_FOXA_LABEL-uHepG2_SORT-ClosestDyad_nearNuc_read1_shuffle_250bp_phase_sort_merge.svg S5/a/
cp $LIBRARY/10phase/FOXA_LABEL-uHepG2_SORT-ClosestDyad/BNase-seq_50U-10min_merge_hg38_FOXA_LABEL-uHepG2_SORT-ClosestDyad_nearNuc_read1_original_250bp_phase_sort_merge.svg S5/a/
# Composites
cp $LIBRARY/10phase/FOXA_LABEL-uHepG2_SORT-ClosestDyad/10xplot/BNase-seq_50U-10min_merge_hg38_FOXA_LABEL-uHepG2_SORT-ClosestDyad_nearNuc_read1_original_phase_2.out S5/a/
cp $LIBRARY/10phase/FOXA_LABEL-uHepG2_SORT-ClosestDyad/10xplot/BNase-seq_50U-10min_merge_hg38_FOXA_LABEL-uHepG2_SORT-ClosestDyad_nearNuc_read1_original_phase_0.out S5/a/
cp $LIBRARY/10phase/FOXA_LABEL-uHepG2_SORT-ClosestDyad/10xplot/BNase-seq_50U-10min_merge_hg38_FOXA_LABEL-uHepG2_SORT-ClosestDyad_nearNuc_read1_shuffle_phase_2.out S5/a/
cp $LIBRARY/10phase/FOXA_LABEL-uHepG2_SORT-ClosestDyad/10xplot/BNase-seq_50U-10min_merge_hg38_FOXA_LABEL-uHepG2_SORT-ClosestDyad_nearNuc_read1_shuffle_phase_0.out S5/a
cp $LIBRARY/10phase/FOXA_phase_aligned/*.out S5/a/


for file in S5/a/YY_*.out S5/a/RR_*.out S5/a/SS_*.out S5/a/WW_*.out ; do
    filename=$(basename "$file" ".out")
    cut -f 1-231 "$file" > S5/a/${filename}_sense_left.out
    echo "$(seq -19 19 | tr '\n' ' ')" > S5/a/${filename}_sense_mid.out
    for i in $(seq 1 $(($(wc -l < "$file") - 1))); do
        echo -e "$(printf '\t%.0s-' {1..39})" >> S5/a/${filename}_sense_mid_placeholder.out
    done
    cat S5/a/${filename}_sense_mid.out S5/a/${filename}_sense_mid_placeholder.out > S5/a/${filename}_sense_mid_combined.out
    rm S5/a/${filename}_sense_mid.out S5/a/${filename}_sense_mid_placeholder.out
    cut -f 271-501 "$file" > S5/a/${filename}_sense_right.out
    paste S5/a/${filename}_sense_left.out S5/a/${filename}_sense_mid_combined.out S5/a/${filename}_sense_right.out > S5/a/${filename}_sense_modified.out
    rm S5/a/${filename}_sense_left.out S5/a/${filename}_sense_mid_combined.out S5/a/${filename}_sense_right.out
done


# ===============================================================================================================================
[ -d S5/b ] || mkdir S5/b
# Composites
BED=FOXA_LABEL-K562_SORT-ClosestDyad_GROUP-NFR_1000bp
cp $LIBRARY/$BED/Composites/K562_FOXA2_BX_rep1_hg38_${BED}_5read1_NCIS.out S5/b/
cp $LIBRARY/$BED/Composites/K562_FOXA2_BX_rep1_hg38_${BED}_5read2_NCIS.out S5/b/
BED=FOXA_LABEL-K562_SORT-ClosestDyad_GROUP-NoOverlap_1000bp
cp $LIBRARY/$BED/Composites/K562_FOXA2_BX_rep1_hg38_${BED}_5read1_NCIS.out S5/b/
cp $LIBRARY/$BED/Composites/K562_FOXA2_BX_rep1_hg38_${BED}_5read2_NCIS.out S5/b/



# ===============================================================================================================================
[ -d S5/c ] || mkdir S5/c

# Heatmaps
BED=FOXA_LABEL-HepG2_SORT-NucleosomeEngagement_500bp
cp $LIBRARY/$BED/SVG/HepG2_FOXA1_BX_rep1_hg38_${BED}_5read1_NCIS_merge_label.svg S5/c/
cp $LIBRARY/$BED/SVG/HepG2_FOXA1_BX_rep1_hg38_${BED}_5read2_NCIS_merge_label.svg S5/c/

# Composites
BED=FOXA_LABEL-HepG2_SORT-NucleosomeEngagement_GROUP-Engaged_1000bp
cp $LIBRARY/$BED/Composites/HepG2_FOXA1_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS.out S5/c/
cp $LIBRARY/$BED/Composites/HepG2_FOXA1_BX_rep1_hg38_${BED}_5read2-MIN100_NCIS.out S5/c/
BED=FOXA_LABEL-HepG2_SORT-NucleosomeEngagement_GROUP-LessEngaged_1000bp
cp $LIBRARY/$BED/Composites/HepG2_FOXA1_BX_rep1_hg38_${BED}_5read1-MIN100_NCIS.out S5/c/
cp $LIBRARY/$BED/Composites/HepG2_FOXA1_BX_rep1_hg38_${BED}_5read2-MIN100_NCIS.out S5/c/

[ -d S5/c ] || mkdir S5/c

# Heatmaps
BED=FOXA_HepG2_SORT-ClosestHNF4A_1000bp
cp $LIBRARY/$BED/SVG/HepG2_FOXA1_BX_rep1_hg38_${BED}_5read1_NCIS_merge_label.svg S5/d/
cp $LIBRARY/$BED/SVG/HepG2_FOXA1_BX_rep1_hg38_${BED}_5read2_NCIS_merge_label.svg S5/d/
cp $LIBRARY/$BED/SVG/HepG2_HNF4A_BX_rep1_hg38_${BED}_5read1_NCIS_merge_label.svg S5/d/

## motif alignment

[ -d S5/e ] || mkdir S5/e
# Composites
BED=FOXA_bound-HepG2_SORT-ClosestHNF4A-bound_Group-Distal_1000bp
cp $LIBRARY/$BED/Composites/HepG2_FOXA1_BX_rep1_hg38_${BED}_5read1_NCIS.out S5/e/
cp $LIBRARY/$BED/Composites/HepG2_HNF4A_BX_rep1_hg38_${BED}_5read1_NCIS.out S5/e/

BED=FOXA_bound-HepG2_SORT-ClosestHNF4A-bound_Group-Proximal_1000bp
cp $LIBRARY/$BED/Composites/HepG2_FOXA1_BX_rep1_hg38_${BED}_5read1_NCIS.out S5/e/
cp $LIBRARY/$BED/Composites/HepG2_HNF4A_BX_rep1_hg38_${BED}_5read1_NCIS.out S5/e/

BED=HNF4A_bound-HepG2_SORT-ClosestFOXA-bound_Group-Distal_1000bp
cp $LIBRARY/$BED/Composites/HepG2_FOXA1_BX_rep1_hg38_${BED}_5read1_NCIS.out S5/e/
cp $LIBRARY/$BED/Composites/HepG2_HNF4A_BX_rep1_hg38_${BED}_5read1_NCIS.out S5/e/

BED=HNF4A_bound-HepG2_SORT-ClosestFOXA-bound_Group-Proximal_1000bp
cp $LIBRARY/$BED/Composites/HepG2_FOXA1_BX_rep1_hg38_${BED}_5read1_NCIS.out S5/e/
cp $LIBRARY/$BED/Composites/HepG2_HNF4A_BX_rep1_hg38_${BED}_5read1_NCIS.out S5/e/


[ -d S5/f ] || mkdir S5/f
BED=FOXA_bound-HepG2_SORT-ClosestHNF4A-bound_Group-Proximal_1000bp
cp $LIBRARY/$BED/Composites/HepG2_FOXA1_BX_rep1_hg38_${BED}_5read2-MIN100_NCIS.out S5/f/
cp $LIBRARY/$BED/Composites/HepG2_HNF4A_BX_rep1_hg38_${BED}_5read2-MIN100_NCIS.out S5/f/

##

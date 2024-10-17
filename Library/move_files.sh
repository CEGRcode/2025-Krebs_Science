#!/bin/bash

# Organize data from 0X_Bulk_Processing into Z_Figures

### CHANGE ME

WRK=/storage/group/bfp2/default/hxc585_HainingChen/Library
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar

[ -d SF3/c ] || mkdir SF3/c
cp $LIBRARY/ZKSCAN1_Occupancy_flip_shift2_1000bp/Composites/K562_ZKSCAN1_BX_rep1_hg38_ZKSCAN1_Occupancy_flip_shift2_MIN100_read2_Normalized.out SF3/c
cp $LIBRARY/ZKSCAN1_Occupancy_flip_shift2_1000bp/Composites/K562_ZKSCAN1_BX_rep1_hg38_ZKSCAN1_Occupancy_flip_shift2_MIN100_read1_Normalized.out SF3/c
cp $LIBRARY/ZKSCAN1_Occupancy_flip_shift2_1000bp/Composites/K562_-_BI_rep1_hg38_ZKSCAN1_Occupancy_flip_shift2_read2.out SF3/c
cp $LIBRARY/ZKSCAN1_Occupancy_flip_shift2_1000bp/Composites/K562_-_BI_rep1_hg38_ZKSCAN1_Occupancy_flip_shift2_read1.out SF3/c

[ -d F4/a ] || mkdir F4/a
cp $LIBRARY/WebLogos/CTCF_M1_logo.eps F4/a/
cp $LIBRARY/CTCF_Occupancy_1000bp/Composites/K562_CTCF_BX_rep1_hg38_CTCF_Occupancy_1000bp_read1_Normalized.out F4/a/
cp $LIBRARY/CTCF_Occupancy_1000bp/Composites/K562_CTCF_BX_rep1_hg38_CTCF_Occupancy_1000bp_read2_Normalized.out F4/a/
cp $LIBRARY/CTCF_Occupancy_1000bp/Composites/K562_-_BI_rep1_hg38_CTCF_Occupancy_1000bp_midpoint.out F4/a/

[ -d F4/b ] || mkdir F4/b
cp $LIBRARY/CTCF_Occupancy_1000bp/FourColor/CTCF_Occupancy_1000bp_31bp.svg F4/b/
cp $LIBRARY/CTCF_Occupancy_1000bp/SVG/Input_CTCF_Occupancy_1000bp_midpoint_combined_treeview_label.svg F4/b/
cp $LIBRARY/CTCF_Occupancy_1000bp/SVG/K562_CTCF_BX_rep1_hg38_CTCF_Occupancy_1000bp_read1-MIN100_Normalized_merge_label.svg F4/b/
cp $LIBRARY/CTCF_Occupancy_1000bp/SVG/K562_CTCF_BX_rep1_hg38_CTCF_Occupancy_1000bp_read2-MIN100_Normalized_merge_label.svg F4/b/
cp $LIBRARY/CTCF_Occupancy_1000bp/SVG/K562_CTCF_BX_rep1_hg38_CTCF_Occupancy_1000bp_read2_Normalized_merge_label F4/b/
cp $LIBRARY/CTCF_Occupancy_1000bp/SVG/K562_CTCF_BX_rep1_hg38_CTCF_Occupancy_1000bp_read1_Normalized_merge_label.svg F4/b/

[ -d F3/c ] || mkdir F4/c
cp $LIBRARY/WebLogos/CTCF_logo.eps F4/a/
cp $LIBRARY/CTCF_Occupancy_1000bp/Composites/K562_CTCF_BX_rep1_hg38_CTCF_Occupancy_1000bp_read1-MIN100_Normalized.out F4/c/
cp $LIBRARY/CTCF_Occupancy_1000bp/Composites/K562_CTCF_BX_rep1_hg38_CTCF_Occupancy_1000bp_read2-MIN100_Normalized.out F4/c/
cp $LIBRARY/CTCF_Occupancy_1000bp/Composites/K562_-_BI_rep1_hg38_CTCF_Occupancy_1000bp_midpoint.out F4/c/

[ -d F2/d ] || mkdir F4/d
cp $LIBRARY/CTCF_Occupancy_1000bp/Composites/K562_SMC3_BX_rep1_hg38_CTCF_Occupancy_1000bp_read1_Normalized.out F4/d/
cp $LIBRARY/CTCF_Occupancy_1000bp/Composites/K562_RAD21_BX_rep1_hg38_CTCF_Occupancy_1000bp_read2_Normalized.out F4/d/
cp $LIBRARY/CTCF_Occupancy_1000bp/Composites/K562_CTCF_BX_rep1_hg38_CTCF_Occupancy_1000bp_read1_Normalized.out F4/d/


[ -d F4/e ] || mkdir F4/e
cp $LIBRARY/WebLogos/CTCF_logo.eps F2/a/
cp $LIBRARY/CTCF_Occupancy_1000bp/Composites/K562_SMC3_BX_rep1_hg38_CTCF_Occupancy_1000bp_read1-MIN100_Normalized.out F4/e/
cp $LIBRARY/CTCF_Occupancy_1000bp/Composites/K562_SMC3_BX_rep1_hg38_CTCF_Occupancy_1000bp_read2-MIN100_Normalized.out F4/e/
cp $LIBRARY/CTCF_Occupancy_1000bp/Composites/K562_RAD21_BX_rep1_hg38_CTCF_Occupancy_1000bp_read1-MIN100_Normalized.out F4/e/
cp $LIBRARY/CTCF_Occupancy_1000bp/Composites/K562_RAD21_BX_rep1_hg38_CTCF_Occupancy_1000bp_read2-MIN100_Normalized.out F4/e/
cp $LIBRARY/CTCF_Occupancy_1000bp/Composites/K562_CTCF_BX_rep1_hg38_CTCF_Occupancy_1000bp_read1-MIN100_Normalized.out F4/e/
cp $LIBRARY/CTCF_Occupancy_1000bp/Composites/K562_CTCF_BX_rep1_hg38_CTCF_Occupancy_1000bp_read2-MIN100_Normalized.out F4/e/
cp $LIBRARY/CTCF_Occupancy_1000bp/Composites/K562_-_BI_rep1_hg38_CTCF_Occupancy_1000bp_midpoint.out F4/e/
cp $LIBRARY/CTCF_Occupancy_1000bp/Composites/K562_-_BI_rep1_hg38_CTCF_Occupancy_1000bp_read1.out F4/e/
cp $LIBRARY/CTCF_Occupancy_1000bp/Composites/K562_-_BI_rep1_hg38_CTCF_Occupancy_1000bp_read2.out F4/e/


# ===============================================================================================================================

[ -d F5 ] || mkdir F5

[ -d F5/a ] || mkdir F5/a
cp $LIBRARY/WebLogos/FOXA2_M1_logo.eps F5/a/
cp $LIBRARY/FOXA_all_1000bp/FourColor/FOXA_all_1000bp_31bp.svg F5/a/
cp $LIBRARY/FOXA_all_1000bp/SVG/K562_IgG_BX_merge_hg38_FOXA_all_1000bp_read1_Normalized_merge_lable.svg F5/a/
cp $LIBRARY/FOXA_all_1000bp/SVG/K562_FOXA2_BX_rep1_hg38_FOXA_all_1000bp_read1_Normalized_merge_lable.svg F5/a/
cp $LIBRARY/FOXA_all_1000bp/SVG/K562_FOXA1_BX_rep1_hg38_FOXA_all_1000bp_read1_Normalized_merge_lable.svg F5/a/
cp $LIBRARY/FOXA_all_1000bp/SVG/K562_-_BI_rep1_hg38_FOXA_all_1000bp_midpoint_combined_lable.svg F5/a/
cp $LIBRARY/FOXA_all_1000bp/SVG/HepG2_IgG_BX_merge_hg38_FOXA_all_1000bp_read1_Normalized_merge_lable.svg F5/a/
cp $LIBRARY/FOXA_all_1000bp/SVG/HepG2_FOXA1_BX_rep1_hg38_FOXA_all_1000bp_read1_Normalized_merge_lable.svg F5/a/

[ -d F5/b ] || mkdir F5/b
cp $LIBRARY/FOXA_uniq_HepG2_NucSort_1000bp/Composites/K562_-_BI_rep1_hg38_FOXA_uniq_HepG2_NucSort_1000bp_read2.out F5/b/
cp $LIBRARY/FOXA_uniq_HepG2_NucSort_1000bp/Composites/K562_-_BI_rep1_hg38_FOXA_uniq_HepG2_NucSort_1000bp_read1.out F5/b/
cp $LIBRARY/FOXA_K562_NucSort_1000bp/Composites/K562_-_BI_rep1_hg38_FOXA_K562_NucSort_1000bp_read2.out F5/b/
cp $LIBRARY/FOXA_K562_NucSort_1000bp/Composites/K562_-_BI_rep1_hg38_FOXA_K562_NucSort_1000bp_read1.out F5/b/

[ -d F5/c ] || mkdir F5/c
cp $LIBRARY/FOXA_K562_NucSort-OVERLAP_1000bp/Composites/K562_FOXA1_BX_rep1_hg38_FOXA_K562_NucSort-OVERLAP_1000bp_read2_Normalized.out F5/c/
cp $LIBRARY/FOXA_K562_NucSort-OVERLAP_1000bp/Composites/K562_FOXA1_BX_rep1_hg38_FOXA_K562_NucSort-OVERLAP_1000bp_read1_Normalized.out F5/c/
cp $LIBRARY/FOXA_K562_NucSort-NFR_1000bp/Composites/K562_FOXA1_BX_rep1_hg38_FOXA_K562_NucSort-NFR_1000bp_read2_Normalized.out F5/c/
cp $LIBRARY/FOXA_K562_NucSort-NFR_1000bp/Composites/K562_FOXA1_BX_rep1_hg38_FOXA_K562_NucSort-NFR_1000bp_read1_Normalized.out F5/c/

[ -d SF4/a ] || mkdir SF4/a
cp $LIBRARY/FOXA_K562_NucSort-OVERLAP_1000bp/Composites/K562_FOXA2_BX_rep1_hg38_FOXA_K562_NucSort-OVERLAP_1000bp_read2_Normalized.out SF4/a
cp $LIBRARY/FOXA_K562_NucSort-OVERLAP_1000bp/Composites/K562_FOXA2_BX_rep1_hg38_FOXA_K562_NucSort-OVERLAP_1000bp_read1_Normalized.out SF4/a
cp $LIBRARY/FOXA_K562_NucSort-NFR_1000bp/Composites/K562_FOXA2_BX_rep1_hg38_FOXA_K562_NucSort-NFR_1000bp_read2_Normalized.out SF4/a
cp $LIBRARY/FOXA_K562_NucSort-NFR_1000bp/Composites/K562_FOXA2_BX_rep1_hg38_FOXA_K562_NucSort-NFR_1000bp_read1_Normalized.out SF4/a

[ -d F5/d ] || mkdir F5/d
cp $LIBRARY/FOXA_K562_NucengageSort_500bp/SVG/K562_FOXA1_BX_rep1_hg38_FOXA_K562_NucengageSort_500bp_read2_Normalized_merge_lable.svg F5/d/
cp $LIBRARY/FOXA_K562_NucengageSort_500bp/SVG/K562_FOXA1_BX_rep1_hg38_FOXA_K562_NucengageSort_500bp_read1_Normalized_merge_lable.svg F5/d/
cp $LIBRARY/FOXA_K562_Nucengage_1000bp/Composites/K562_FOXA1_BX_rep1_hg38_FOXA_K562_Nucengage_1000bp_read2-MIN100_Normalized.out F5/d/
cp $LIBRARY/FOXA_K562_Nucengage_1000bp/Composites/K562_FOXA1_BX_rep1_hg38_FOXA_K562_Nucengage_1000bp_read1-MIN100_Normalized.out F5/d/
cp $LIBRARY/FOXA_K562_Nucengage_1000bp/Composites/K562_-_BI_rep1_hg38_FOXA_K562_Nucengage_1000bp_read2.out F5/d/
cp $LIBRARY/FOXA_K562_Nucengage_1000bp/Composites/K562_-_BI_rep1_hg38_FOXA_K562_Nucengage_1000bp_read1.out F5/d/
cp $LIBRARY/FOXA_K562_NoNuc_1000bp/Composites/K562_FOXA1_BX_rep1_hg38_FOXA_K562_NoNuc_1000bp_read2-MIN100_Normalized.out F5/d/
cp $LIBRARY/FOXA_K562_NoNuc_1000bp/Composites/K562_FOXA1_BX_rep1_hg38_FOXA_K562_NoNuc_1000bp_read1-MIN100_Normalized.out F5/d/
cp $LIBRARY/FOXA_K562_NoNuc_1000bp/Composites/K562_-_BI_rep1_hg38_FOXA_K562_NoNuc_1000bp_read1.out F5/d/
cp $LIBRARY/FOXA_K562_NoNuc_1000bp/Composites/K562_-_BI_rep1_hg38_FOXA_K562_NoNuc_1000bp_read1.out F5/d/


[ -d SF4/b ] || mkdir SF4/b
cp $LIBRARY/FOXA_uniq_HepG2_NucengageSort_500bp/SVG/HepG2_FOXA1_BX_rep1_hg38_FOXA_uniq_HepG2_NucengageSort_500bp_read2_Normalized_merge_lable.svg SF4/b
cp $LIBRARY/FOXA_uniq_HepG2_NucengageSort_500bp/SVG/epG2_FOXA1_BX_rep1_hg38_FOXA_uniq_HepG2_NucengageSort_500bp_read1_Normalized_merge_lable.svg SF4/b
cp $LIBRARY/FOXA_HepG2_Nucengage_1000bp/Composites/HepG2_FOXA1_BX_rep1_hg38_FOXA_HepG2_Nucengage_1000bp_read2-MIN100_Normalized.out SF4/b
cp $LIBRARY/FOXA_HepG2_Nucengage_1000bp/Composites/HepG2_FOXA1_BX_rep1_hg38_FOXA_HepG2_Nucengage_1000bp_read1-MIN100_Normalized.out SF4/b
cp $LIBRARY/FOXA_HepG2_NoNuc_1000bp/Composites/HepG2_FOXA1_BX_rep1_hg38_FOXA_HepG2_NoNuc_1000bp_read2-MIN100_Normalized.out SF4/b
cp $LIBRARY/FOXA_HepG2_NoNuc_1000bp/Composites/HepG2_FOXA1_BX_rep1_hg38_FOXA_HepG2_NoNuc_1000bp_read1-MIN100_Normalized.out SF4/b

# ===============================================================================================================================

[ -d F6/a ] || mkdir F6/a
cp $LIBRARY/NFIA_downNuc_500bp/Composites/dbSnp153_snv_NFIA_downNuc_500bp.out F6/a
cp $LIBRARY/NFIA_downNuc_500bp/Composites/hg38.phyloP30way.NFIA_downNuc_500bp.out F6/a
cp $LIBRARY/WebLogos/NFIA_M1_logo.eps F6/a
cp $LIBRARY/NFIA_downNuc_4_1000bp/Composites/K562_-_BI_rep1_hg38_NFIA_downNuc_4_1000bp_read2.out F6/a
cp $LIBRARY/NFIA_downNuc_4_1000bp/Composites/K562_-_BI_rep1_hg38_NFIA_downNuc_4_1000bp_read1.out F6/a

[ -d F6/b ] || mkdir F6/b
cp $LIBRARY/NFIA_downNuc_500bp/FourColor/NFIA_downNuc_500bp_31bp.svg F6/b/
cp $LIBRARY/NFIA_downNuc_500bp/SVG/K562_NFIA_BX_rep1_hg38_NFIA_downNuc_500bp_read1_Normalized_merge_label.svg F6/b/
cp $LIBRARY/NFIA_downNuc_500bp/SVG/K562_IgG_BX_merge_hg38_NFIA_downNuc_500bp_read1_Normalized_merge_label.svg F6/b/
cp $LIBRARY/NFIA_downNuc_down125_250bp/SVG/K562_NFIA_BX_rep1_hg38_NFIA_downNuc_down125_250bp_read1-MIN100_Normalized_merge_label.svg F6/b/


[ -d F6/c ] || mkdir F6/c
cp $LIBRARY/NFIA_engageNucdown_NucSort_1000bp/SVG/K562_-_BI_rep1_hg38_NFIA_engageNucdown_NucSort_1000bp_midpoint_combined_treeview_label.svg  F6/c
cp $LIBRARY/NFIA_engageNucdown_NucSort_1000bp/SVG/K562_NFIA_BX_rep1_hg38_NFIA_engageNucdown_NucSort_1000bp_read1_merge_treeview_label.svg  F6/c
cp $LIBRARY/NFIA_NucSort-DOWNSTREAM_1000/Composites/K562_NFIA_BX_rep1_hg38_NFIA_NucSort-DOWNSTREAM_1000bp_read1-MIN100_Normalized.out F6/c
cp $LIBRARY/NFIA_NucSort-DOWNSTREAM_1000/Composites/K562_NFIA_BX_rep1_hg38_NFIA_NucSort-DOWNSTREAM_1000bp_read2-MIN100_Normalized.out F6/c
cp $LIBRARY/NFIA_NucSort-OVERLAP_1000/Composites/K562_NFIA_BX_rep1_hg38_NFIA_NucSort-OVERLAP_1000bp_read1-MIN100_Normalized.out F6/c
cp $LIBRARY/NFIA_NucSort-OVERLAP_1000/Composites/K562_NFIA_BX_rep1_hg38_NFIA_NucSort-OVERLAP_1000bp_read2-MIN100_Normalized.out F6/c
cp $LIBRARY/NFIA_NucSort-UPSTREAM_1000/Composites/K562_NFIA_BX_rep1_hg38_NFIA_NucSort-UPSTREAM_1000bp_read1-MIN100_Normalized.out F6/c
cp $LIBRARY/NFIA_NucSort-UPSTREAM_1000/Composites/K562_NFIA_BX_rep1_hg38_NFIA_NucSort-UPSTREAM_1000bp_read2-MIN100_Normalized.out F6/c


# ===============================================================================================================================

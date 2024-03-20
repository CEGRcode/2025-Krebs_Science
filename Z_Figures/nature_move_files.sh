#!/bin/bash

# Organize data from 0X_Bulk_Processing into Z_Figures

### CHANGE ME
WRK=/path/to/2024-Chen_Nature/Z_Figures
WRK=/storage/home/owl5022/scratch/2024-Chen_Nature/Z_Figures
###

LIBRARY=$WRK/../0X_Bulk_Processing/Library
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.14.jar

[ -d F2 ] || mkdir F2
REFPT=$LIBRARY/NFIA_Occupancy_1000bp/FourColor/NFIA_Occupancy_1000bp.svg

[ -d F2/b ] || mkdir F2/a
cp $LIBRARY/WebLogos/NFIA_logo.eps F2/a/

[ -d F2/b ] || mkdir F2/b
cp $LIBRARY/NFIA_Occupancy_1000bp/FourColor/NFIA_Occupancy_1000bp_31bp.svg F2/b/
cp $LIBRARY/NFIA_Occupancy_1000bp/SVG/K562_NFIA_BX_rep1_hg19_NFIA_Occupancy_1000bp_read1_Normalized_merge_label.svg F2/b/
cp $LIBRARY/NFIA_Occupancy_1000bp/SVG/K562_IgG_BX_merge_hg19_NFIA_Occupancy_1000bp_read1_Normalized_merge_label.svg F2/b/

[ -d F2/c ] || mkdir F2/c
cp $LIBRARY/NFIA_Occupancy_1000bp/Composites/primates.phyloP46way_NFIA_Occupancy_1000bp.out F2/c/1_primates.phyloP46way_NFIA_Occupancy_1000bp.out
cp $LIBRARY/NFIA_Occupancy_1000bp/Composites/dbSnp153_snv_NFIA_Occupancy_1000bp.out F2/c/2_dbSnp153_snv_NFIA_Occupancy_1000bp.out

[ -d F2/d ] || mkdir F2/d
cp $LIBRARY/NFIA_NucSort_1000bp/Composites/K562_NFIA_BX_rep1_hg19_NFIA_NucSort_1000bp_read1-MIN100_Normalized.out F2/d/1_K562_NFIA_BX_rep1_hg19_NFIA_NucSort_1000bp_read1-MIN100_Normalized.out
cp $LIBRARY/NFIA_NucSort_1000bp/Composites/K562_NFIA_BX_rep1_hg19_NFIA_NucSort_1000bp_read2-MIN100_Normalized.out F2/d/2_K562_NFIA_BX_rep1_hg19_NFIA_NucSort_1000bp_read2-MIN100_Normalized.out
cp $LIBRARY/NFIA_Occupancy_1000bp/Composites/K562_IgG_BX_merge_hg19_NFIA_Occupancy_1000bp_read1_Normalized.out F2/d/3_K562_IgG_BX_merge_hg19_NFIA_Occupancy_1000bp_read1_Normalized.out

# Separate script for F2/e

[ -d F2/f ] || mkdir F2/f
cp $LIBRARY/NFIA_NucSort_1000bp/SVG/K562_-_BI_rep1_hg19_NFIA_NucSort_1000bp_midpoint_combined_treeview_label.svg F2/f/
cp $LIBRARY/NFIA_NucSort_1000bp/SVG/K562_H3K27ac_BX_rep1_hg19_NFIA_NucSort_1000bp_midpoint_combined_treeview_label.svg F2/f/
cp $LIBRARY/NFIA_NucSort_1000bp/SVG/K562_NFIA_BX_rep1_hg19_NFIA_NucSort_1000bp_read1_Normalized_merge_label.svg F2/f/
# midpoint or endo?

[ -d F2/g ] || mkdir F2/g
cp $LIBRARY/NFIA_NucSort_1000bp/Composites/K562_-_BI_rep1_hg19_NFIA_NucSort_1000bp_read2.out F2/g/1_K562_-_BI_rep1_hg19_NFIA_NucSort_1000bp_read2.out

[ -d F2/h ] || mkdir F2/h
cp $LIBRARY/NFIA_NucSort-DOWNSTREAM_1000bp/Composites/K562_NFIA_BX_rep1_hg19_NFIA_NucSort-DOWNSTREAM_1000bp_read1-MIN100_Normalized.out F2/h/1_K562_NFIA_BX_rep1_hg19_NFIA_NucSort-DOWNSTREAM_1000bp_read1-MIN100_Normalized.out
cp $LIBRARY/NFIA_NucSort-DOWNSTREAM_1000bp/Composites/K562_NFIA_BX_rep1_hg19_NFIA_NucSort-DOWNSTREAM_1000bp_read2-MIN100_Normalized.out F2/h/2_K562_NFIA_BX_rep1_hg19_NFIA_NucSort-DOWNSTREAM_1000bp_read2-MIN100_Normalized.out
cp $LIBRARY/NFIA_Occupancy_1000bp/Composites/K562_IgG_BX_merge_hg19_NFIA_Occupancy_1000bp_read1_Normalized.out F2/h/3_K562_IgG_BX_merge_hg19_NFIA_Occupancy_1000bp_read1_Normalized.out
cp $LIBRARY/NFIA_NucSort-OVERLAP_1000bp/Composites/K562_NFIA_BX_rep1_hg19_NFIA_NucSort-OVERLAP_1000bp_read1-MIN100_Normalized.out F2/h/1_K562_NFIA_BX_rep1_hg19_NFIA_NucSort-OVERLAP_1000bp_read1-MIN100_Normalized.out
cp $LIBRARY/NFIA_NucSort-OVERLAP_1000bp/Composites/K562_NFIA_BX_rep1_hg19_NFIA_NucSort-OVERLAP_1000bp_read2-MIN100_Normalized.out F2/h/2_K562_NFIA_BX_rep1_hg19_NFIA_NucSort-OVERLAP_1000bp_read2-MIN100_Normalized.out
cp $LIBRARY/NFIA_NucSort-UPSTREAM_1000bp/Composites/K562_NFIA_BX_rep1_hg19_NFIA_NucSort-UPSTREAM_1000bp_read1-MIN100_Normalized.out F2/h/1_K562_NFIA_BX_rep1_hg19_NFIA_NucSort-UPSTREAM_1000bp_read1-MIN100_Normalized.out
cp $LIBRARY/NFIA_NucSort-UPSTREAM_1000bp/Composites/K562_NFIA_BX_rep1_hg19_NFIA_NucSort-UPSTREAM_1000bp_read2-MIN100_Normalized.out F2/h/2_K562_NFIA_BX_rep1_hg19_NFIA_NucSort-UPSTREAM_1000bp_read2-MIN100_Normalized.out

[ -d F2/i ] || mkdir F2/i
# Aggregate occupancy scores across both replicates (Down, Overlap, Up)


# ===============================================================================================================================

[ -d F3 ] || mkdir F3

[ -d F3/a ] || mkdir F3/a
cp $LIBRARY/WebLogos/CTCF_logo.eps F3/a/

[ -d F3/b ] || mkdir F3/b
cp $LIBRARY/CTCF_Occupancy_1000bp/FourColor/CTCF_Occupancy_1000bp_31bp.svg F3/b/
cp $LIBRARY/CTCF_Occupancy_1000bp/SVG/K562_CTCF_BX_rep1_hg19_CTCF_Occupancy_1000bp_read1_Normalized_merge_label.svg F3/b/
cp $LIBRARY/CTCF_Occupancy_1000bp/SVG/K562_IgG_BX_merge_hg19_CTCF_Occupancy_1000bp_read1_Normalized_merge_label.svg F3/b/
cp $LIBRARY/CTCF_Occupancy_1000bp/SVG/K562_-_BI_rep1_hg19_CTCF_Occupancy_1000bp_midpoint_combined_treeview_label.svg F3/b/
# Aggregate occupancy scores across both replicates (Top, Mid, Bottom)

[ -d F3/c ] || mkdir F3/c
cp $LIBRARY/CTCF_Occupancy_1000bp/Composites/K562_CTCF_BX_rep1_hg19_CTCF_Occupancy_1000bp_read1_Normalized.out F3/c/1_K562_CTCF_BX_rep1_hg19_CTCF_Occupancy_1000bp_read1_Normalized.out
cp $LIBRARY/CTCF_Occupancy_1000bp/Composites/K562_CTCF_BX_rep1_hg19_CTCF_Occupancy_1000bp_read2_Normalized.out F3/c/2_K562_CTCF_BX_rep1_hg19_CTCF_Occupancy_1000bp_read2_Normalized.out
cp $LIBRARY/CTCF_Occupancy_1000bp/Composites/K562_-_BI_rep1_hg19_CTCF_Occupancy_1000bp_midpoint.out F3/c/3_K562_-_BI_rep1_hg19_CTCF_Occupancy_1000bp_midpoint.out
# cp $LIBRARY/CTCF_Occupancy_1000bp/Composites/K562_IgG_BX_merge_hg19_CTCF_Occupancy_1000bp_read1_Normalized.out F3/c/3_K562_IgG_BX_merge_hg19_CTCF_Occupancy_1000bp_read1_Normalized.out

[ -d F3/d ] || mkdir F3/d
cp $LIBRARY/CTCF_Occupancy_1000bp/Composites/K562_-_BI_rep1_hg19_CTCF_Occupancy_1000bp_read2.out F3/d/1_K562_-_BI_rep1_hg19_CTCF_Occupancy_1000bp_read2.out

[ -d F3/e ] || mkdir F3/e
cp $LIBRARY/CTCF_Occupancy_1000bp/Composites/K562_CTCF_BX_rep1_hg19_CTCF_Occupancy_1000bp_read1-MIN100_Normalized.out F3/e/1_K562_CTCF_BX_rep1_hg19_CTCF_Occupancy_1000bp_read1-MIN100_Normalized.out
cp $LIBRARY/CTCF_Occupancy_1000bp/Composites/K562_CTCF_BX_rep1_hg19_CTCF_Occupancy_1000bp_read2-MIN100_Normalized.out F3/e/2_K562_CTCF_BX_rep1_hg19_CTCF_Occupancy_1000bp_read2-MIN100_Normalized.out
cp $LIBRARY/CTCF_Occupancy_1000bp/Composites/K562_IgG_BX_merge_hg19_CTCF_Occupancy_1000bp_read1_Normalized.out F3/e/3_K562_IgG_BX_merge_hg19_CTCF_Occupancy_1000bp_read1_Normalized.out
cp $LIBRARY/CTCF_Occupancy_1000bp/Composites/K562_-_BI_rep1_hg19_CTCF_Occupancy_1000bp_midpoint.out F3/e/4_K562_-_BI_rep1_hg19_CTCF_Occupancy_1000bp_midpoint.out

[ -d F3/f ] || mkdir F3/f
cp $LIBRARY/CTCF_Occupancy_1000bp/Composites/K562_CTCF_BX_rep1_hg19_CTCF_Occupancy_1000bp_read1_Normalized.out F3/f/CTCF_K562_CTCF_BX_rep1_hg19_CTCF_Occupancy_1000bp_read1_Normalized.out
cp $LIBRARY/CTCF_Occupancy_1000bp/Composites/K562_RAD21_BX_rep1_hg19_CTCF_Occupancy_1000bp_read1_Normalized.out F3/f/RAD21_K562_RAD21_BX_rep1_hg19_CTCF_Occupancy_1000bp_read1_Normalized.out
cp $LIBRARY/CTCF_Occupancy_1000bp/Composites/K562_SMC3_BX_rep1_hg19_CTCF_Occupancy_1000bp_read1_Normalized.out F3/f/SMC3_K562_SMC3_BX_rep1_hg19_CTCF_Occupancy_1000bp_read1_Normalized.out

[ -d F3/g ] || mkdir F3/g
cp $LIBRARY/CTCF_Occupancy_1000bp/Composites/K562_RAD21_BX_rep1_hg19_CTCF_Occupancy_1000bp_read1-MIN100_Normalized.out F3/g/1_K562_RAD21_BX_rep1_hg19_CTCF_Occupancy_1000bp_read1-MIN100_Normalized.out
cp $LIBRARY/CTCF_Occupancy_1000bp/Composites/K562_RAD21_BX_rep1_hg19_CTCF_Occupancy_1000bp_read2-MIN100_Normalized.out F3/g/2_K562_RAD21_BX_rep1_hg19_CTCF_Occupancy_1000bp_read2-MIN100_Normalized.out
cp $LIBRARY/CTCF_Occupancy_1000bp/Composites/K562_SMC3_BX_rep1_hg19_CTCF_Occupancy_1000bp_read1-MIN100_Normalized.out F3/g/1_K562_SMC3_BX_rep1_hg19_CTCF_Occupancy_1000bp_read1-MIN100_Normalized.out
cp $LIBRARY/CTCF_Occupancy_1000bp/Composites/K562_SMC3_BX_rep1_hg19_CTCF_Occupancy_1000bp_read2-MIN100_Normalized.out F3/g/2_K562_SMC3_BX_rep1_hg19_CTCF_Occupancy_1000bp_read2-MIN100_Normalized.out
cp $LIBRARY/CTCF_Occupancy_1000bp/Composites/K562_-_BI_rep1_hg19_CTCF_Occupancy_1000bp_midpoint.out F3/g/3_K562_-_BI_rep1_hg19_CTCF_Occupancy_1000bp_midpoint.out
# Filter min?


# ===============================================================================================================================

[ -d F4 ] || mkdir F4

[ -d F4/a ] || mkdir F4/a
cp $LIBRARY/WebLogos/WDR5_logo.eps F4/a/

[ -d F4/b ] || mkdir F4/b
# peak-align analysis
# make bar plots

[ -d F4/c ] || mkdir F4/c
cp $LIBRARY/TSS_DistWDR5_Filter-SameStrand_1000bp/Composites/K562_Pol2_BX_rep1_hg19_TSS_DistWDR5_Filter-SameStrand_1000bp_read1.out F4/c/1_K562_Pol2_BX_rep1_hg19_TSS_DistWDR5_Filter-SameStrand_1000bp_read1.out
cp $LIBRARY/TSS_DistWDR5_Filter-SameStrand_1000bp/Composites/K562_RBBP5_BX_rep1_hg19_TSS_DistWDR5_Filter-SameStrand_1000bp_read1.out F4/c/2_K562_RBBP5_BX_rep1_hg19_TSS_DistWDR5_Filter-SameStrand_1000bp_read1.out
cp $LIBRARY/TSS_DistWDR5_Filter-SameStrand_1000bp/Composites/K562_WDR5_BX_rep1_hg19_TSS_DistWDR5_Filter-SameStrand_1000bp_read1.out F4/c/3_K562_WDR5_BX_rep1_hg19_TSS_DistWDR5_Filter-SameStrand_1000bp_read1.out
# cp $LIBRARY/TSS_DistWDR5_Filter-SameStrand_1000bp/Composites/K562_IgG_BX_merge_hg19_TSS_DistWDR5_Filter-SameStrand_1000bp_read1.out F4/c/4_K562_IgG_BX_merge_hg19_TSS_DistWDR5_Filter-SameStrand_1000bp_read1.out
cp $LIBRARY/TSS_DistWDR5_Filter-SameStrand_1000bp/Composites/K562_H3K4me3_BX_rep1_hg19_TSS_DistWDR5_Filter-SameStrand_1000bp_midpoint.out F4/c/5_K562_H3K4me3_BX_repq_hg19_TSS_DistWDR5_Filter-SameStrand_1000bp_midpoint.out

[ -d F4/d ] || mkdir F4/d
cp $LIBRARY/WDR5_Occupancy_1000bp/Composites/K562_-_BI_rep1_hg19_WDR5_Occupancy_1000bp_read2.out F4/d/1_K562_-_BI_rep1_hg19_WDR5_Occupancy_1000bp_read2.out
# read2? or midpoint?

[ -d F4/e ] || mkdir F4/e
cp $LIBRARY/WDR5_Occupancy_1000bp/Composites/K562_RBBP5_BX_rep1_hg19_WDR5_Occupancy_1000bp_read1_Normalized.out F4/e/1-RBBP5_K562_RBBP5_BX_rep1_hg19_WDR5_Occupancy_1000bp_read1_Normalized.out
cp $LIBRARY/WDR5_Occupancy_1000bp/Composites/K562_RBBP5_BX_rep1_hg19_WDR5_Occupancy_1000bp_read2_Normalized.out F4/e/2-RBBP5_K562_RBBP5_BX_rep1_hg19_WDR5_Occupancy_1000bp_read2_Normalized.out
cp $LIBRARY/WDR5_Occupancy_1000bp/Composites/K562_WDR5_BX_rep1_hg19_WDR5_Occupancy_1000bp_read1_Normalized.out F4/e/1-WDR5_K562_WDR5_BX_rep1_hg19_WDR5_Occupancy_1000bp_read1_Normalized.out
cp $LIBRARY/WDR5_Occupancy_1000bp/Composites/K562_WDR5_BX_rep1_hg19_WDR5_Occupancy_1000bp_read2_Normalized.out F4/e/2-WDR5_K562_WDR5_BX_rep1_hg19_WDR5_Occupancy_1000bp_read2_Normalized.out
cp $LIBRARY/WDR5_Occupancy_1000bp/Composites/K562_IgG_BX_merge_hg19_WDR5_Occupancy_1000bp_read1_Normalized.out F4/e/3_K562_IgG_BX_merge_hg19_WDR5_Occupancy_1000bp_read1_Normalized.out
cp $LIBRARY/WDR5_Occupancy_1000bp/Composites/K562_H3K4me3_BX_rep1_hg19_WDR5_Occupancy_1000bp_midpoint.out F4/e/4_K562_H3K4me3_BX_rep1_hg19_WDR5_Occupancy_1000bp_midpoint.out
# Figure out zoom in--extra? Filter min?

[ -d F4/f ] || mkdir F4/f
cp $LIBRARY/WDR5_Occupancy_1000bp/Composites/K562_-_BI_rep1_hg19_WDR5_Occupancy_1000bp_read2.out F4/f/1_K562_-_BI_rep1_hg19_WDR5_Occupancy_1000bp_read2.out

[ -d F4/g ] || mkdir F4/g
cp $LIBRARY/WDR5_Occupancy_1000bp/Composites/K562_WDR5_BX_rep1_hg19_WDR5_Occupancy_1000bp_read2_Normalized.out F4/g/1_K562_WDR5_BX_rep1_hg19_WDR5_Occupancy_1000bp_read2_Normalized.out
cp $LIBRARY/WDR5_Occupancy_1000bp/Composites/K562_-_BI_rep1_hg19_WDR5_Occupancy_1000bp_read2.out F4/g/2_K562_-_BI_rep1_hg19_WDR5_Occupancy_1000bp_read2.out
# Filter min?

[ -d F4/i ] || mkdir F4/i
cp $LIBRARY/WDR5_Occupancy_1000bp/Composites/K562_WDR5_BX_rep1_hg19_WDR5_Occupancy_1000bp_read2_Normalized.out F4/i/1_K562_WDR5_BX_rep1_hg19_WDR5_Occupancy_1000bp_read2_Normalized.out
cp $LIBRARY/WDR5_Occupancy_1000bp/Composites/K562_-_BI_rep1_hg19_WDR5_Occupancy_1000bp_read2.out F4/i/2_K562_-_BI_rep1_hg19_WDR5_Occupancy_1000bp_read2.out
cp $LIBRARY/WDR5_Occupancy_1000bp/Composites/K562_H3K4me3_BX_rep1_hg19_WDR5_Occupancy_1000bp_midpoint.out F4/i/3_K562_H3K4me3_BX_rep1_hg19_WDR5_Occupancy_1000bp_midpoint.out

[ -d F4/j ] || mkdir F4/j
cp $LIBRARY/WDR5_Occupancy_1000bp/Composites/K562_WDR5_BX_rep1_hg19_WDR5_Occupancy_1000bp_read1_Normalized.out F4/j/WDR5_K562_WDR5_BX_rep1_hg19_WDR5_Occupancy_1000bp_read1_Normalized.out
cp $LIBRARY/WDR5_Occupancy_1000bp/Composites/K562_Pol2_BX_rep1_hg19_WDR5_Occupancy_1000bp_read1_Normalized.out F4/j/Pol2_K562_Pol2_BX_rep1_hg19_WDR5_Occupancy_1000bp_read1_Normalized.out


# ===============================================================================================================================

[ -d F5 ] || mkdir F5

[ -d F5/a ] || mkdir F5/a
cp $LIBRARY/K562_CoPRO-expressed_Gene-refSeqTSS_1000bp/Composites/K562_GABPA_BX_rep1_hg19_K562_CoPRO-expressed_Gene-refSeqTSS_1000bp_read1.out F5/a/1_K562_GABPA_BX_rep1_hg19_K562_CoPRO-expressed_Gene-refSeqTSS_1000bp_read1.out
cp $LIBRARY/K562_CoPRO-expressed_Gene-refSeqTSS_1000bp/Composites/K562_SP1_BX_rep1_hg19_K562_CoPRO-expressed_Gene-refSeqTSS_1000bp_read1.out F5/a/2_K562_SP1_BX_rep1_hg19_K562_CoPRO-expressed_Gene-refSeqTSS_1000bp_read1.out
cp $LIBRARY/K562_CoPRO-expressed_Gene-refSeqTSS_1000bp/Composites/K562_IgG_BX_merge_hg19_K562_CoPRO-expressed_Gene-refSeqTSS_1000bp_read1.out F5/a/3_K562_IgG_BX_merge_hg19_K562_CoPRO-expressed_Gene-refSeqTSS_1000bp_read1.out
cp $LIBRARY/K562_CoPRO-expressed_Gene-refSeqTSS_1000bp/Composites/K562_EP300_BX_rep1_hg19_K562_CoPRO-expressed_Gene-refSeqTSS_1000bp_read1.out F5/a/4_K562_EP300_BX_rep1_hg19_K562_CoPRO-expressed_Gene-refSeqTSS_1000bp_read1.out
cp $LIBRARY/K562_CoPRO-expressed_Gene-refSeqTSS_1000bp/Composites/K562_TBP_BX_rep1_hg19_K562_CoPRO-expressed_Gene-refSeqTSS_1000bp_read1.out F5/a/5_K562_TBP_BX_rep1_hg19_K562_CoPRO-expressed_Gene-refSeqTSS_1000bp_read1.out
cp $LIBRARY/K562_CoPRO-expressed_Gene-refSeqTSS_1000bp/Composites/K562_GTF2B_BX_rep1_hg19_K562_CoPRO-expressed_Gene-refSeqTSS_1000bp_read1.out F5/a/6_K562_GTF2B_BX_rep1_hg19_K562_CoPRO-expressed_Gene-refSeqTSS_1000bp_read1.out
cp $LIBRARY/K562_CoPRO-expressed_Gene-refSeqTSS_1000bp/Composites/K562_Pol2_BX_rep1_hg19_K562_CoPRO-expressed_Gene-refSeqTSS_1000bp_read1.out F5/a/7_K562_Pol2_BX_rep1_hg19_K562_CoPRO-expressed_Gene-refSeqTSS_1000bp_read1.out
cp $LIBRARY/K562_CoPRO-expressed_Gene-refSeqTSS_1000bp/Composites/K562_H3K4me3_BX_rep1_hg19_K562_CoPRO-expressed_Gene-refSeqTSS_1000bp_midpoint.out F5/a/8_K562_H3K4me3_BX_rep1_hg19_K562_CoPRO-expressed_Gene-refSeqTSS_1000bp_midpoint.out
cp $LIBRARY/K562_CoPRO-expressed_Gene-refSeqTSS_1000bp/Composites/K562_-_BI_rep1_hg19_K562_CoPRO-expressed_Gene-refSeqTSS_1000bp_midpoint.out F5/a/9_K562_-_BI_rep1_hg19_K562_CoPRO-expressed_Gene-refSeqTSS_1000bp_midpoint.out

[ -d F5/b ] || mkdir F5/b
cp $LIBRARY/GABPA_Occupancy_1000bp/Composites/K562_GABPA_BX_rep1_hg19_GABPA_Occupancy_1000bp_read1_Normalized.out F5/b/1_K562_GABPA_BX_rep1_hg19_GABPA_Occupancy_1000bp_read1_Normalized.out
cp $LIBRARY/GABPA_Occupancy_1000bp/Composites/K562_SP1_BX_rep1_hg19_GABPA_Occupancy_1000bp_read1_Normalized.out F5/b/2_K562_SP1_BX_rep1_hg19_GABPA_Occupancy_1000bp_read1_Normalized.out
cp $LIBRARY/GABPA_Occupancy_1000bp/Composites/K562_EP300_BX_rep1_hg19_GABPA_Occupancy_1000bp_read1_Normalized.out F5/b/3_K562_EP300_BX_rep1_hg19_GABPA_Occupancy_1000bp_read1_Normalized.out
cp $LIBRARY/GABPA_Occupancy_1000bp/Composites/K562_IgG_BX_merge_hg19_GABPA_Occupancy_1000bp_read1_Normalized.out F5/b/4_K562_IgG_BX_merge_hg19_GABPA_Occupancy_1000bp_read1_Normalized.out

[ -d F5/c ] || mkdir F5/c
cp $LIBRARY/GABPA_Occupancy_1000bp/Composites/K562_GABPA_BX_rep1_hg19_GABPA_Occupancy_1000bp_read1_Normalized.out F5/c/1_K562_GABPA_BX_rep1_hg19_GABPA_Occupancy_1000bp_read1_Normalized.out
cp $LIBRARY/GABPA_Occupancy_1000bp/Composites/K562_TBP_BX_rep1_hg19_GABPA_Occupancy_1000bp_read1_Normalized.out F5/c/2_K562_TBP_BX_rep1_hg19_GABPA_Occupancy_1000bp_read1_Normalized.out
cp $LIBRARY/GABPA_Occupancy_1000bp/Composites/K562_GTF2B_BX_rep1_hg19_GABPA_Occupancy_1000bp_read1_Normalized.out F5/c/3_K562_GTF2B_BX_rep1_hg19_GABPA_Occupancy_1000bp_read1_Normalized.out

[ -d F5/d ] || mkdir F5/d
cp $LIBRARY/GABPA_Occupancy_1000bp/Composites/K562_GABPA_BX_rep1_hg19_GABPA_Occupancy_1000bp_read1_Normalized.out F5/d/1_K562_GABPA_BX_rep1_hg19_GABPA_Occupancy_1000bp_read1_Normalized.out
cp $LIBRARY/GABPA_Occupancy_1000bp/Composites/K562_Pol2_BX_rep1_hg19_GABPA_Occupancy_1000bp_read1_Normalized.out F5/d/2_K562_Pol2_BX_rep1_hg19_GABPA_Occupancy_1000bp_read1_Normalized.out
cp $LIBRARY/GABPA_Occupancy_1000bp/Composites/K562_H3K4me3_BX_rep1_hg19_GABPA_Occupancy_1000bp_read1.out F5/d/3_K562_H3K4me3_BX_rep1_hg19_GABPA_Occupancy_1000bp_read1.out

[ -d F5/e ] || mkdir F5/e
cp $LIBRARY/SP1_Occupancy_1000bp/Composites/K562_SP1_BX_rep1_hg19_SP1_Occupancy_1000bp_read1_Normalized.out F5/e/1_K562_SP1_BX_rep1_hg19_SP1_Occupancy_1000bp_read1_Normalized.out
cp $LIBRARY/SP1_Occupancy_1000bp/Composites/K562_GABPA_BX_rep1_hg19_SP1_Occupancy_1000bp_read1_Normalized.out F5/e/2_K562_GABPA_BX_rep1_hg19_SP1_Occupancy_1000bp_read1_Normalized.out
cp $LIBRARY/SP1_Occupancy_1000bp/Composites/K562_EP300_BX_rep1_hg19_SP1_Occupancy_1000bp_read1_Normalized.out F5/e/3_K562_EP300_BX_rep1_hg19_SP1_Occupancy_1000bp_read1_Normalized.out
cp $LIBRARY/SP1_Occupancy_1000bp/Composites/K562_IgG_BX_merge_hg19_SP1_Occupancy_1000bp_read1_Normalized.out F5/e/4_K562_IgG_BX_merge_hg19_SP1_Occupancy_1000bp_read1_Normalized.out

[ -d F5/f ] || mkdir F5/f
cp $LIBRARY/SP1_Occupancy_1000bp/Composites/K562_SP1_BX_rep1_hg19_SP1_Occupancy_1000bp_read1_Normalized.out F5/f/1_K562_SP1_BX_rep1_hg19_SP1_Occupancy_1000bp_read1_Normalized.out
cp $LIBRARY/SP1_Occupancy_1000bp/Composites/K562_TBP_BX_rep1_hg19_SP1_Occupancy_1000bp_read1_Normalized.out F5/f/2_K562_TBP_BX_rep1_hg19_SP1_Occupancy_1000bp_read1_Normalized.out
cp $LIBRARY/SP1_Occupancy_1000bp/Composites/K562_GTF2B_BX_rep1_hg19_SP1_Occupancy_1000bp_read1_Normalized.out F5/f/3_K562_GTF2B_BX_rep1_hg19_SP1_Occupancy_1000bp_read1_Normalized.out

[ -d F5/g ] || mkdir F5/g
cp $LIBRARY/SP1_Occupancy_1000bp/Composites/K562_SP1_BX_rep1_hg19_SP1_Occupancy_1000bp_read1_Normalized.out F5/g/1_K562_SP1_BX_rep1_hg19_SP1_Occupancy_1000bp_read1_Normalized.out
cp $LIBRARY/SP1_Occupancy_1000bp/Composites/K562_Pol2_BX_rep1_hg19_SP1_Occupancy_1000bp_read1_Normalized.out F5/g/2_K562_Pol2_BX_rep1_hg19_SP1_Occupancy_1000bp_read1_Normalized.out
cp $LIBRARY/SP1_Occupancy_1000bp/Composites/K562_H3K4me3_BX_rep1_hg19_SP1_Occupancy_1000bp_read1.out F5/g/3_K562_H3K4me3_BX_rep1_hg19_SP1_Occupancy_1000bp_read1.out




# ===============================================================================================================================

[ -d S1 ] || mkdir S1

[ -d S1/a ] || mkdir S1/a
cp $LIBRARY/WebLogos/scrambledNFIA_logo.eps S1/a/

[ -d S1/b ] || mkdir S1/b
cp $LIBRARY/scrambledNFIA_NucSort_1000bp/FourColor/scrambledNFIA_NucSort_1000bp_31bp.svg S1/b/
cp $LIBRARY/scrambledNFIA_NucSort_1000bp/SVG/K562_NFIA_BX_rep1_hg19_scrambledNFIA_NucSort_1000bp_read1_Normalized_merge_label.svg S1/b/
cp $LIBRARY/scrambledNFIA_NucSort_1000bp/SVG/K562_-_BI_rep1_hg19_scrambledNFIA_NucSort_1000bp_midpoint_combined_treeview_label.svg S1/b/

[ -d S1/c ] || mkdir S1/c
cp $LIBRARY/scrambledNFIA_NucSort_1000bp/Composites/K562_-_BI_rep1_hg19_scrambledNFIA_NucSort_1000bp_read2.out S1/c/1_K562_-_BI_rep1_hg19_scrambledNFIA_NucSort_1000bp_read2.out


# ===============================================================================================================================


[ -d S2 ] || mkdir S2

[ -d S2/a ] || mkdir S2/a
cp $LIBRARY/scrambledCTCF_NucSort_1000bp/Composites/primates.phyloP46way_scrambledCTCF_NucSort_1000bp.out S2/a/1_primates.phyloP46way_scrambledCTCF_NucSort_1000bp.out
cp $LIBRARY/scrambledCTCF_NucSort_1000bp/Composites/dbSnp153_snv_scrambledCTCF_NucSort_1000bp.out S2/a/2_dbSnp153_snv_scrambledCTCF_NucSort_1000bp.out

[ -d S2/b ] || mkdir S2/b
cp $LIBRARY/WebLogos/scrambledCTCF_logo.eps S2/b/

[ -d S2/c ] || mkdir S2/c
cp $LIBRARY/scrambledCTCF_NucSort_1000bp/FourColor/scrambledCTCF_NucSort_1000bp_31bp.svg S2/c/
cp $LIBRARY/scrambledCTCF_NucSort_1000bp/SVG/K562_CTCF_BX_rep1_hg19_scrambledCTCF_NucSort_1000bp_read1_Normalized_merge_label.svg S2/c/
cp $LIBRARY/scrambledCTCF_NucSort_1000bp/SVG/K562_-_BI_rep1_hg19_scrambledCTCF_NucSort_1000bp_midpoint_combined_treeview_label.svg S2/c/

[ -d S2/d ] || mkdir S2/d
cp $LIBRARY/scrambledCTCF_NucSort_1000bp/Composites/K562_-_BI_rep1_hg19_scrambledCTCF_NucSort_1000bp_read2.out S2/d/1_K562_-_BI_rep1_hg19_scrambledCTCF_NucSort_1000bp_read2.out

# ===============================================================================================================================

[ -d S3 ] || mkdir S3

[ -d S3/a ] || mkdir S3/a
cp $LIBRARY/WDR5_Occupancy_1000bp/FourColor/WDR5_Occupancy_1000bp_31bp.svg S3/a/
cp $LIBRARY/WDR5_Occupancy_1000bp/SVG/K562_WDR5_BX_rep1_hg19_WDR5_Occupancy_1000bp_read1_Normalized_merge_label.svg S3/a/
cp $LIBRARY/WDR5_Occupancy_1000bp/SVG/K562_RBBP5_BX_rep1_hg19_WDR5_Occupancy_1000bp_read1_Normalized_merge_label.svg S3/a/
cp $LIBRARY/WDR5_Occupancy_1000bp/SVG/K562_IgG_BX_merge_hg19_WDR5_Occupancy_1000bp_read1_Normalized_merge_label.svg S3/a/


[ -d S3/b ] || mkdir S3/b
# build barplot and pie charts


[ -d S3/c ] || mkdir S3/c
cp $LIBRARY/WebLogos/WDR5_logo.eps S3/c/
cp $LIBRARY/WDR5_Occupancy_1000bp/Composites/primates.phyloP46way_WDR5_Occupancy_1000bp.out S3/c/1_primates.phyloP46way_WDR5_Occupancy_1000bp.out
cp $LIBRARY/WDR5_Occupancy_1000bp/Composites/dbSnp153_snv_WDR5_Occupancy_1000bp.out S3/c/2_dbSnp153_snv_WDR5_Occupancy_1000bp.out

# ===============================================================================================================================

[ -d S4 ] || mkdir S4

[ -d S4/a ] || mkdir S4/a
cp $LIBRARY/WebLogos/GABPA_logo.eps S4/a/

[ -d S4/b ] || mkdir S4/b/
cp $LIBRARY/GABPA_Occupancy_1000bp/FourColor/GABPA_Occupancy_1000bp_31bp.svg S4/b/
cp $LIBRARY/GABPA_Occupancy_1000bp/SVG/K562_GABPA_BX_rep1_hg19_GABPA_Occupancy_1000bp_read1_Normalized_merge_label.svg S4/b/
cp $LIBRARY/GABPA_Occupancy_1000bp/SVG/K562_IgG_BX_merge_hg19_GABPA_Occupancy_1000bp_read1_Normalized_merge_label.svg S4/b/

[ -d S4/c ] || mkdir S4/c
cp $LIBRARY/WebLogos/SP1_logo.eps S4/c/

[ -d S4/d ] || mkdir S4/d
cp $LIBRARY/SP1_Occupancy_1000bp/FourColor/SP1_Occupancy_1000bp_31bp.svg S4/d/
cp $LIBRARY/SP1_Occupancy_1000bp/SVG/K562_SP1_BX_rep1_hg19_SP1_Occupancy_1000bp_read1_Normalized_merge_label.svg S4/d/
cp $LIBRARY/SP1_Occupancy_1000bp/SVG/K562_IgG_BX_merge_hg19_SP1_Occupancy_1000bp_read1_Normalized_merge_label.svg S4/d/
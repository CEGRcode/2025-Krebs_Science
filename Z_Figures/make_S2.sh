#!/bin/bash

# Organize data from X_Bulk_Processing into Z_Figures for F4

### CHANGE ME
WRK=/storage/group/bfp2/default/hxc585_HainingChen/2025_Chen_TF-Nuc/Z_Figures
#WRK=/storage/home/owl5022/scratch/2024-Krebs_Science/Z_Figures
###

set -exo
module load anaconda3
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

LIBRARY=$WRK/../X_Bulk_Processing/Library
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar

[ -d S2 ] || mkdir S2
[ -d S2/a ] || mkdir S2/a
## check with Frank that NRF1 and JUND SP1 SPI1 ZKSCAN1 ZNF263 are redundent
# ===============================================================================================================================

cp $LIBRARY/WebLogos/ATF2_MA1632.2_logo.eps	 S2/a
cp $LIBRARY/WebLogos/ATF3_MA0605.3_logo.eps S2/a
cp $LIBRARY/WebLogos/ATF4_MA0833.3_logo.eps S2/a
cp $LIBRARY/WebLogos/ATF7_MA0834.1_logo.eps S2/a
cp $LIBRARY/WebLogos/BHLHE40_MA0464.3_logo.eps	 S2/a
cp $LIBRARY/WebLogos/CREB1_MA0018.1_logo.eps S2/a
cp $LIBRARY/WebLogos/CREM_MA0609.3_logo.eps	 S2/a
cp $LIBRARY/WebLogos/CTCF_MA0139.1_logo.eps	 S2/a
cp $LIBRARY/WebLogos/EGR1_MA0162.5_logo.eps	 S2/a
cp $LIBRARY/WebLogos/ELF1_MA0473.2_logo.eps	 S2/a
cp $LIBRARY/WebLogos/ELF4_MA0641.1_logo.eps	 S2/a
cp $LIBRARY/WebLogos/ELK1_MA0028.3_logo.eps	 S2/a
cp $LIBRARY/WebLogos/ESRRA_MA0592.1_logo.eps S2/a
cp $LIBRARY/WebLogos/JUND_MA0492.2_logo.eps S2/a
cp $LIBRARY/WebLogos/MAFF_MA0495.2_logo.eps	S2/a
cp $LIBRARY/WebLogos/MAFK_MA0496.2_logo.eps	S2/a
cp $LIBRARY/WebLogos/MAX_MA0058.3_logo.eps	S2/a
cp $LIBRARY/WebLogos/MAZ_MA1522.1_logo.eps	S2/a
cp $LIBRARY/WebLogos/MEIS2_MA1640.1_logo.eps S2/a
cp $LIBRARY/WebLogos/MITF_MA0620.2_logo.eps S2/a
cp $LIBRARY/WebLogos/MYC_MA0147.4_logo.eps	S2/a
cp $LIBRARY/WebLogos/NFIC_MA1527.2_logo.eps S2/a
cp $LIBRARY/WebLogos/NR2F1_MA1537.1_logo.eps 	S2/a
cp $LIBRARY/WebLogos/NRF1_MA0506.3_logo.eps S2/a
cp $LIBRARY/WebLogos/PKNOX1_MA0782.3_logo.eps	S2/a
cp $LIBRARY/WebLogos/REST_MA0138.3_logo.eps		S2/a
cp $LIBRARY/WebLogos/RFX1_MA0509.3_logo.eps	S2/a
cp $LIBRARY/WebLogos/SP1_MA0079.5_logoRC.eps	S2/a
cp $LIBRARY/WebLogos/SPI1_MA0080.7_logo.eps	S2/a
cp $LIBRARY/WebLogos/USF2_MA0526.5_logo.eps		S2/a
cp $LIBRARY/WebLogos/ZBTB33_MA0527.1_logo.eps	S2/a
cp $LIBRARY/WebLogos/ZNF24_MA1124.1_logo.eps	S2/a	
cp $LIBRARY/WebLogos/ZNF263_MA0528.1_logo.eps	S2/a


# Composites
cp $LIBRARY/BI_Pileups/MA1632/Composites/*.out S2/a
cp $LIBRARY/BI_Pileups/MA0605/Composites/*.out S2/a
cp $LIBRARY/BI_Pileups/MA0833/Composites/*.out S2/a
cp $LIBRARY//BI_Pileups/MA0834/Composites/*.out S2/a
cp $LIBRARY/BI_Pileups/MA0464/Composites/*.out S2/a
cp $LIBRARY/BI_Pileups/MA0018/Composites/*.out S2/a
cp $LIBRARY/BI_Pileups/MA0609/Composites/*.out S2/a
cp $LIBRARY/BI_Pileups/MA0138/Composites/*.out S2/a
cp $LIBRARY/BI_Pileups/MA0162/Composites/*.out S2/a
cp $LIBRARY/BI_Pileups/MA0473/Composites/*.out S2/a
cp $LIBRARY/BI_Pileups/MA0641/Composites/*.out S2/a
cp $LIBRARY/BI_Pileups/MA0028/Composites/*.out S2/a
cp $LIBRARY/BI_Pileups/MA0592/Composites/*.out S2/a
cp $LIBRARY/BI_Pileups/MA0492/Composites/*.out S2/a
cp $LIBRARY/BI_Pileups/MA0495/Composites/*.out S2/a
cp $LIBRARY/BI_Pileups/MA0496/Composites/*.out S2/a
cp $LIBRARY/BI_Pileups/MA0058/Composites/*.out S2/a
cp $LIBRARY/BI_Pileups/MA1522/Composites/*.out S2/a
cp $LIBRARY/BI_Pileups/MA1640/Composites/*.out S2/a
cp $LIBRARY/BI_Pileups/MA0620/Composites/*.out S2/a
cp $LIBRARY/BI_Pileups/MA0147/Composites/*.out S2/a
cp $LIBRARY//BI_Pileups/MA1527/Composites/*.out S2/a
cp $LIBRARY/BI_Pileups/MA1537/Composites/*.out S2/a
cp $LIBRARY/BI_Pileups/MA0509/Composites/*.out S2/a
cp $LIBRARY/BI_Pileups/MA0782/Composites/*.out S2/a
cp $LIBRARY/BI_Pileups/MA0138/Composites/*.out S2/a
cp $LIBRARY/BI_Pileups/MA0509/Composites/*.out S2/a
cp $LIBRARY/BI_Pileups/MA0079/Composites/*.out S2/a
cp $LIBRARY/BI_Pileups/MA0080/Composites/*.out S2/a
cp $LIBRARY/BI_Pileups/MA0526/Composites/*.out S2/a
cp $LIBRARY/BI_Pileups/MA0527/Composites/*.out S2/a
cp $LIBRARY//BI_Pileups/MA1124/Composites/*.out S2/a
cp $LIBRARY/BI_Pileups/MA0528/Composites/*.out S2/a

# ===============================================================================================================================
[ -d S2 ] || mkdir S2
[ -d S2/b ] || mkdir S2/b
cp $LIBRARY/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile1_1000bp/Composites/*.out S2/b
cp $LIBRARY/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile2_1000bp/Composites/*.out S2/b
cp $LIBRARY/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile3_1000bp/Composites/*.out S2/b
cp $LIBRARY/CTCF_MA1929.1_SORT-TFnucRatio_GROUP-Quartile4_1000bp/Composites/*.out S2/b
cp $LIBRARY/MA1929_1_mm10_intersected_MPE-seq10min_164bp_category1_1000bp/Composites/*.out S2/b
cp $LIBRARY/MA1929_1_mm10_intersected_MPE-seq10min_164bp_category2_1000bp/Composites/*.out S2/b
cp $LIBRARY/MA1929_1_mm10_intersected_MPE-seq10min_164bp_category3_1000bp/Composites/*.out S2/b
cp $LIBRARY/MA1929_1_mm10_intersected_MPE-seq10min_164bp_category4_1000bp/Composites/*.out S2/b


# ===============================================================================================================================
[ -d S2 ] || mkdir S2
[ -d S2/c ] || mkdir S2/c
cp $LIBRARY/ATF2_MA1632.2_K562-specific-Unbound_1000bp/Composites/*.out S2/c
cp  $LIBRARY/ATF2_MA1632.2_Unbound_1000bp/Composites/*.out S2/c
cp  $LIBRARY/CREM_MA0609.3_K562-specific-Unbound_1000bp/Composites/*.out S2/c
cp  $LIBRARY/CREM_MA0609.3_Unbound_1000bp/Composites/*.out S2/c
cp  $LIBRARY/CTCF_MA1929.1_K562-specific-Unbound_1000bp/Composites/*.out S2/c
cp  $LIBRARY/CTCF_MA1929.1_Unbound_1000bp/Composites/*.out S2/c
cp  $LIBRARY/ZNF263_MA0528.1_K562-specific-Unbound_1000bp/Composites/*.out S2/c
cp  $LIBRARY/ZNF263_MA0528.1_Unbound_1000bp/Composites/*.out S2/c
cp  $LIBRARY/ZKSCAN1_MA1585.1_K562-specific-Unbound_1000bp/Composites/*.out S2/c
cp  $LIBRARY/ZKSCAN1_MA1585.1_Unbound_1000bp/Composites/*.out S2/c
	









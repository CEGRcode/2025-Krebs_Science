#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=48gb
#SBATCH -t 3:00:00
#SBATCH -A open
#SBATCH -o logs/9_MNase_Peak_Calling.log.out
#SBATCH -e logs/9_MNase_Peak_Calling.log.err

# Perform peak calling of various nucleosome particle sized fragment populations.
# see New_scripts_JEK_240731/01b_shift_scIDX_240709_v1_MNase.sh
# see New_scripts_JEK_240731/04_MNase_peaks__240710.sh

# data/RefPT-Krebs
#   |--MNase-Nucleosomes.bed

# Dependencies
# - java
# - python2
# - numpy (py2)

set -exo
source activate genetrack

# Inputs and outputs
BAMFILE=../data/BAM/MNase-seq_ENCODE_merge_hg38.bam
BLACKLIST=../data/hg38_files/ENCFF356LFX_hg38_exclude.bed
KREBS=../data/RefPT-Krebs

# Script shortcuts
GENETRACK=../bin/genetrack_v2.py
PICARD=../bin/picard.jar
SCRIPTMANAGER=../bin/ScriptManager-v0.15.jar

# Create output directories if they don't exist
[ -d logs ] || mkdir logs
[ -d SCIDX ] || mkdir SCIDX
[ -d mnase ] || mkdir mnase
[ -d AllParticles ] || mkdir AllParticles

# BAM to scIDX
echo "BAM to scIDX..."
time1=$(date +%s)
java -jar $SCRIPTMANAGER bam-format-converter bam-to-scidx -1 --shift 82 $BAMFILE -o SCIDX/mnase.tab
time2=$(date +%s)
echo "Elapsed Time (mnase): $(($time2-$time1)) seconds"
# X1 sec

# Reformat combined value to "forward" and zero out "reverse"
echo "Reformat..."
time1=$(date +%s)
awk 'BEGIN {getline; print $0; OFS="\t"; FS="\t";}{
	print $1, $2, $3+$4, "0", $5
}' SCIDX/mnase.tab > mnase/formatted.tab
time2=$(date +%s)
echo "Elapsed Time (mnase): $(($time2-$time1)) seconds"
# X2 sec

# GeneTrack
echo "GeneTrack..."
time1=$(date +%s)
python2 $GENETRACK -s 40 -e 80 -F 32 mnase/formatted.tab
time2=$(date +%s)
echo "Elapsed Time (mnase): $(($time2-$time1)) seconds"
# X3 sec

# Convert to BED
echo "Convert GFF to BED..."
time1=$(date +%s)
java -jar $SCRIPTMANAGER coordinate-manipulation gff-to-bed mnase/genetrack_s40e80F32/formatted_s40e80F32.gff -o mnase/genetrack_output.bed
time2=$(date +%s)
echo "Elapsed Time (mnase): $(($time2-$time1)) seconds"
# X4 sec

# Expand BED (max DNA frag size for indicated nuc peak except supraoctasomes expanded to 200bp)
echo "Expand BED..."
time1=$(date +%s)
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 164 mnase/genetrack_output.bed -o mnase/expanded.bed
time2=$(date +%s)
echo "Elapsed Time (mnase): $(($time2-$time1)) seconds"
# X5 sec

# Filter to exclude lines expanded into negative value coordinate indexes & add 7th column with particle type indicated
echo "Filter out negative coordinates and build ID..."
time1=$(date +%s)
awk '{OFS="\t"}{FS="\t"}{if ($2>0) {$4="Subtetra-" NR "-" $1"_"$2"_"$3"_"$4; print;}}' mnase/expanded.bed > $KREBS/MNase-Nucleosomes.bed
time2=$(date +%s)
echo "Elapsed Time (mnase): $(($time2-$time1)) seconds"
# X6 sec

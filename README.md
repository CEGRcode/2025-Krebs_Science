# Genome-wide rotational and translational phasing of nucleosomes with human transcription factors

### Haining Chen<sup>1,4</sup>, Jordan E. Krebs<sup>2,4</sup>, Olivia W. Lang<sup>1</sup>, Jiayi Hu<sup>1</sup>, Devin C. Mellini<sup>1</sup>, Judith Hyle<sup>3</sup>, Chunliang Li<sup>3</sup>, William K. M. Lai<sup>1,5</sup>, B. Franklin Pugh<sup>1,5,6,*</sup>

<sup>1</sup>Department of Molecular Biology and Genetics, Cornell University, Ithaca, New York, 14853, USA.
<sup>2</sup>MD/PhD Medical Scientist Training Program, Penn State College of Medicine, Hershey, PA, USA.
<sup>3</sup>Department of Tumor Cell Biology, St. Jude Children’s Research Hospital, Memphis, TN 38105, USA
<sup>4</sup>These authors contributed equally
<sup>5</sup>Senior author
<sup>6</sup>Lead contact
<sup>*</sup>Correspondence: fp265@cornell.edu

### PMID : [XXXXXXXX](https://pubmed.ncbi.nlm.nih.gov/XXXXXXXX/)
### GEO ID : [GSE266547](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266547) [GSE267711](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267711)

## Abstract

How transcription factors (TFs) and their binding sites organize and engage nucleosomes at natural genomic locations remains poorly understood. Here we develop Benzonase-seq to measure the rotational phasing of nucleosomes in human cells, and enhance ChIP-exo (v6) to measure rotational phasing on the same DNA molecule bound by a TF. Unbound CTCF sites were found to be rotationally accessible on nucleosomes and this rotational accessibility is encoded by classical dinucleotide periodicities. CTCF binding results in nucleosome displacement to adjacent DNA phasing sequences. Examining 40 TF classes, their unbound sites were phased inward, outward, or lacked phasing. In all examined cases, TF binding (e.g. NFIA and FoxA) results in adjacent rotational and translational phasing, which is not dinucleotide encoded. Benzonase-seq also more robustly maps nucleosome and subnucleosome positions in hard-to-map CpG islands. These findings provide a clearer view of how TFs engage and position nucleosomes to shape the natural chromatin landscape.


## Directions
To recreate the figures for this manuscript, please execute the scripts in each directory in numerical order. Each directory's README includes more specific details on execution. To be more explicit, run the scripts in each directory in the following order: `00_Download_and_Preprocessing`, `01_Run_GenoPipe`, `02_Call_Nucleosomes`, `03_Call_JASPAR`, `04_Call_Motifs`, `X_Bulk_Processing`, and then finally `Z_Figures`.

## Dependencies

> [!NOTE]  
> The scripts in this repo `source` the conda environment names used here within the scripts so you will want to modify them to point to your own environments if you attempt to run these scripts yourself.

### conda bx

Use the following [anaconda](https://anaconda.org/) environment initialization for setting up dependencies

```sh
conda create -n bx -c bioconda -c conda-forge bedtools bowtie2 bowtie=1.2.3 bwa meme pybigwig pysam samtools kneed opencv openjdk scipy seaborn wget
conda activate bx
```

### conda genetrack

For genetrack-executing scripts, a python2 environment needed to be created. The create command for that env is as follows:

```sh
conda create -n genetrack -c conda-forge -c bioconda python=2.7 numpy bedtools
conda activate genetrack
```

### SRA Toolkit

You will also need to install the SRA toolkit to retrieve raw data from SRA (cannot be too old of a version since we used `fasterq-dump` unless you refactor our scripts to use the older, slower, `fastq-dump`). We did not use a conda install due to the inconsistent success we've had with using the `sra-tools` recipes. Don't forget to [configure the toolkit](https://github.com/ncbi/sra-tools/wiki/03.-Quick-Toolkit-Configuration) before using (`vdb-config -i`).

## Table of Contents

### 00_Download_and_Preprocessing
Perform the preprocessing steps including alignment of raw sequencing data from both novel and previously published data.

### 01_Run_GenoPipe
Perform quality control for genetic background on these data by running GenoPipe on the aligned BAMs.

### 02_Call_Nucleosomes
Call nucleosome positions and identify TSS and +1 nucleosome reference points with different sorts.

### 03_Call_JASPAR
Call JASPAR motifs and subset to "bound" sites using ENCODE peak data.

### 04_Call_Motifs
Build de novo sequence-specific transcription factor (ssTF) motif reference points using Benzonase ChIP-exo data.

### X_Bulk_Processing
With the BAM and BED files built from the scripts in the above directories, perform bulk read pileups for heatmaps and composites.

### Z_Figures
Copy/organize results from bulk processing into figure-specific directories corresponding to subfigures in the manuscript. Also includes custom/one-off scripts for analysis that didn't need bulk-style execution.

### AI_files
all figures in paper

### data
Store large files to be globally accessed by the scripts in each directory

### bin
Generalized scripts and executables for global access by each of the numbered directories.

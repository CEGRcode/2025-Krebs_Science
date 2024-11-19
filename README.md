# Genome-wide rotational and translational setting of transcription factors with nucleosomes

### Jordan E. Krebs<sup>&#x2020; 1,2</sup>, Haining Chen<sup>&#x2020; 2</sup>, Olivia W. Lang<sup>2</sup>, William K. M. Lai<sup>2</sup>, B. Franklin Pugh<sup>2</sup>

&#x2020; Co-first author

<sup>1</sup>MD/PhD Medical Scientist Training Program, Penn State College of Medicine, Hershey, PA, USA.

<sup>2</sup>Department of Molecular Biology and Genetics, Cornell University, Ithaca, New York, 14853, USA

### Correspondence:fp265@cornell.edu

### PMID : [XXXXXXXX](https://pubmed.ncbi.nlm.nih.gov/XXXXXXXX/)
### GEO ID : [XXXXXXXX](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=XXXXXXXX)

## Abstract
Genes are regulated by transcription factors (TFs) bound to DNA sites (TFBS). Unbound sites are often nucleosomal and inaccessible. Current assays cannot address whether TFs encounter phasing of the DNA helix on nucleosomal surfaces across a natural genomic context, and how this may change upon TF binding. Here we use an endonuclease, either alone or coupled to ChIP-exo to measure the genomic rotational and translational (positional) phasing of nucleosomal DNA at unbound and TF-bound TFBSs in human cells. Unbound sites had a preferred rotational phase, but generally lacked a translational phase. In contrast, TF/TFBS complexes were often engaged with an adjacent translationally and rotationally phased nucleosome. Thus, a few molecular themes may govern how TFs engage nucleosomes.


## Directions
To recreate the figures for this manuscript, please execute the scripts in each directory in numerical order. Each directory's README includes more specific details on execution. To be more explicit, run the scripts in each directory in the following order: `00_Download_and_Preprocessing`, `01_Run_GenoPipe`, `02_Call_Nucleosomes`, `03_Call_JASPAR`, `04_Call_Motifs`, `X_Bulk_Processing`, and then finally `Z_Figures`.

## Dependencies
Use the following [anaconda](https://anaconda.org/) environment initialization for setting up dependencies

```
conda create -n bx -c bioconda -c conda-forge bedtools bowtie2 bwa cutadapt meme opencv pandas samtools scipy sra-tools wget pybigwig
```

For genetrack-executing script, a python2 environment needed to be created. The create command for that env is as follows:

```
conda create -n genetrack -c conda-forge -c bioconda python=2.7 numpy
```


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

### data
Store large files to be globally accessed by the scripts in each directory

### bin
Generalized scripts and executables for global access by each of the numbered directories.

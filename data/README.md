
This directory stores PEGR data ( BAM,MEME files) and generally used reference files (BED, RefPT-Motif).

### data/BAM/
merged and renamed BAM files go here

#### data/BAM/NormalizationFactors/
Your normalization factors stored in `.txt` files go here.

### data/MEME/
MEME motif discovery results from BX K562 go here

### data/RefPT-Other/
Nucleosome and TSS centered reference files with various sorts and expansions

### data/RefPT-Motif/
Motif-centered reference files with various sorts and expansions.

### data/FASTQ/
Raw sequencing datasets go here.

### data/hg19_files
Run the following commands in the terminal from the `data/hg19_files` directory to download and index the reference genome. (Required to run alignments in `00_Download_and_Preprocessing`, strain checks in `01_Run_GenoPipe`, and for running other sequence analyses and figure generation)
```
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz -O hg19.fa.gz
gzip -d hg19.fa.gz         # hg19.fa
bwa index hg19.fa          # hg19.fa.amb, hg19.fa.ann, hg19.fa.bwt, hg19.fa.pac, hg19.fa.sa
samtools faidx hg19.fa     # hg19.fa.fai
```

...and there are some other reference files already provided for you within this repository.
```
hg19.info                   # used by ChExMix
hg19_exclude_contig.txt     # used by ChExMix
hg19_background_model.txt   # used by ChExMix
hg19_exclude.bed            # used by 00_Download_and_Preprocessing/3_normalize_samples.pbs
                            # and 02_Call_RefPT/2_FIMO_Motifs_from_Genome.pbs
```



This directory stores generally used scripts and executables.

`* Please note the marked binaries in bin that are specific for executing on Linux machines. Please overwrite with appropriate binaries if you are executing from a different OS (e.g. MacOS)`

### ScriptManager-v0.15.jar
Download the Java binary executable for ScriptManager that includes a collection of tools including TagPileup which is used to count tags and calculate coverage of samples around reference points.
```
wget https://github.com/CEGRcode/scriptmanager/releases/download/v0.15/ScriptManager-v0.15.jar
```

### picard.jar
Download the Java binary executable for Picard.
```
wget https://github.com/broadinstitute/picard/releases/download/2.27.3/picard.jar
```

### chexmix.v0.52.public.jar
Download the ChExMix binary executable for peak calling.
```
wget https://github.com/seqcode/chexmix/releases/download/v0.52/chexmix.v0.52.public.jar
```

### UCSC binaries
Download the appropriate binary for your OS from UCSC. (current binaries are for `linux.x86_64`)

#### bedGraphToBigWig (Linux)*
```
wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
chmod 755 bedGraphToBigWig
```

#### bigBedToBed (Linux)*
```
wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed
chmod 755 bigBedToBed
```

### calculate_BED_ScoreRatio.pl

Calculate a ratio between the score columns of two different BED files to create a new BED file with the ratio for the score.

```
usage:		perl calculate_BED_ScoreRatio.pl	Numerator_BED_File	Denominator_BED_File	Output_BED
Example:	perl calculate_BED_ScoreRatio.pl numerator.bed denominator.bed ratio.bed
	BED information inherited from numerator BED file
```

### dinucleotide_CDT_from_FASTA.py
Scan a FASTA file for positional dinucleotide content
```
python dinucleotide_CDT_from_FASTA.py  -h
usage: dinucleotide_CDT_from_FASTA.py [-h] -i fasta_fn -s dinucleotides_str -o tsv_fn

============
Get 0/1 matrix (CDT format) of dinucleotides
============

optional arguments:
 -h, --help            show this help message and exit
 -i fasta_fn, --input fasta_fn
                       the FASTA file to analyze
 -s dinucleotides_str, --seq dinucleotides_str
                       the "-" delimited set of dinucleotides to check for
 -o tsv_fn, --output tsv_fn
                       the output CDT formatted 0/1 matrix of dinucleotide matches
```

### kmer_tally_to_pwm.py
Can feed this script the tally output from `upstream_seq_tally.py` to summarize the kmer tallies' positional nucleotide content.
```
usage: kmer_tally_to_pwm.py [-h] -i bam_fn -o tsv_fn [-c COLUMN]
kmer_tally_to_pwm.py: error: the following arguments are required: -i/--input, -o/--output
```

### make_violin_plot.py
Make violin plots with presets formatted for this manuscript.
```
usage: make_violin_plot.py [-h] [-i two_col_file] [--width width] [--height height] [--title title] [--xlabel xlabel] [--ylabel ylabel] [--preset1] [--preset2] [-o output_svg]

optional arguments:
  -h, --help            show this help message and exit
  -i two_col_file, --input two_col_file
                        tab-delimited file made of two columns: first column y values to plot (must all be numeric values), second column is the grouping (which violin group along x-axis to contribute to)
  --width width         width of figure
  --height height       height of figure
  --title title         title of figure
  --xlabel xlabel       x-axis label
  --ylabel ylabel       y-axis label
  --preset1             use proximal/distal presets for nucleosome intervals (Fig 4d)
  --preset2             use proximal/distal presets for half-nucleosome intervals (Fig 4d)
  -o output_svg, --output output_svg
                        name of SVG filepath to save figure to (if none provided, figure pops up in new window)
```

### upstream_seq_tally.py
Tally up kmers upstream of each read (can use proper-pair flag for paired-end data). Reverse-stand mapped reads are reverse complemented to orient kmers with the read on the right.
```
usage: upstream_seq_tally.py [-h] -i bam_fn -g fasta_fn -o tsv_fn [-p]
                             [-k KMER]
upstream_seq_tally.py: error: the following arguments are required: -i/--input, -g/--genome, -o/--output
```

### resize_png.py
```
Usage:
This script will compress a PNG file using INTER_AREA interpolation

python resize_png.py -i <input image> -o <output image> -r <row num after compress> -c <col num after compress>'

Example:
python resize_png.py -i 17385_memeCluster1_10.png -o test.png -r 400 -c 400
```

### sum_Col_CDT.pl
This script sums the columns of a CDT matrix file by column values (CDT to composite).
```
usage:		perl sum_Col_CDT.pl	Input_CDT_File	Output_TAB_File
Example:	perl sum_Col_CDT.pl input.cdt composite.out
```

### sum_each_CDT.py
```
usage: sum_each_CDT.py [-h] -1 cdt_fn -2 cdt_fn -o cdt_fn

This script will element-wise sum two CDT matrices. Checks for matching YORF
and NAME Example: python sum_each_CDT.py -1 FIRST.cdt -2 SECOND.cdt -o SUM.cdt

optional arguments:
  -h, --help            show this help message and exit
  -1 cdt_fn, --file1 cdt_fn
                        the first CDT file to sum
  -2 cdt_fn, --file2 cdt_fn
                        the second CDT file to sum
  -o cdt_fn, --output cdt_fn
                        the summed CDT file
```

## convert_wig_to_bedgraph.py
```
python convert_wig_to_bedgraph.py -i test.wig -o test.bg
```

## pileup_BedGraph_on_RefPT.py
```
python pileup_BedGraph_on_RefPT.py -i test.bg -r test.bed -o blah.cdt
```

## pileup_BigWig_on_RefPT.py
```
python pileup_BigWig_on_RefPT.py -i test.bg -r test.bed -o blah.cdt
```

## pileup_BigWig_on_RefPT_stranded.py
```
python pileup_BigWig_on_RefPT_stranded.py -i test.bg -r test.bed -o blah.cdt
```

This directory stores generally used scripts and executables.

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

#### bigWigToBedGraph

#### bedGraphToBigWig

#### bigBedToBed

#### bedToBigBed


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
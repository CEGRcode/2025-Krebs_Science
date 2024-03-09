import sys
import argparse
import re


def getParams():
    '''Parse parameters from the command line'''

    parser = argparse.ArgumentParser(description='''
This script converts a wiggle-formatted file to a BedGraph format.

Example: python convert_wig_to_scidx.py -i WIGGLE.wig -o BEDGRAPH.bedgraph''')
    parser.add_argument('-i','--input', metavar='wigfile', required=True, help='a wig file to convert')
    parser.add_argument('-o','--output', metavar='bgfile', required=True, help='the output BedGraog file')

    args = parser.parse_args()
    return(args)

# track type=wiggle_0 name="100 Vert. Cons" description="100 vertebrates Basewise Conservation by PhyloP"
# #   output date: 2024-02-19 19:33:16 UTC
# #   chrom specified: chr6
# #   position specified: 42847238-42848238
# #   Subtrack merge, primary table = phyloP100wayAll (100 vertebrates Basewise Conservation by PhyloP)
# #   Subtrack merge operation: product of phyloP100wayAll and selected subtracks:
# variableStep chrom=chr6 span=1
# 42847239    -0.257307
# 42847240    -2.82533
# 42847241    -0.879858
# 42847242    -0.490764
# 42847243    -0.72422
# 42847244    -0.0238504
# 42847245    -0.101669
# 42847246    -0.646402
# 42847247    -0.179488
# 42847248    0.598701
# 42847249    0.598701
# 42847250    0.209606
# 42847251    -1.26895
# 42847252    0.131787
# 42847253    -0.101669


if __name__ == "__main__":
    '''Convert a wiggle format to a BedGraph format.'''

    args = getParams()

    sys.stderr.write( 'Begin parsing wiggle file: %s\n' % args.input )
    this_chr = 'NoChr'
    this_span = 1

    # Initialize writer
    writer = open(args.output, 'w')

    # Store Signal information from Signal file
    reader = open(args.input,'r')
    for line in reader:
        print(line)
        # skip comments
        if( line.strip().find('#')==0 ):
            continue
        elif( line.find('track')==0 ):
            continue
        # Find and set chr name
        elif( line.find('variableStep')==0 ):
            try:
                this_chr = re.search(" chrom=(.+) ", line).group(1)
            except AttributeError:
                print("Unable to parse chromosome info from line: " % line)
            try:
                this_span = int(re.search(" span=(.+) ", line).group(1))
            except AttributeError:
                this_span = 1
            continue
        # tokenize and write to output
        tokens = line.strip().split('\t')
        start = int(tokens[0]) - 1
        writer.write("%s\t%i\t%i\t%s\n" % (this_chr, start, start+this_span, tokens[1]))
    # Close file i/o
    reader.close()
    writer.close()

import sys
import argparse
import pandas as pd

def getParams():
    '''Parse parameters from the command line'''
    parser = argparse.ArgumentParser(description='''
This script will element-wise sum two CDT matrices. Checks for matching YORF and NAME

Example: python sum_each_CDT.py -1 FIRST.cdt -2 SECOND.cdt -o SUM.cdt''')
    parser.add_argument('-1','--file1', metavar='cdt_fn', required=True, help='the first CDT file to sum')
    parser.add_argument('-2','--file2', metavar='cdt_fn', required=True, help='the second CDT file to sum')
    parser.add_argument('-o','--output', metavar='cdt_fn', required=True, help='the summed CDT file')

    args = parser.parse_args()
    return(args)

if __name__ == "__main__":
    '''Calls max signal within CDT pileup and returns BED-formatted coordinates for the max signal positions.'''

    # Load reference BED coordinates
    args = getParams()

    mat1 = pd.read_csv(args.file1, sep="\t", header=0, index_col=[0,1])
    mat2 = pd.read_csv(args.file2, sep="\t", header=0, index_col=[0,1])

    mat1.add(mat2)
    mat1.to_csv(args.output, sep="\t")


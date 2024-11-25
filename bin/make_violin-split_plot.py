#!/usr/bin/env python
from os.path import splitext
import sys, re, argparse, random
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Python 3.6+
# relies on dict insertion order

# Check Matplotlib colors when building your config files: https://matplotlib.org/stable/gallery/color/named_colors.html

def getParams():
	'''Parse parameters from the command line'''
	parser = argparse.ArgumentParser(description='')

	parser.add_argument('-i','--input', metavar='two_col_file', dest='data_file', required=False, default=None, help='tab-delimited file made of three columns: first column y values to plot (must all be numeric values), second column is the grouping (which violin group along x-axis to contribute to), third column is "Proximal" or "Distal" category to split by (Proximal left, Distal right)')
	parser.add_argument('--width', metavar='width', required=False, type=int, default=8, help='width of figure')
	parser.add_argument('--height', metavar='height', required=False, type=int, default=4, help='height of figure')
	parser.add_argument('--title', metavar='title', required=False, default=None, help='title of figure')
	parser.add_argument('--xlabel', metavar='xlabel', required=False, default=None, help='x-axis label')
	parser.add_argument('--ylabel', metavar='ylabel', required=False, default=None, help='y-axis label')
	parser.add_argument('--preset1', action='store_true', help='use proximal/distal presets for nucleosome intervals (Fig 4d)')
	parser.add_argument('--preset2', action='store_true', help='use proximal/distal presets for half-nucleosome intervals (Fig 4d)')
	# parset.add_argument('--cut', a)
	parser.add_argument('-o','--output', metavar='output_svg', dest='output_svg', required=False, default=None, help='name of SVG filepath to save figure to (if none provided, figure pops up in new window)')

	args = parser.parse_args()

	return(args)


preset1_order = [
	"H2AZ-H2B",
	"H3K4me3-H3",
	"H3K9ac-H3",
	"H3K27ac-H3"
]

preset2_order = [
	"H3K4me3-H3",
	"H3K9ac-H3",
	"H3K27ac-H3"
]


# Example Data:
# 0.550716	category1
# 0.493109	category3
# 0.401034	category2
# 0.498233	category1
# 0.579172	category3
# 0.658480	category1
# 0.386018	category3
# 0.464670	category1
# 0.481569	category2

if __name__ == "__main__":
	'''Plot violin plot'''
	args = getParams()

	# Load data
	data = pd.read_csv(args.data_file, sep="\t", header=None)

	# Relabel columns (third label if split)
	data.rename(columns={data.columns[0]: 'y', data.columns[1]: 'category'}, inplace=True)
	data.rename(columns={data.columns[2]: 'split-category'}, inplace=True)
	print(data)

	# Initialize plot base with user-specified size
	fig, ax= plt.subplots()
	plt.figure(figsize=(args.width, args.height))

	# Hardcode split options
	custom_hue_order = ["Proximal", "Distal"]
	custom_pal = {"Proximal":"#808080", "Distal":"#E6E6E6"}

	# Hardcode order presets
	custom_order = sorted(list(data['category'].unique()))
	if (args.preset1):
		custom_order = preset1_order
	elif (args.preset2):
		custom_order = preset2_order

	# Scatter Plot or other type of plot here!!!!!
	# could swap out for 'swarmplot' and google seaborn library for others
	ax = sns.violinplot(x="category", y="y", hue="split-category", split=True, data=data, palette=custom_pal, order=custom_order, hue_order=custom_hue_order, inner="quart", linecolor='black', linewidth=0.5)

	# Format x-axis tickmarks
	ax.tick_params(axis='x', rotation=0)

	# Add title
	if (args.title!=None):
		ax.set_title(args.title)
	# Label x-axis
	if (args.xlabel!=None):
		plt.xlabel(args.xlabel)
	# Label y-axis
	if (args.ylabel!=None):
		plt.ylabel(args.ylabel)

	# Output...
	if(args.output_svg==None):
		plt.show()
	elif(splitext(args.output_svg)[-1]!=".svg"):
		sys.stderr.write("Please use SVG file extension to save output!\n")
		plt.savefig(args.output_svg, transparent=True)
		plt.show()
	else:
		plt.savefig(args.output_svg, transparent=True)
		sys.stderr.write("SVG written to %s\n" % args.output_svg)

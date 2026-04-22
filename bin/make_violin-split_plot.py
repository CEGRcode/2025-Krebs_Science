#!/usr/bin/env python
from os.path import splitext
import sys, argparse
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Use text, not shapes of text in SVG
plt.rcParams['svg.fonttype'] = 'none'

def getParams():
    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-i','--input', dest='data_file', required=False, default=None,
                        help='tab-delimited file: col1=y, col2=category, col3=split group')
    parser.add_argument('--width', type=int, default=8)
    parser.add_argument('--height', type=int, default=4)
    parser.add_argument('--title', default=None)
    parser.add_argument('--xlabel', default=None)
    parser.add_argument('--ylabel', default=None)
    parser.add_argument('--preset1', action='store_true')
    parser.add_argument('--preset2', action='store_true')
    parser.add_argument('-o','--output', dest='output_svg', default=None)

    return parser.parse_args()


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


if __name__ == "__main__":
    args = getParams()

    # Load data
    data = pd.read_csv(args.data_file, sep="\t", header=None)

    # Rename columns
    data.rename(columns={
        data.columns[0]: 'y',
        data.columns[1]: 'category',
        data.columns[2]: 'split-category'
    }, inplace=True)

    print(data)

    # Initialize plot
    fig, ax = plt.subplots(figsize=(args.width, args.height))

    # ---- KEY PART: force left = blue, right = red ----
    # First in hue_order = LEFT, second = RIGHT
    custom_hue_order = ["Proximal", "Distal"]
    custom_pal = {
        "Proximal": "#0000ff",   # blue (LEFT)
        "Distal": "#ff0000"      # red (RIGHT)
    }

    # Category order
    custom_order = sorted(list(data['category'].unique()))
    if args.preset1:
        custom_order = preset1_order
    elif args.preset2:
        custom_order = preset2_order

    # Plot
    ax = sns.violinplot(
        x="category",
        y="y",
        hue="split-category",
        split=True,
        data=data,
        palette=custom_pal,
        order=custom_order,
        hue_order=custom_hue_order,
        inner="quart",
        linecolor='black',
        linewidth=0.5
    )

    # Formatting
    ax.tick_params(axis='x', rotation=0)

    if args.title:
        ax.set_title(args.title)
    if args.xlabel:
        ax.set_xlabel(args.xlabel)
    if args.ylabel:
        ax.set_ylabel(args.ylabel)

    # Output
    if args.output_svg is None:
        plt.show()
    elif splitext(args.output_svg)[-1] != ".svg":
        sys.stderr.write("Please use SVG file extension!\n")
        plt.savefig(args.output_svg, transparent=True)
        plt.show()
    else:
        plt.savefig(args.output_svg, transparent=True)
        sys.stderr.write(f"SVG written to {args.output_svg}\n")
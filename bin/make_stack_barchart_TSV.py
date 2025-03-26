#!/usr/bin/env python
import argparse
import sys
from os.path import isfile, join, splitext
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Use text, not shapes of text in SVG
plt.rcParams['svg.fonttype'] = 'none'

# Python 3.6+
# relies on dict insertion order

# Check Matplotlib colors when building your config files: https://matplotlib.org/stable/gallery/color/named_colors.html

def getParams():
    '''Parse parameters from the command line'''
    parser = argparse.ArgumentParser(description='Generate a stacked bar chart from input data.')

    parser.add_argument('-i','--input', metavar='tsv_file', required=True, default=None, help='tab-delimited file where each column label is the x-axis category and each row is a different group (color)')
    parser.add_argument('-o','--output', metavar='output_svg', required=False, default=None, help='name of SVG filepath to save figure to (if none provided, figure pops up in new window)')

    parser.add_argument('--palette', metavar='seaborn_palette_name', required=False, default=None, help='name of Seaborn palette to use')
    parser.add_argument('--title', metavar='string', required=False, default=None, help='title to add to the plot')
    parser.add_argument('--entropy', action='store_true', help='add entropy line')

    parser.add_argument('--height', metavar='string', required=False, default=6, help='figure height')
    parser.add_argument('--width', metavar='string', required=False, default=10, help='figure width')

    args = parser.parse_args()
    return(args)

# Custom base colors
def set_base_colors():
    """Return specific colors for bases"""
    return {'A': 'red', 'C': 'blue', 'G': '#FFC834', 'T': '#42E500'}

def calculate_entropy(df):
    """Calculate entropy for each position."""
    # Turn counts into frequencies
    freq = df.div(df.sum())

    # Take log2 of frequencies
    log2f = -1 * np.log2(freq)

    # f * np.log2(f)
    flog2f = freq * log2f

    # Calculate the entropy for each position
    entropy = flog2f.sum()

    return(entropy)

if __name__ == "__main__":
    '''Plot stacked bar chart'''
    args = getParams()

    # Parse input data to a DataFrame
    df = pd.read_csv(args.input, sep='\t', index_col=0, header=0)
    print(df)

    # Get base-specific colors
    base_colors = set_base_colors()

    # Plotting
    fig, ax1 = plt.subplots(figsize=(args.width, args.height))

    # Initialize the bottom array for stacking
    bottom = np.zeros(len(df.columns))

    # Plot each category as a stacked bar
    for i, category in enumerate(df.index):
        # Set base color for each base in the category
        base_color = base_colors.get(category, 'gray')  # default to gray if not found

        # Plot bar from bottom and update
        ax1.bar(df.columns, df.loc[category, :], color=base_color, label=category, bottom=bottom)
        bottom += df.loc[category, :]  # Update the bottom for the next category

    # Add labels and title for the primary y-axis
    ax1.set_xlabel('X_Values')
    ax1.set_ylabel('Count (Primary Y-Axis)')
    if args.title is not None:
        ax1.set_title(args.title)
    ax1.set_ylim(bottom=0, top=max(df.sum()))

    # Store legend info
    lines, labels = ax1.get_legend_handles_labels()

    if (args.entropy):
        # Create a secondary y-axis for the line plot
        ax2 = ax1.twinx()

        # Calculate and plot entropy
        ax2.plot(df.columns, calculate_entropy(df).T, color='orchid', marker='o', label='Entropy', linewidth=2)

        # Add label for the secondary y-axis
        ax2.set_ylabel('Entropy', color='orchid')
        ax2.tick_params(axis='y', labelcolor='orchid')  # Set secondary y-axis tick color to purple
        ax2.set_ylim(bottom=0, top=2)

        # Add secondary legend info
        lines2, labels2 = ax2.get_legend_handles_labels()
        lines = lines + lines2
        labels = labels + labels2

    # Combine legends from both axes
    ax1.legend(lines, labels, title='Category', bbox_to_anchor=(1.1, 1), loc='upper left')

    # Show the plot
    plt.tight_layout()

    # Output...
    if(args.output==None):
        plt.show()
    else:
        try:
            plt.savefig(args.output, transparent=True)
            sys.stderr.write("Image written to %s\n" % args.output)
        except:
            sys.stderr.write("Please use SVG/PNG file extension to save output!\n")
            plt.show()

#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# -----------------------
# INPUT
# -----------------------
infile = sys.argv[1]
outfile = sys.argv[2] if len(sys.argv) > 2 else "violin.svg"

# -----------------------
# LOAD DATA
# -----------------------
data = pd.read_csv(infile, sep="\t", header=None, names=["value", "group"])

# remove zeros (log safety)
data = data[data["value"] > 0]

# log10 transform
data["value"] = np.log10(data["value"])

# -----------------------
# ORDER + COLORS
# -----------------------
order = [
    "MNase-Overlap",
    "MNase-NoOverlap",
    "BNase-Overlap",
    "BNase-NoOverlap"
]

palette = {
    "MNase-Overlap": "#0000FF",
    "MNase-NoOverlap": "#FF0000",
    "BNase-Overlap": "#0000FF",
    "BNase-NoOverlap": "#FF0000"
}

# -----------------------
# PLOT
# -----------------------
plt.figure(figsize=(5, 4))

ax = sns.violinplot(
    data=data,
    x="group",
    y="value",
    order=order,
    palette=palette,
    inner="quartile",
    cut=0,
    linewidth=0.8
)

# -----------------------
# FORMATTING
# -----------------------
ax.set_xlabel("")
ax.set_ylabel("Tag Occupancy (log10)")
ax.set_ylim(0, 4)

ax.tick_params(axis='x', rotation=45)

plt.tight_layout()

# -----------------------
# OUTPUT
# -----------------------
plt.savefig(outfile, dpi=300, transparent=True)
print(f"Saved: {outfile}")

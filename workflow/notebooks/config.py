# Config

import pandas as pd
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------
# VARIABLES
# ------------------------------------------------------------------------

# Pandas
# pd.set_option("display.max_rows", None, "display.max_columns", None)
pd.set_option("display.max_rows", 10, "display.max_columns", None)

# Branch Support Thresholds (from IQTREE docs)
ALRT_THRESH = 80
UFBOOT_THRESH = 95

# Significant digits for writing newick files
BRANCH_LEN_SIG_DIG = 12

# Data parsing
NO_DATA_CHAR = "NA"

# How to color branch supports
LOW_COL = "black"
HIGH_COL = "red"
TERM_COL = "grey"
THRESH_COL = "blue"

# Continuous data color palette
CONT_COLOR_PAL = "rainbow"

# Nextstrain / augur / auspice
JSON_INDENT = 2

# Plotting Graphics
figsize = (6.4, 4.8)
figsize_alt = (9.6, 4.8)
dpi = 400

SM_FONT = 5
MED_FONT = 8
LG_FONT = 10

plt.rc("font", size=SM_FONT)  # controls default text sizes
plt.rc("figure", titlesize=LG_FONT)  # fontsize of the figure title
# plt.rc('axes', labelsize=MED_FONT)    # fontsize of the x and y labels
plt.rc("lines", linewidth=0.5)

FMT = "svg"

# ------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------

# This code is from the biopython Phylo module


def get_x_positions(tree):
    """Create a mapping of each clade to its horizontal position.
    Dict of {clade: x-coord}
    """
    depths = tree.depths()
    # If there are no branch lengths, assume unit branch lengths
    if not max(depths.values()):
        depths = tree.depths(unit_branch_lengths=True)
    return depths


# This code is from the biopython Phylo module
def get_y_positions(tree):
    """Create a mapping of each clade to its vertical position.
    Dict of {clade: y-coord}.
    Coordinates are negative, and integers for tips.
    """
    maxheight = tree.count_terminals()
    # Rows are defined by the tips
    heights = {
        tip: maxheight - i for i, tip in enumerate(reversed(tree.get_terminals()))
    }

    # Internal nodes: place at midpoint of children
    def calc_row(clade):
        for subclade in clade:
            if subclade not in heights:
                calc_row(subclade)
        # Closure over heights
        heights[clade] = (heights[clade.clades[0]] + heights[clade.clades[-1]]) / 2.0

    if tree.root.clades:
        calc_row(tree.root)
    return heights

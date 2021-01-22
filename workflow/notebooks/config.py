# Config

import pandas as pd


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

# Nextstrain / augur / auspice
JSON_INDENT = 2


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

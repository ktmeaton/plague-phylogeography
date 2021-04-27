#!/usr/bin/env python3
import sys
from Bio import Phylo

tree_input = sys.argv[1]
tree_output = sys.argv[2]

tree = Phylo.read(tree_input, "newick")

# Get minimum branch length for output precision
min_branch_length = 1e10

for c in tree.find_clades():
    # Remove comments
    c.comment = None
    if c.branch_length and c.branch_length < min_branch_length:
        min_branch_length = c.branch_length

# This is a terrible way to get sig digits
sig_dig = 5
# If we're using scientific notation, it's really small
if "e-" in str(min_branch_length):
    sig_dig = int(str(min_branch_length).split("e-")[-1]) + 2

Phylo.write(tree, tree_output, "nexus", format_branch_length="%1.{}f".format(sig_dig))

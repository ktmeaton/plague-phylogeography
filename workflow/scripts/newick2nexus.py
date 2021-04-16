#!/usr/bin/env python3
import sys
from Bio import Phylo

tree_input = sys.argv[1]
tree_output = sys.argv[2]

tree = Phylo.read(tree_input, "newick")
Phylo.write(tree, tree_output, "nexus")

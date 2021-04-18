#!/usr/bin/env python3

# Usage
# ./multi2bi.py \
#   --tree beast.2-MED.nex \
#   --out beast.2-MED.bifurc.nex

import click
from io import StringIO
from Bio import Phylo
import os
from ete3 import Tree

# -----------------------------------------------------------------------------
# Command Line Arguments


@click.command()
@click.help_option("--help", "-h")
@click.option(
    "-t",
    "--tree",
    help="Input multifurcating tree.",
    type=click.Path(exists=True, dir_okay=False, allow_dash=True),
    required=True,
)
@click.option(
    "-o",
    "--out",
    help="Output bifurcating tree.",
    type=click.Path(dir_okay=False, allow_dash=True),
    required=True,
)
def main(
    tree: str, out: str,
):
    """This script converts a multifurcating tree to a bifurcating tree."""

    tree_path = tree
    out_path = out

    # -------------------------------------------------------------------------
    # Identify tree type
    tree_ext = os.path.splitext(tree_path)[1]
    if tree_ext == ".nex" or tree_ext == ".nexus":
        tree_type = "nexus"
    elif tree_ext == ".nwk" or tree_ext == ".newick":
        tree_type = "newick"
    else:
        print("Unrecognized tree format:", tree_ext)
        quit()

    # Import tree
    tree = Phylo.read(tree_path, tree_type)

    # Strip comments and confidence
    for c in tree.find_clades():
        c.comment = None

    # Conver to bifurcating
    tree_multi = Tree(tree.format("newick"))
    tree_multi.standardize()

    tree_bi = Phylo.read(StringIO(tree_multi.write(format=0)), "newick")

    for c in tree_bi.find_clades():
        c.confidence = None

    Phylo.write(tree_bi, out_path, "nexus")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3

# Usage
# ./filter_alignment.py \
#   --tree beast.2-MED.nex \
#   --aln snippy-multi.snps.aln \
#   --out beast.2-MED.fasta

import click
from Bio import Phylo, AlignIO, Align, SeqIO
import os
import copy

# -----------------------------------------------------------------------------
# Command Line Arguments


@click.command()
@click.help_option("--help", "-h")
@click.option(
    "-t",
    "--tree",
    help="Input tree defining taxon set.",
    type=click.Path(exists=True, dir_okay=False, allow_dash=True),
    required=True,
)
@click.option(
    "-a",
    "--aln",
    help="Input alignment.",
    type=click.Path(exists=True, dir_okay=False, allow_dash=True),
    required=True,
)
@click.option(
    "--prune-tips",
    help="File of tips to prune, one tip per line",
    type=click.Path(dir_okay=False, allow_dash=True),
    required=False,
)
@click.option(
    "--prune-tree",
    help="File for output pruned tree.",
    type=click.Path(dir_okay=False, allow_dash=True),
    required=False,
)
@click.option(
    "-o",
    "--out",
    help="Output alignment that is filtered.",
    type=click.Path(dir_okay=False, allow_dash=True),
    required=True,
)
def main(
    tree: str, aln: str, out: str, prune_tips: str, prune_tree: str,
):
    """This script filters an alignment based on taxa in a tree"""

    tree_path = tree
    aln_path = aln
    out_path = out
    prune_tips_path = prune_tips
    prune_tree_path = prune_tree

    # -------------------------------------------------------------------------
    # Identify tree type
    tree_ext = os.path.splitext(tree_path)[1]
    if tree_ext == ".nex" or tree_ext == ".nexus":
        tree_type = "nexus"
    elif tree_ext == ".nwk" or tree_ext == ".newick" or tree_ext == ".treefile":
        tree_type = "newick"
    else:
        print("Unrecognized tree format:", tree_ext)
        quit()

    # Import tree
    tree = Phylo.read(tree_path, tree_type)

    # Parse taxa to prune
    prune_taxa = []
    if prune_tips:
        with open(prune_tips_path) as infile:
            prune_taxa = infile.read().strip().split("\n")

    # Store the taxa in the tree
    tree_taxa = []
    tree_edit = copy.deepcopy(tree)
    for c in tree.find_clades():
        # Prune taxa
        if c.name in prune_taxa:
            tree_edit.prune(c.name)

        if c.name in prune_taxa and prune_tips:
            continue
        else:
            tree_taxa.append(c.name)

    # Write the pruned tree
    if prune_tree:
        Phylo.write(tree_edit, prune_tree_path, "newick")

    # Import alignment
    aln = AlignIO.read(aln_path, "fasta")
    filter_seq = [rec for rec in aln if rec.id in tree_taxa]
    filter_aln = Align.MultipleSeqAlignment(filter_seq)

    with open(out_path, "w") as outfile:
        count = SeqIO.write(filter_aln, outfile, "fasta")
        print("\n", count, "alignments written.")


if __name__ == "__main__":
    main()

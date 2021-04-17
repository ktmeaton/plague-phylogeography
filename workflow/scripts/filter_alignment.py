#!/usr/bin/env python3

# Usage
# ./filter_alignment.py \
#   --tree beast.2-MED.nex \
#   --aln snippy-multi.snps.aln \
#   --out beast.2-MED.fasta

import click
from Bio import Phylo, AlignIO, Align, SeqIO
import os

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
    "-o",
    "--out",
    help="Output alignment that is filtered.",
    type=click.Path(dir_okay=False, allow_dash=True),
    required=True,
)
def main(
    tree: str, aln: str, out: str,
):
    """This script filters an alignment based on taxa in a tree"""

    tree_path = tree
    aln_path = aln
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

    # Store the taxa in the tree
    tree_taxa = []
    for c in tree.find_clades():
        tree_taxa.append(c.name)

    # Import alignment
    aln = AlignIO.read(aln_path, "fasta")
    filter_seq = [rec for rec in aln if rec.id in tree_taxa]
    filter_aln = Align.MultipleSeqAlignment(filter_seq)

    with open(out_path, "w") as outfile:
        count = SeqIO.write(filter_aln, outfile, "fasta")
        print("\n", count, "alignments written.")


if __name__ == "__main__":
    main()

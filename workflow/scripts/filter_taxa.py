#!/usr/bin/env python3

# Usage
# ./filter_alignment.py \
#   --tree beast.2-MED.nex \
#   --aln snippy-multi.snps.aln \
#   --out beast.2-MED.fasta

# USE CASE #1
# Filter inputs based on taxa in a tree.
#   inputs: aln, tsv
# USE CASE #2
# Filter inputs based on taxa in text files
#    inputs: aln, tsv, tree

import click
from Bio import Phylo, AlignIO, Align, SeqIO
import os
import copy
import pandas as pd

# -----------------------------------------------------------------------------
# Command Line Arguments


@click.command()
@click.help_option("--help", "-h")
@click.option(
    "-m",
    "--metadata",
    help="Input metadata.",
    type=click.Path(dir_okay=False, allow_dash=True),
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
    "--outdir",
    help="Output directory.",
    type=click.Path(dir_okay=True, allow_dash=True),
    required=True,
)
@click.option(
    "-t",
    "--tree",
    help="Input tree defining taxon set.",
    type=click.Path(exists=True, dir_okay=False, allow_dash=True),
    required=False,
)
@click.option(
    "--prune-tips",
    help="File of tips to prune, one tip per line",
    type=click.Path(dir_okay=False, allow_dash=True),
    required=False,
)
@click.option(
    "--keep-tips",
    help="File of tips to keep, one tip per line",
    type=click.Path(dir_okay=False, allow_dash=True),
    required=False,
)
def main(
    tree: str, aln: str, prune_tips: str, keep_tips: str, outdir: str, metadata: str,
):
    """This script filters an alignment based on taxa in a tree"""

    tree_path = tree
    aln_path = aln
    prune_tips_path = prune_tips
    keep_tips_path = keep_tips
    metadata_path = metadata

    # -------------------------------------------------------------------------
    # Identify Mode
    if prune_tips_path and keep_tips_path:
        print("ERROR: Only one of --prune-tips or --keep-tips can be selected.")
        quit()

    # Mandatory
    aln_basename = os.path.basename(aln_path)
    out_path_aln = os.path.join(outdir, aln_basename)

    metadata_basename = os.path.basename(metadata_path)
    out_path_metadata = os.path.join(outdir, metadata_basename)

    # Using an input tree
    if tree_path:
        tree_basename = os.path.basename(tree_path)
        tree_ext = os.path.splitext(tree_basename)[1]
        if tree_ext == ".nwk" or ".newick" or ".treefile":
            tree_type = "newick"
        elif tree_ext == ".nex" or ".nexus":
            tree_type = "nexus"
        if prune_tips_path or keep_tips_path:
            filter_tree_path = os.path.join(outdir, tree_basename)

    # -------------------------------------------------------------------------
    # Load Metadata
    metadata_df = pd.read_csv(metadata_path, sep="\t")
    metadata_df.fillna("NA", inplace=True)
    metadata_df.set_index(metadata_df.columns[0], inplace=True)

    # Import alignment
    aln = AlignIO.read(aln_path, "fasta")

    # Parse filtering taxa files
    keep_taxa = [rec.id for rec in aln]

    if prune_tips:
        with open(prune_tips_path) as infile:
            prune_taxa = infile.read().strip().split("\n")
            for taxa in prune_taxa:
                keep_taxa.remove(taxa)

    elif keep_tips:
        with open(keep_tips_path) as infile:
            keep_taxa = infile.read().strip().split("\n")

    # -------------------------------------------------------------------------
    # Identify tree type
    if tree_path:
        tree_ext = os.path.splitext(tree_path)[1]
        if tree_ext == ".nex" or tree_ext == ".nexus":
            tree_type = "nexus"
        elif tree_ext == ".nwk" or tree_ext == ".newick" or tree_ext == ".treefile":
            tree_type = "newick"
        else:
            print("Unrecognized tree format:", tree_ext)
            quit()

        # Import tree
        trees = Phylo.parse(tree_path, tree_type)

        # Parse out the tree
        for tree in trees:
            clades = [c for c in tree.find_clades()]
            if len(clades) > 1:
                break

        # Store the taxa in the tree
        tree_taxa = []
        tree_edit = copy.deepcopy(tree)

        # Get minimum branch length for output precision
        min_branch_length = 1e10

        for c in tree.get_terminals():
            # Prune taxa
            if c.name not in keep_taxa:
                print("Pruning taxa", c.name, "from the tree.")
                tree_edit.prune(c.name)

            else:
                tree_taxa.append(c.name)
                if c.branch_length and c.branch_length < min_branch_length:
                    min_branch_length = c.branch_length

        # This is a terrible way to get sig digits
        sig_dig = 5
        # If we're using scientific notation, it's really small
        if "e-" in str(min_branch_length):
            sig_dig = int(str(min_branch_length).split("e-")[-1]) + 2

    # -----------------------------------------------------
    # FILTER

    # Alignment
    filter_seq = [rec for rec in aln if rec.id in keep_taxa]
    filter_aln = Align.MultipleSeqAlignment(filter_seq)
    with open(out_path_aln, "w") as outfile:
        count = SeqIO.write(filter_aln, outfile, "fasta")
        print("\n", count, "alignments written.")

    # Metadata
    filter_df = copy.deepcopy(metadata_df)
    for sample in metadata_df.index:
        if sample not in keep_taxa:
            filter_df.drop(sample, inplace=True)
    filter_df.to_csv(out_path_metadata, sep="\t")

    # Write the pruned tree
    if tree_path and (prune_tips or keep_tips):
        Phylo.write(
            tree_edit,
            filter_tree_path,
            tree_type,
            format_branch_length="%1.{}f".format(sig_dig),
        )


if __name__ == "__main__":
    main()

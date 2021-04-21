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
import subprocess
import pandas as pd

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
    "-o",
    "--outdir",
    help="Output directory.",
    type=click.Path(dir_okay=True, allow_dash=True),
    required=True,
)
@click.option(
    "-m",
    "--metadata",
    help="Input metadata.",
    type=click.Path(dir_okay=False, allow_dash=True),
    required=True,
)
def main(
    tree: str, aln: str, prune_tips: str, outdir: str, metadata: str,
):
    """This script filters an alignment based on taxa in a tree"""

    tree_path = tree
    aln_path = aln
    prune_tips_path = prune_tips
    metadata_path = metadata

    # -------------------------------------------------------------------------
    # Output file names
    aln_basename = os.path.basename(aln_path)
    aln_ext = os.path.splitext(aln_basename)[1]
    tree_basename = os.path.basename(tree_path)
    tree_prefix = os.path.splitext(tree_basename)[0]
    metadata_basename = os.path.basename(metadata_path)

    # Output alignment takes its name from the tree
    out_path_aln = os.path.join(outdir, tree_prefix + ".filter" + aln_ext)
    out_path_metadata = os.path.join(outdir, metadata_basename)

    # Path for pruned tree
    prune_tree_prefix = os.path.join(outdir, tree_prefix + ".filter")

    # -------------------------------------------------------------------------
    # Load Metadata
    metadata_df = pd.read_csv(metadata_path, sep="\t")
    metadata_df.fillna("NA", inplace=True)
    metadata_df.set_index(metadata_df.columns[0], inplace=True)

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

    # -----------------------------------------------------
    # Parse Tree
    # Import tree
    trees = Phylo.parse(tree_path, tree_type)

    # Parse out the tree
    for tree in trees:
        clades = [c for c in tree.find_clades()]
        if len(clades) > 1:
            break

    # Parse taxa to prune
    prune_taxa = []
    if prune_tips:
        with open(prune_tips_path) as infile:
            prune_taxa = infile.read().strip().split("\n")

    # Store the taxa in the tree
    tree_taxa = []
    tree_edit = copy.deepcopy(tree)

    # Get minimum branch length for output precision
    min_branch_length = 1e10

    for c in tree.find_clades():
        # Prune taxa
        if c.name in prune_taxa:
            print("Pruning taxa", c.name, "from the tree.")
            tree_edit.prune(c.name)

        if c.name in prune_taxa and prune_tips:
            continue
        else:
            tree_taxa.append(c.name)
            if c.branch_length and c.branch_length < min_branch_length:
                min_branch_length = c.branch_length

    # This is a terrible way to get sig digits
    sig_dig = 5
    # If we're using scientific notation, it's really small
    if "e-" in str(min_branch_length):
        sig_dig = int(str(min_branch_length).split("e-")[-1]) + 2

    # Write the pruned tree
    if prune_tips:
        Phylo.write(
            tree_edit,
            prune_tree_prefix + ".nwk",
            "newick",
            format_branch_length="%1.{}f".format(sig_dig),
        )
        Phylo.write(
            tree_edit,
            prune_tree_prefix + ".nex",
            "nexus",
            format_branch_length="%1.{}f".format(sig_dig),
        )

    # -----------------------------------------------------
    # Filter Alignment

    # Import alignment
    aln = AlignIO.read(aln_path, "fasta")
    filter_seq = [rec for rec in aln if rec.id in tree_taxa]
    filter_aln = Align.MultipleSeqAlignment(filter_seq)

    with open(out_path_aln + ".tmp", "w") as outfile:
        count = SeqIO.write(filter_aln, outfile, "fasta")
        print("\n", count, "alignments written.")

    # Remove constant sites
    subprocess.run(["snp-sites", "-m", "-o", out_path_aln, out_path_aln + ".tmp"])

    # Remove temporary file
    os.remove(out_path_aln + ".tmp")

    # -----------------------------------------------------
    # Filter Metadata
    filter_df = copy.deepcopy(metadata_df)
    for sample in metadata_df.index:
        if sample not in tree_taxa:
            filter_df.drop(sample, inplace=True)

    filter_df.to_csv(out_path_metadata, sep="\t")


if __name__ == "__main__":
    main()

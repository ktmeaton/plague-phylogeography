#!/usr/bin/env python3

# Usage
# ./beast_nexus.py \
#   -m metadata.tsv

import click
from Bio import AlignIO  # Multiple Alignment parsing
from Bio import Nexus
from Bio import Phylo
import pandas as pd

"""
from Bio import SeqIO  # Nucleotide parsing
from Bio.Seq import Seq  # Nucleotide parsing
from Bio.SeqRecord import SeqRecord  # Sequencing create for writing
from Bio.Alphabet import generic_dna
from Bio import Nexus
"""

CALIBRATE_BASE_STR = "CALIBRATE {} = {}"

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
    help="Input fasta alignment.",
    type=click.Path(dir_okay=False, allow_dash=True),
    required=True,
)
@click.option(
    "--nwk",
    help="Input newick tree.",
    type=click.Path(dir_okay=False, allow_dash=True),
    required=True,
)
@click.option(
    "--nex",
    help="Output nexus.",
    type=click.Path(dir_okay=False, allow_dash=True),
    required=True,
)
def main(
    metadata, aln, nex, nwk,
):
    """Create a comprehensive nexus for BEAST2."""

    tree = Phylo.read(nwk, "newick")
    metadata_df = pd.read_csv(metadata, sep="\t", index_col="sample")

    # Step 1. Create Nexus with alignment
    alignment = AlignIO.read(open(aln), "fasta")
    n = Nexus.Nexus.Nexus()

    for rec in alignment:
        n.add_sequence(sequence=str(rec.seq), name=rec.id)
    n.write_nexus_data(filename=nex)

    # -----------------
    # Step 2. Add assumptions
    calibrations = []

    for c in tree.get_terminals():
        date_mean = metadata_df["date_bp_mean"][c.name]
        date_err = metadata_df["date_err"][c.name]

        prior = "fixed({})".format(date_mean)
        if date_err > 1:
            # By default, use uncertainty divided by 2 as std
            prior = "normal({},{})".format(date_mean, date_err / 2)

        calibrations.append(CALIBRATE_BASE_STR.format(c.name, prior))

    # Add the formatting char
    assumptions = "\t" + ",\n\t".join(calibrations) + ";"
    assumptions_block = (
        "begin ASSUMPTIONS;"
        + "\n\tOPTIONS SCALE = years;"
        + "\n\n{}\n\nend;".format(assumptions)
    )

    with open(nex, "a") as nex_file:
        nex_file.write("\n")
        nex_file.write(assumptions_block)

    # -----------------
    # Step 3. Add tree
    writer = Phylo.NewickIO.Writer(trees=[tree])
    nwk_str = ""
    for tree_str in writer.to_strings(format_branch_length="%1.10f"):
        nwk_str = tree_str

    trees_block = "begin Trees;\n\tTree tree1={}\nend;".format(nwk_str)

    with open(nex, "a") as nex_file:
        nex_file.write("\n\n")
        nex_file.write(trees_block)


if __name__ == "__main__":
    main()

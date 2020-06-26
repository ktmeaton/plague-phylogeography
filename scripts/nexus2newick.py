#!/usr/bin/env python3
"""
@author: Katherine Eaton

Convert Nexus tree to Newick Format.

./nexus2newick.py --nexus tree.nexus --newick tree.nwk
"""

# -----------------------------------------------------------------------#
#                         Modules and Packages                          #
# -----------------------------------------------------------------------#
import argparse  # Command-line argument parsing
import dendropy  # Read and write trees
import os
import sys

# This program should only be called from the command-line
if __name__ != "__main__":
    quit()

# -----------------------------------------------------------------------#
#                            Argument Parsing                           #
# -----------------------------------------------------------------------#

parser = argparse.ArgumentParser(
    description="Convert Nexus tree to Newick Format.", add_help=True
)

# Argument groups for the program

parser.add_argument(
    "--nexus",
    help="Path to the input nexus tree.",
    action="store",
    dest="nexusPath",
    required=True,
)

parser.add_argument(
    "--newick",
    help="Path to the output newick tree.",
    action="store",
    dest="newickPath",
    required=True,
)

# Retrieve user parameters
args = vars(parser.parse_args())

nexus_path = args["nexusPath"]
newick_path = args["newickPath"]

# ------------------------------------------------------------------------------#
#                                Error Catching                                #
# ------------------------------------------------------------------------------#
post_trees = dendropy.TreeList()

# Check if  NEXUS file exists
if not os.path.exists(nexus_path):
    print("An error occurred while trying to open", nexus_path)
    sys.exit(1)

post_trees.read(file=open(nexus_path, "r"), schema="nexus")

post_trees.write(path=newick_path, schema="newick")

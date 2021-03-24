#!/usr/bin/env python3
"""
@author: Katherine Eaton

Convert treetime clock output to auspice json.

./treetime_dates_json.py
"""

# -----------------------------------------------------------------------#
#                         Modules and Packages                          #
# -----------------------------------------------------------------------#
import argparse  # Command-line argument parsing
import os

import json
from Bio import Phylo

# This program should only be called from the command-line
if __name__ != "__main__":
    quit()

# -----------------------------------------------------------------------#
#                            Argument Parsing                           #
# -----------------------------------------------------------------------#

parser = argparse.ArgumentParser(
    description="Convert treetime clock output to auspice json..", add_help=True
)

# Argument groups for the program

parser.add_argument(
    "--time",
    help="Path to the time tree.",
    action="store",
    dest="timePath",
    required=True,
)

parser.add_argument(
    "--dates",
    help="Path to the dates tsv.",
    action="store",
    dest="datesPath",
    required=True,
)

parser.add_argument(
    "--delim",
    help="Dates file delimiter.",
    action="store",
    dest="datesDelim",
    required=False,
    default="\t",
)

parser.add_argument(
    "--json", help="Output json file.", action="store", dest="jsonPath", required=True
)

# Retrieve user parameters
args = vars(parser.parse_args())

time_path = args["timePath"]
dates_path = args["datesPath"]
dates_delim = args["datesDelim"]
json_path = args["jsonPath"]

filename, file_extension = os.path.splitext(time_path)
if file_extension == ".nexus" or file_extension == ".nex":
    format = "nexus"
elif file_extension == ".nwk" or file_extension == ".newick":
    format = "newick"
else:
    exit(1)

tree = Phylo.read(time_path, format)
dates_tsv = open(dates_path, "r")
output_json = open(json_path, "w")

# node   date    numeric date    lower bound     upper bound
NAME_INDEX = 0
DATE_INDEX = 1
NUMERIC_INDEX = 2
LOWER_INDEX = 3
UPPER_INDEX = 4

node_dict = {"nodes": {}}

# Tree parsing
for c in tree.find_clades():
    node_dict["nodes"][c.name] = {}
    node_dict["nodes"][c.name]["branch_length"] = c.branch_length
    node_dict["nodes"][c.name]["clock_length"] = c.branch_length

# Dates confidence parsing
read_line = dates_tsv.readline().strip()
# Skip the header lines
while read_line[0] == "#":
    read_line = dates_tsv.readline().strip()

while read_line:
    split_line = read_line.split(dates_delim)
    date_name = split_line[NAME_INDEX]
    date_date = split_line[DATE_INDEX]
    date_numeric = float(split_line[NUMERIC_INDEX])
    date_lower = float(split_line[LOWER_INDEX])
    date_upper = float(split_line[UPPER_INDEX])

    node_dict["nodes"][date_name]["date"] = date_date
    node_dict["nodes"][date_name]["numdate"] = date_numeric
    node_dict["nodes"][date_name]["num_date_confidence"] = [date_lower, date_upper]

    read_line = dates_tsv.readline().strip()

json_string = json.dumps(
    node_dict, default=lambda o: o.__dict__, sort_keys=True, indent=2
)
output_json.write(json_string)

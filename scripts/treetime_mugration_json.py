#!/usr/bin/env python3
"""
@author: Katherine Eaton

Convert treetime mugration output to auspice json.

./treetime_mugration_json.py
"""

# This program should only be called from the command-line
if __name__ != "__main__": quit()

#-----------------------------------------------------------------------#
#                         Modules and Packages                          #
#-----------------------------------------------------------------------#
import argparse # Command-line argument parsing
import os # Checking files and file extensions

import json # Write json files
from Bio import Phylo # Read and parse trees

#-----------------------------------------------------------------------#
#                            Argument Parsing                           #
#-----------------------------------------------------------------------#

parser = argparse.ArgumentParser(description="Convert treetime mugration output to auspice json.",
                                 add_help=True)

# Argument groups for the program

parser.add_argument('--tree',
                    help = 'Path to the tree file.',
                    action = 'store',
                    dest = 'treePath',
                    required = True)

parser.add_argument('--json',
                    help = 'Output json file.',
                    action = 'store',
                    dest = 'jsonPath',
                    required = True)

parser.add_argument('--trait',
                    help = 'Name of the mugration trait.',
                    action = 'store',
                    dest = 'mugTrait',
                    required = True)

parser.add_argument('--conf',
                    help = 'Path to the mugration confidence csv.',
                    action = 'store',
                    dest = 'mugConfPath',
                    required = True)

parser.add_argument('--model',
                    help = 'Path to the mugration model txt.',
                    action = 'store',
                    dest = 'mugModelPath',
                    required = False)

# Retrieve user parameters
args = vars(parser.parse_args())

tree_path = args['treePath']
json_path = args['jsonPath']
mug_trait = args['mugTrait']
mug_conf_path = args['mugConfPath']
mug_model_path = args['mugModelPath']

filename, file_extension = os.path.splitext(tree_path)

if file_extension == ".nexus" or file_extension == ".nex": format = "nexus"
elif file_extension == ".nwk" or file_extension == ".newick": format = "newick"
else: exit(1)

tree = Phylo.read(tree_path, format)
output_json = open(json_path, "w")
mug_conf = open(mug_conf_path, "r")

node_dict = {"nodes": {}, "models": {} }

# Tree parsing
for c in tree.find_clades():
    node_dict["nodes"][c.name] = {}

# Mug confidence parsing
mug_header_list = []
read_line = mug_conf.readline().strip()

# Parse out the header (possible mug trait values)
while read_line[0] == "#":
    split_line = read_line.split(",")
    for i in range(1,len(split_line)):
        mug_header_list.append(split_line[i].strip())
    read_line = mug_conf.readline().strip()

# Parse the rest of the mug confidence file
while read_line:
    split_line = read_line.split(",")
    mug_node_name = split_line[0].strip()
    mug_conf_dict = {}
    mug_conf_dict_name = mug_trait + "_" + "confidence"
    for i in range(1,len(split_line)):
        mug_conf_key = mug_header_list[i-1]
        mug_conf_val = float(split_line[i])
        mug_conf_dict[mug_conf_key] = mug_conf_val
    max_value_keys = [
            k for k in mug_conf_dict.keys()
            if mug_conf_dict[k] == max(mug_conf_dict.values())]
    mug_trait_assign = max_value_keys[0]
    # Assign the confidence dict to the node
    node_dict["nodes"][mug_node_name][mug_trait] = mug_trait_assign
    node_dict["nodes"][mug_node_name][mug_conf_dict_name] = mug_conf_dict
    read_line = mug_conf.readline().strip()

json_string = json.dumps(node_dict,
                         default=lambda o: o.__dict__,
                         sort_keys=True,
                         indent=2)
output_json.write(json_string)

## Cleanup
output_json.close()
mug_conf.close()

#!/usr/bin/env python3

# USAGE:
"""
tree2network.py \
  --tree ../../../plague-phylogeography-projects/main/mugration/all/chromosome_filter5/mugration_model_timetree.nwk \
  --tsv ../../../plague-phylogeography-projects/main/mugration/all/chromosome_filter5/mugration_model.tsv \
  --output test.js
"""

from Bio import Phylo
import pandas as pd
import argparse

# ------------------------------------------------------------------------------#
# Argument Parsing


parser = argparse.ArgumentParser(
    description=("Convert a phylogeny and metadata to javascript network."),
    add_help=True,
)

# Argument groups for the program

parser.add_argument(
    "--tree", help="Newick Phylogeny", action="store", dest="treePath", required=True,
)

parser.add_argument(
    "--tsv", help="Metadata TSV", action="store", dest="tsvPath", required=True,
)

parser.add_argument(
    "--output",
    help="Ouput Javascript Network.",
    action="store",
    dest="outPath",
    required=True,
)


# Retrieve user parameters
args = vars(parser.parse_args())
tree_path = args["treePath"]
tsv_path = args["tsvPath"]
out_path = args["outPath"]

# ------------------------------------------------------------------------------#
# Functions


def get_parent(tree, child_clade):
    node_path = [tree.root] + tree.get_path(child_clade)
    try:
        return node_path[-2]
    except IndexError:
        return None


def tree_network(tree):
    """
    Return a list of node connections.
    """
    network = []
    for c in tree.find_clades(order="postorder"):
        parent = get_parent(tree, c)
        if not parent:
            continue
        connection = [parent, c]
        if connection not in network:
            network.append(connection)

    # Put root back at start
    network.reverse()
    return network


GEO = "Country"

# ------------------------------------------------------------------------------#
# Nodes

# Import metadata as dataframe
tree_df = pd.read_csv(tsv_path, sep="\t")
tree_df.fillna("NA", inplace=True)
tree_df.set_index("Name", inplace=True)

nodes = {}
i = 0
for rec in tree_df.iterrows():
    node_name = rec[0]
    if "NODE" in node_name:
        continue
    country = rec[1]["Country"]

    if country not in nodes:
        # continent = rec[1]["Continent"]
        continent = "NA"
        lat = tree_df["Mugration_Country_Lat"][node_name]
        lon = tree_df["Mugration_Country_Lon"][node_name]

        nodes[country] = {
            "size": 0,
            "id": i,
            "country": country,
            "continent": continent,
            "lat": lat,
            "lon": lon,
        }
        i += 1
    nodes[country]["size"] += 1

# ------------------------------------------------------------------------------#
# Links

# Import Tree and convert to connections
tree = Phylo.read(tree_path, "newick")
phylo_network = tree_network(tree)

# Get network links
links = {}
for connection in phylo_network:
    # Parse the source node
    source = connection[0].name
    source_country = tree_df["Mugration_Country"][source]
    source_id = nodes[source_country]["id"]
    target = connection[1].name
    target_country = tree_df["Mugration_Country"][target]
    target_id = nodes[target_country]["id"]
    connection_name = "{}_{}".format(source_country, target_country)
    # Skip if it's staying in place
    if source_id == target_id:
        continue
    if connection_name not in links:
        links[connection_name] = {"source": source_id, "target": target_id, "size": 1}
    else:
        links[connection_name]["size"] += 1

# Write Javascript out
with open(out_path, "w") as outfile:
    outfile.write("var network = {" + "\n\tnodes:[" + "\n")
    i = 0
    for country in nodes:
        outfile.write(
            "\t\t{"
            + "id:{}, country: '{}', continent: '{}', size:{}, lat:{}, lon: {}".format(
                nodes[country]["id"],
                nodes[country]["country"],
                nodes[country]["continent"],
                nodes[country]["size"],
                nodes[country]["lat"],
                nodes[country]["lon"],
            )
            + "}"
        )
        if i < len(nodes) - 1:
            outfile.write(",")
        outfile.write("\n")
        i += 1

    outfile.write("\t]," "\n\tlinks:[" + "\n")

    i = 0
    for connection_name in links:
        outfile.write(
            "\t\t{"
            + "source:{}, target:{}, size:{}".format(
                links[connection_name]["source"],
                links[connection_name]["target"],
                links[connection_name]["size"],
            )
            + "}"
        )
        if i < len(links) - 1:
            outfile.write(",")
        outfile.write("\n")
        i += 1

    outfile.write("\t]" + "\n}" + "\n")

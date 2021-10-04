# -------------------------------------------------------------------------
# Dependencies
# mamba create -n newick2nextstrain -c conda-forge -c bioconda \
#   click pandas numpy biopython augur

import click
import os
from Bio import Phylo
import pandas as pd
import numpy as np
from augur import utils, export_v2
import json
import time

# -------------------------------------------------------------------------
# CONSTANTS

JSON_INDENT = 2
NO_DATA_CHAR = "NA"

# -------------------------------------------------------------------------
# Functions


def branch_attributes(tree_dict, sub_dict, df, label_col):
    """
    Add branch attributes to an auspice tree.
    """
    root = sub_dict
    if "children" not in root:
        return root
    else:
        for child in root["children"]:
            node = branch_attributes(tree_dict, child, df, label_col)
            node_type = df["node_type"][node["name"]]
            if node_type != "internal":
                continue
            branch_labels = {}
            for col in label_col:
                col_pretty = col.replace("_", " ").title()
                # clade needs to stay lowercase
                if col == "clade":
                    col_pretty = "clade"
                branch_labels[col_pretty] = df[col][node["name"]]

            node["branch_attrs"]["labels"] = branch_labels

        branch_labels = {}
        for col in label_col:
            col_pretty = col.replace("_", " ").title()
            # clade needs to stay lowercase
            if col == "clade":
                col_pretty = "clade"
            branch_labels[col_pretty] = df[col][node["name"]]
        root["branch_attrs"]["labels"] = branch_labels
        return root


def augur_export(
    tree_path=None,
    aln_path=None,
    tree=None,
    tree_df=None,
    color_keyword_exclude=None,
    type_convert=None,
):
    """
    Export the Augur JSON of tree node data.
    """
    augur_data = {}
    augur_data["alignment"] = aln_path
    augur_data["input_tree"] = tree_path
    augur_data["nodes"] = {}

    # Iterate through all nodes in the tree
    for c in tree.find_clades():
        # Add the node to the dictionary
        augur_data["nodes"][c.name] = {}

        # Iterate through all attributes in the dataframe
        for attr in tree_df.columns:
            attr_val = tree_df[attr][c.name]

            # Check if this attribute should be excluded
            exclude = False
            for keyword in color_keyword_exclude:
                if keyword in attr.lower():
                    exclude = True

            if not exclude:
                # Fix numpy attr to float/int
                if type(attr_val) == np.float64:
                    attr_val = float(attr_val)
                elif type(attr_val) == np.int64:
                    attr_val = int(attr_val)

                # Perform requested type conversions, using lambda functions
                if attr in type_convert:
                    attr_val = type_convert[attr](attr_val)

                # Convert boolean variables to str
                if type(attr_val) == np.bool_:
                    attr_val = str(attr_val)

                # Convert to lowercase
                augur_data["nodes"][c.name][attr.lower()] = attr_val

    return augur_data


def auspice_export(
    tree=None,
    augur_json_paths=None,
    auspice_config_path=None,
    auspice_colors_path=None,
    auspice_latlons_path=None,
    auspice_geo_res=None,
):
    """
    Export the full Auspice JSON v2 for visualization
    """

    if not auspice_geo_res:
        auspice_geo_res = ["country", "province"]
    export_v2.configure_warnings()

    # Initialize the auspice json
    data_json = {"version": "v2", "meta": {"updated": time.strftime("%Y-%m-%d")}}

    # parse input files
    (
        node_data,
        node_attrs,
        node_data_names,
        metadata_names,
    ) = export_v2.parse_node_data_and_metadata(tree, augur_json_paths, None)

    # data_json["tree"] = convert_tree_to_json_structure(T.root, node_attrs)

    # Validate and load config file (could put this in try except)
    export_v2.validate_auspice_config_v2(auspice_config_path)
    print("Validation success.")
    config = export_v2.read_config(auspice_config_path)

    # set metadata structures
    export_v2.set_title(data_json, config, config["title"])
    export_v2.set_display_defaults(data_json, config)
    export_v2.set_maintainers(data_json, config, config["maintainers"])
    export_v2.set_build_url(data_json, config, config["build_url"])
    export_v2.set_annotations(data_json, node_data)

    # Set Colors
    export_v2.set_colorings(
        data_json=data_json,
        config=export_v2.get_config_colorings_as_dict(config),
        node_data_colorings=node_data_names,
        provided_colors=utils.read_colors(auspice_colors_path),
        node_attrs=node_attrs,
        command_line_colorings=None,
        metadata_names=metadata_names,
    )

    # Set Filters
    export_v2.set_filters(data_json, config)

    # Set Tree
    data_json["tree"] = export_v2.convert_tree_to_json_structure(tree.root, node_attrs)
    export_v2.set_node_attrs_on_tree(data_json, node_attrs)
    export_v2.set_geo_resolutions(
        data_json,
        config,
        auspice_geo_res,
        utils.read_lat_longs(auspice_latlons_path),
        node_attrs,
    )

    # Set panels
    export_v2.set_panels(data_json, config, cmd_line_panels=None)

    return data_json


@click.command()
@click.help_option("--help", "-h")
@click.option(
    "-m",
    "--metadata",
    help="Input metadata.",
    type=click.Path(dir_okay=False, allow_dash=True),
    required=False,
)
@click.option(
    "-o",
    "--outdir",
    help="Output directory.",
    type=click.Path(dir_okay=True, allow_dash=True),
    required=False,
    default=os.getcwd(),
)
@click.option(
    "-t",
    "--tree",
    help="Input newick tree.",
    type=click.Path(exists=True, dir_okay=False, allow_dash=True),
    required=True,
)
@click.option(
    "--config",
    help="Auspice config file.",
    type=click.Path(exists=True, dir_okay=False, allow_dash=True),
    required=False,
)
@click.option(
    "--colors",
    help="Auspice colors tsv file.",
    type=click.Path(exists=True, dir_okay=False, allow_dash=True),
    required=False,
)
def main(
    tree: str, outdir: str, metadata: str, config: str, colors: str,
):
    """This script converts a newick tree to an auspice nextstrain json."""

    tree_path = tree
    metadata_path = metadata
    out_dir = outdir
    auspice_config_path = config
    auspice_colors_path = colors

    # -------------------------------------------------------------------------
    # Import Tree
    tree = Phylo.read(tree_path, "newick")
    tree.ladderize(reverse=False)

    node_i = 0
    tree_dict = {}

    for c in tree.find_clades():

        # Rename internal nodes
        if not c.name:
            c.name = "NODE{}".format(node_i)
            node_i += 1

        node_type = "internal"
        if c.is_terminal():
            node_type = "terminal"

        tree_dict[c.name] = {}
        tree_dict[c.name]["branch_length"] = c.branch_length
        tree_dict[c.name]["node_type"] = node_type

    auspice_df = pd.DataFrame.from_dict(tree_dict, orient="index")

    # -------------------------------------------------------------------------
    # Import Metadata
    if metadata_path:

        metadata_df = pd.read_csv(metadata_path, sep="\t")
        metadata_df.set_index(metadata_df.columns[0], inplace=True)

        # Add metadata to dict
        for sample in metadata_df.index:
            # Skip this sample if its not in the tree
            if sample not in tree_dict:
                continue
            for column in metadata_df.columns:
                auspice_df.at[sample, column] = str(metadata_df[column][sample])

    auspice_df.fillna(NO_DATA_CHAR, inplace=True)

    # Create a default config if one wasn't supplied
    if not auspice_config_path:
        config_dict = {
            "title": "Newick to nextstrain.",
            "build_url": "N/A",
            "maintainers": [{"name": "N/A", "url": "N/A"}],
            "panels": ["tree", "map"],
            "display_defaults": {"distance_measure": "div"},
        }
        auspice_config_path = os.path.join(out_dir, "auspice_config.json")
        with open(auspice_config_path, "w") as outfile:
            json.dump(config_dict, outfile, indent=4)

    export_v2.validate_auspice_config_v2(auspice_config_path)

    # -------------------------------------------------------------------------
    # Augur JSON
    augur_dict = augur_export(
        tree_path=None,
        aln_path=None,
        tree=tree,
        tree_df=auspice_df,
        color_keyword_exclude=["geometry"],
        type_convert={"branch_number": (lambda x: str(x))},
    )

    tree_basename = os.path.splitext(os.path.basename(tree_path))[0]
    out_path_json = os.path.join(out_dir, tree_basename + ".json")

    utils.write_json(data=augur_dict, file_name=out_path_json, indent=JSON_INDENT)

    # -------------------------------------------------------------------------
    # Auspice JSON
    auspice_dict = auspice_export(
        tree=tree,
        augur_json_paths=[out_path_json],
        auspice_config_path=auspice_config_path,
        auspice_colors_path=auspice_colors_path,
        auspice_latlons_path=None,
    )

    label_col = list(auspice_df.columns)

    # Recursively add branch attrs
    branch_attributes(
        tree_dict=auspice_dict["tree"],
        sub_dict=auspice_dict["tree"],
        df=auspice_df,
        label_col=label_col,
    )

    # Write outputs
    utils.write_json(
        data=auspice_dict,
        file_name=out_path_json,
        indent=JSON_INDENT,
        include_version=False,
    )
    export_v2.validate_data_json(out_path_json)
    print("Validation successful for json: {}\n".format(out_path_json))


if __name__ == "__main__":
    main()

import numpy as np
import copy
from augur import utils, export_v2
import time

AUSPICE_GEO_RES = ["country", "province"]


# Get X And Y Positions
def get_x_positions(tree):
    """Create a mapping of each clade to its horizontal position.
    Dict of {clade: x-coord}
    """
    depths = tree.depths()
    # If there are no branch lengths, assume unit branch lengths
    if not max(depths.values()):
        depths = tree.depths(unit_branch_lengths=True)
    return depths


# This code is from the biopython Phylo module
def get_y_positions(tree):
    """Create a mapping of each clade to its vertical position.
    Dict of {clade: y-coord}.
    Coordinates are negative, and integers for tips.
    """
    maxheight = tree.count_terminals()
    # Rows are defined by the tips
    heights = {
        tip: maxheight - i for i, tip in enumerate(reversed(tree.get_terminals()))
    }

    # Internal nodes: place at midpoint of children
    def calc_row(clade):
        for subclade in clade:
            if subclade not in heights:
                calc_row(subclade)
        # Closure over heights
        heights[clade] = (heights[clade.clades[0]] + heights[clade.clades[-1]]) / 2.0

    if tree.root.clades:
        calc_row(tree.root)
    return heights


def parse_comment(comment):
    """
    Parse nexus comment into dict
    """
    comment_dict = {}
    comment = comment.strip("[]").lstrip("&").replace('"', "")
    attr_name = ""
    attr_val = ""
    comment_sep = ","
    attr_sep = "="

    parsing_name = True
    parsing_val = False

    for i in range(0, len(comment)):

        # Parsing attribute names
        if parsing_name:
            if comment[i] != attr_sep:
                attr_name += comment[i]
            else:
                comment_dict[attr_name] = ""
                parsing_name = False
                parsing_val = True

        elif parsing_val:
            if comment[i] == comment_sep:
                # The val does not contain an interval, truly the end
                if "{" not in attr_val or ("{" in attr_val and "}" in attr_val):
                    comment_dict[attr_name] = attr_val
                    attr_val = ""
                    attr_name = ""
                    parsing_val = False
                    parsing_name = True
                else:
                    attr_val += comment[i]
            elif i == len(comment) - 1:
                attr_val += comment[i]
                comment_dict[attr_name] = attr_val
            else:
                attr_val += comment[i]

    return comment_dict


# My own conversion function
def convert_timetree_ticks(tree, step):
    """
    Return a dict of axis locations and labels for an input timetree tree.
    """
    # Step 1: Figure out offset to convert year branch length to calendar date
    if not tree.root.branch_length:
        tree.root.branch_length = 0
    min_date = tree.root.numdate - tree.root.branch_length
    max_date = np.max([n.numdate for n in tree.get_terminals()])
    offset = min_date
    date_range = max_date - min_date

    # Step 2: Relabel xticks and space them differently
    # Distance between ticks
    dtick = step
    # Minimum tick value
    min_tick = step * (offset // step)

    # Extra tick increment
    extra = dtick if dtick < date_range else dtick
    # New tick values
    tick_vals = np.arange(min_tick, min_tick + date_range + extra, dtick)
    # New tick locations
    tick_locs = list(tick_vals - offset)
    # New tick labels
    tick_labels = ["%d" % (int(x)) for x in tick_vals]
    return {"tick_locs": tick_locs, "tick_labels": tick_labels}


def metadata_to_comment(tree, tree_df):
    """
    Add metadata from a dataframe as comments in a phylogeny (modifies tree in place).
    """
    for c in tree.find_clades():
        c.comment = ""
        for col in tree_df.columns:
            col_val = tree_df[col][c.name]

            # Test for list
            if type(col_val) == list:
                col_val = "{" + "{},{}".format(col_val[0], col_val[1]) + "}"

            elif type(col_val) == str:

                # Test for not ascii
                if not col_val.isascii():
                    col_val = col_val.encode("unicode_escape")

                # Test for list as string
                elif "[" in col_val and "]" in col_val:
                    col_val = col_val.replace("[", "").replace("]", "")

                    if ":" in col_val:
                        col_val = col_val.split(":")
                    elif "," in col_val:
                        col_val = col_val.split(",")

                    # Convert to a dictionary for figree bars
                    col_val = "{" + "{},{}".format(col_val[0], col_val[1]) + "}"

            # Add comment
            if not hasattr(c, "comment") or not c.comment:
                c.comment = "&{}={}".format(col, col_val)
            else:
                c.comment += ",{}={}".format(col, col_val)


def extract_subtree(tree, tips, df, color_branches=False):
    """ Extract a subtree defined by tip names"""
    subtree_mrca = tree.common_ancestor(tips)
    subtree_clade = copy.deepcopy(subtree_mrca)

    # Search for the nodes, finding tips to prune
    for c in subtree_mrca.find_clades(order="postorder"):
        cur_term = [t.name for t in subtree_clade.get_terminals()]
        if c.name not in tips:
            if c.name in cur_term:
                subtree_clade.collapse(target=c.name)

    # Color branches
    # if color_branches:
    #    for c in subtree_clade.find_clades():
    #        state = df["mugration_branch_major"][c.name]
    #        color = colors_dict[state]
    #        c.color = color

    return subtree_clade


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

                # We need the value assigned to the trait by mugration
                if (
                    "mugration" in attr
                    and "confidence" not in attr
                    and "entropy" not in attr
                ):
                    attr = attr.replace("mugration_", "")
                    # Check again for a type conversion
                    if attr in type_convert:
                        attr_val = type_convert[attr](attr_val)
                    augur_data["nodes"][c.name][attr.lower()] = attr_val

                    # Prepare an empty dict for the confidence values
                    attr_conf = attr + "_confidence"
                    augur_data["nodes"][c.name][attr_conf.lower()] = {attr_val: "NA"}

                # Get the mugration entropy associated with the trait
                elif "mugration" in attr and "entropy" in attr:
                    attr = attr.replace("mugration_", "")
                    augur_data["nodes"][c.name][attr.lower()] = attr_val

                # Get the mugration confidence associated with the trait
                elif "mugration" in attr and "confidence" in attr:
                    attr = attr.replace("mugration_", "")
                    attr_assoc = attr.replace("_confidence", "")

                    # If original attr is not yet in dict, set to NA
                    try:
                        attr_mug_val = augur_data["nodes"][c.name][attr_assoc.lower()]

                    except KeyError:
                        attr_mug_val = "NA"
                    attr_val = {attr_mug_val: attr_val}

                    augur_data["nodes"][c.name][attr_conf.lower()] = attr_val

                # Remove the timetree prefix
                elif "timetree" in attr:
                    attr = attr.replace("timetree_", "")
                    # Catch other list types (like timetree_num_date_bar)
                    if "confidence" not in attr and type(attr_val) == list:
                        attr_val = ":".join([str(x) for x in attr_val])
                    augur_data["nodes"][c.name][attr.lower()] = attr_val

                else:
                    # Make attribute name in dict lowercase
                    # (ex. Branch_Length -> branch_length)
                    augur_data["nodes"][c.name][attr.lower()] = attr_val

    return augur_data


def auspice_export(
    tree=None,
    augur_json_paths=None,
    auspice_config_path=None,
    auspice_colors_path=None,
    auspice_latlons_path=None,
):
    """
    Export the full Auspice JSON v2 for visualization
    """
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

    # print(data_json["tree"] = convert_tree_to_json_structure(T.root, node_attrs)

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
        AUSPICE_GEO_RES,
        utils.read_lat_longs(auspice_latlons_path),
        node_attrs,
    )

    # Set panels
    export_v2.set_panels(data_json, config, cmd_line_panels=None)

    return data_json


def branch_attributes(tree_dict, sub_dict, df, label_col):
    """
    Add branch attributes to an auspice tree
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
            branch_labels = {col: df[col][node["name"]] for col in label_col}

            node["branch_attrs"]["labels"] = branch_labels

        branch_labels = {col: df[col][root["name"]] for col in label_col}
        root["branch_attrs"]["labels"] = branch_labels
        return root

import numpy as np
import copy


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

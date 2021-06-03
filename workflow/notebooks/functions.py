import numpy as np
import copy
from augur import utils, export_v2
import time
import pandas as pd

# import matplotlib.pyplot as plt
import statsmodels.api as sma
import statsmodels.formula.api as smfa

# import seaborn as sns

AUSPICE_GEO_RES = ["country", "province"]
ALPHA = 0.05


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
        AUSPICE_GEO_RES,
        utils.read_lat_longs(auspice_latlons_path),
        node_attrs,
    )

    # Set panels
    export_v2.set_panels(data_json, config, cmd_line_panels=None)

    return data_json


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
                col_pretty = (
                    col.replace("mugration_", "")
                    .replace("timetree_", "")
                    .replace("_", " ")
                    .title()
                )
                # clade needs to stay lowercase
                if col == "clade":
                    col_pretty = "clade"
                branch_labels[col_pretty] = df[col][node["name"]]

            node["branch_attrs"]["labels"] = branch_labels

        branch_labels = {}
        for col in label_col:
            col_pretty = (
                col.replace("mugration_", "")
                .replace("timetree_", "")
                .replace("_", " ")
                .title()
            )
            # clade needs to stay lowercase
            if col == "clade":
                col_pretty = "clade"
            branch_labels[col_pretty] = df[col][node["name"]]
        root["branch_attrs"]["labels"] = branch_labels
        return root


def linregress_bootstrap(
    x, y, xerr, nboots=100, s=100, plot=False, color="black", label=None, confidence=95
):
    """
    Bootstrap a linear regression.
    """

    # Add constant
    data_df = pd.DataFrame({"x": x, "y": y})

    # Construct a linear model
    ols_model = smfa.ols(formula="y ~ x", data=data_df)
    results = ols_model.fit()
    # print(results.conf_int())

    # get predicted values and residuals
    y_pred = results.predict(data_df["x"])
    resids = results.resid
    y_intercept = results.params[0]
    slope = results.params[1]
    x_intercept = (0 - y_intercept) / slope
    r_squared = results.rsquared_adj
    # print(results.summary())
    # 0 is intercept, 1 is x
    p_value = results.pvalues[1]

    bootstrap_dict = {
        "n": len(x),
        "x": x,
        "y": y,
        "slope": slope,
        "p_value": p_value,
        "r_squared": r_squared,
        "x_intercept": x_intercept,
        "y_intercept": y_intercept,
        "bootstrap_x": [],
        "bootstrap_y": [],
        "bootstrap_slopes": [],
        "bootstrap_slope_peak": [],
        "bootstrap_slope_ci": [],
        "bootstrap_slope_ci_pretty": [],
        "bootstrap_slope_kde": [],
        "bootstrap_x_intercepts": [],
        "bootstrap_x_intercept_peak": [],
        "bootstrap_x_intercept_ci": [],
        "bootstrap_x_intercept_ci_pretty": [],
        "bootstrap_x_intercept_kde": [],
        "bootstrap_y_intercepts": [],
    }

    for _ in range(nboots):
        # sample residuals with replacement
        boot_resids = np.random.choice(resids, len(x), replace=True)

        # Randomly add residuals to y values
        y_temp = [y_pred_i + resid_i for y_pred_i, resid_i in zip(y_pred, boot_resids)]

        # Store the new random coordinates
        sample_df = pd.DataFrame({"x": list(x), "y": y_temp})

        # Fit a new linear regression
        ols_model_temp = smfa.ols(formula="y ~ x", data=sample_df)
        results_temp = ols_model_temp.fit()

        bootstrap_dict["bootstrap_x"] += x
        bootstrap_dict["bootstrap_y"] += y_temp

        # get coefficients
        y_intercept = results_temp.params[0]
        slope = results_temp.params[1]
        x_intercept = (0 - y_intercept) / slope
        bootstrap_dict["bootstrap_slopes"] += [slope]
        bootstrap_dict["bootstrap_y_intercepts"] += [y_intercept]
        bootstrap_dict["bootstrap_x_intercepts"] += [x_intercept]

    # Estimate a kernel density of the slope/rate
    slope_kde = sma.nonparametric.KDEUnivariate(bootstrap_dict["bootstrap_slopes"])
    slope_kde.fit()
    peak_slope_i = np.argmax(slope_kde.density)
    peak_slope = slope_kde.support[peak_slope_i]

    slope_ci = np.array(
        np.percentile(
            np.array(bootstrap_dict["bootstrap_slopes"]),
            (100 - confidence, confidence),
            axis=0,
        )
    )
    slope_ci_pretty = (
        str(["{:.2e}".format(n) for n in slope_ci]).strip("[]").replace("'", "")
    )
    bootstrap_dict["bootstrap_slope_peak"] = peak_slope
    bootstrap_dict["bootstrap_slope_ci"] = slope_ci
    bootstrap_dict["bootstrap_slope_ci_pretty"] = slope_ci_pretty
    bootstrap_dict["bootstrap_slope_kde"] = slope_kde

    # Estimate a kernel density of the x_intercept/mrca
    x_int_kde = sma.nonparametric.KDEUnivariate(
        bootstrap_dict["bootstrap_x_intercepts"]
    )
    x_int_kde.fit()
    peak_x_int_i = np.argmax(x_int_kde.density)
    peak_x_int = x_int_kde.support[peak_x_int_i]
    x_int_ci = np.array(
        np.percentile(
            np.array(bootstrap_dict["bootstrap_x_intercepts"]),
            (100 - confidence, confidence),
            axis=0,
        )
    )
    x_int_ci_pretty = str([round(n) for n in x_int_ci]).strip("[]")
    bootstrap_dict["bootstrap_x_intercept_peak"] = peak_x_int
    bootstrap_dict["bootstrap_x_intercept_ci"] = x_int_ci
    bootstrap_dict["bootstrap_x_intercept_ci_pretty"] = x_int_ci_pretty
    bootstrap_dict["bootstrap_x_intercept_kde"] = x_int_kde

    """
    if plot:

        TARGET_RES = [1280, 240]
        DPI = 200
        FIGSIZE = [TARGET_RES[0] / DPI, TARGET_RES[1] / DPI]
        FONTSIZE = 5
        plt.rc("font", size=FONTSIZE)

        fig, axes = plt.subplots(1, 3, figsize=FIGSIZE, dpi=DPI)
        fig.subplots_adjust(wspace=0.40)
        for ax in axes:
            ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
            for spine in ax.spines:
                ax.spines[spine].set_linewidth(0.5)
        axes[0].set_title("Root-To-Tip Regression", y=1.05)
        axes[1].set_title("Substitution Rate", y=1.05)
        axes[2].set_title("MRCA Date", y=1.05)

        # -----------------------------------------------------
        # Slopes / Rates
        ax = axes[1]
        ax.plot(slope_kde.support, slope_kde.density, color="black", lw=0.5)
        ax.fill_between(slope_kde.support, slope_kde.density, color=color, alpha=0.8)
        ax.set_xlim(bootstrap_dict["bootstrap_slope_ci"])
        ax.set_ylabel("Density")
        ax.set_xlabel("Substitution Rate")

        # -----------------------------------------------------
        # X Intercept / MRCA
        ax = axes[2]
        ax.plot(x_int_kde.support, x_int_kde.density, color="black", lw=0.5)
        ax.fill_between(x_int_kde.support, x_int_kde.density, color=color, alpha=0.8)
        ax.set_xlim(bootstrap_dict["bootstrap_x_intercept_ci"])
        ax.set_ylabel("Density")
        ax.set_xlabel("Date")

        # -----------------------------------------------------
        # Linear Regression
        ax = axes[0]
        sns.regplot(
            ax=ax,
            data=data_df,
            x="x",
            y="y",
            ci=95,
            scatter_kws={"s": 0},
            line_kws={"linewidth": 0.5},
            color="grey",
        )

        ax.errorbar(
            x=x,
            y=y,
            xerr=xerr,
            yerr=None,
            ls="none",
            c=color,
            capsize=1,
            label=None,
            zorder=2,
            lw=0.5,
        )
        ax.scatter(
            x,
            y,
            s=20,
            color=color,
            ec="black",
            lw=0.50,
            label=None,
            alpha=0.8,
            zorder=3,
        )

        x_int_ci_pretty = str(
            [round(n) for n in bootstrap_dict["bootstrap_x_intercept_ci"]]
        ).strip("[]")
        slope_ci_pretty = (
            str(["{:.2e}".format(n) for n in bootstrap_dict["bootstrap_slope_ci"]])
            .strip("[]")
            .replace("'", "")
        )

        ax.annotate(
            (
                "      Clade: {}".format(label)
                + "\n"
                + "\n      R²: {}".format(round(bootstrap_dict["r_squared"], 2))
                + "\n"
                + "\n       p: {:.2e}{}".format(bootstrap_dict["p_value"], p_sig)
                + "\n"
                + "\n  Rate: {:.2e}".format(bootstrap_dict["bootstrap_slope_peak"])
                + "\n           ({})".format(slope_ci_pretty,)
                + "\n"
                + "\nMRCA: {}".format(round(bootstrap_dict["x_intercept_peak"]))
                + "\n           ({})".format(x_int_ci_pretty,)
            ),
            xy=(-1.15, 0.5),
            xycoords="axes fraction",
            size=5,
            ha="left",
            va="center",
            bbox=dict(fc="w", lw=0.5),
        )

        # Extend x-axis by 5% of the range
        xlim = ax.get_xlim()
        perc = 0.05
        x_buff = (xlim[1] - xlim[0]) * perc
        new_xlim = [xlim[0] - x_buff, xlim[1] + x_buff]
        ax.set_xlim(new_xlim)

        ax.set_xlabel("Date")
        ax.set_ylabel("Distance to Clade Root")
    """

    return bootstrap_dict

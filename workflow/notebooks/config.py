# Config

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
import geopandas
import shapely
import time
import datetime
from augur import utils, export_v2
from Bio import Phylo
import os

# ------------------------------------------------------------------------
# VARIABLES
# ------------------------------------------------------------------------

# Pandas
# pd.set_option("display.max_rows", None, "display.max_columns", None)
pd.set_option("display.max_rows", 10, "display.max_columns", None)

# Dates
CURRENT_YEAR = datetime.datetime.utcnow().year

# Branch Support Thresholds (from IQTREE docs)
ALRT_THRESH = 80
UFBOOT_THRESH = 95
SCF_THRESH = 95

# Significant digits for writing newick files
BRANCH_LEN_SIG_DIG = 12

# Data parsing
NO_DATA_CHAR = "NA"

# Treemmer
TARGET_RTL = 0.95

# Mugration Parameters
DATE_COL = "Date"
ATTRIBUTE_LIST = [
    "Branch_Number",
    "Branch_Major",
    "Branch_Minor",
    "Country",
    "Province",
]
# ATTRIBUTE_LIST = ["Branch_Number", "Branch_Major"]
MUG_CONF_THRESH = 0.95
MUG_TINY = 1e-12

# Reference Info
REF_META = {
    "Date": 1992.0,
    "DateBP": 0 - (CURRENT_YEAR - 1992.0),
    "Branch_Number": 1,
    "Branch_Major": "1.ORI",
    "Branch_Minor": "1.ORI1",
    "Country": "United States of America",
    "Province": "Colorado",
    "Biovar": "Orientalis",
    "BioSampleComment": "KEEP: Assembly Modern Reference",
    "CountryLat": 39.7837304,
    "CountryLon": -100.4458825,
    "ProvinceLat": 38.7251776,
    "ProvinceLon": -105.607716,
}
REF_DATE = 1992.0
REF_STRAIN = "CO92"


REF_LEN = 4653728

# Clock models
CONFIDENCE = 0.95
TIME_MARGINAL = True
SEQ_MARGINAL = False
MAX_ITER = 3
RELAXED_CLOCK = {"slack": 1.0, "coupling": 0}
# RELAXED_CLOCK = False
TC = "skyline"

# N_IQD Explanation
# 1_IQD is np.percentile(residuals,75) - np.percentile(residuals,25)
# 3_IQD is 3 *1_IQD
N_IQD = 3
N_STD = 2

# How to color branch supports
LOW_COL = "black"
HIGH_COL = "red"
TERM_COL = "grey"
THRESH_COL = "blue"

# Continuous data color palette
CONT_COLOR_PAL = "rainbow"

# Discrete data color palette
DISC_COLOR_PAL = "tab10"
DISC_CMAP_N = 10
DISC_CMAP = plt.get_cmap(DISC_COLOR_PAL, DISC_CMAP_N)
DISC_CMAPLIST = [DISC_CMAP(i) for i in range(DISC_CMAP.N)]
# Remove the orange from the cmap
# del DISC_CMAPLIST[1]

BLIND_CMAPLIST = ["#7b85d4", "#f37738", "#83c995"]

# Nextstrain / augur / auspice
JSON_INDENT = 2
AUSPICE_GEO_RES = ["country", "province"]
AUSPICE_PREFIX = "plague-phylogeography_"

# Plotting Graphics
figsize = (6.4, 4.8)
figsize_flip = (4.8, 6.4)
figsize_alt = (9.6, 4.8)
figsize_mini = (4.8, 2.4)
figsize_tiny = (2.4, 1.2)
dpi = 400

SM_FONT = 5
MED_FONT = 8
LG_FONT = 10

plt.rc("font", size=SM_FONT)  # controls default text sizes
plt.rc("figure", titlesize=LG_FONT)  # fontsize of the figure title
# plt.rc('axes', labelsize=MED_FONT)    # fontsize of the x and y labels
plt.rc("lines", linewidth=0.5)

plt.rc("font", size=SM_FONT)  # controls default text sizes
plt.rc("figure", titlesize=LG_FONT)  # fontsize of the figure title
plt.rc("legend", title_fontsize=MED_FONT)  # fontsize of the legend title
plt.rc("legend", frameon=False)  # legend frame
plt.rc("axes", labelsize=MED_FONT)  # fontsize of the x and y labels
plt.rc("axes", titlesize=LG_FONT)  # fontsize of axis titles
plt.rc("lines", linewidth=0.5)
plt.rc("legend", labelspacing=0.75)

FMT = "svg"

# ------------------------------------------------------------------------
# Geospatial
# ------------------------------------------------------------------------

CRS = "epsg:4326"
WEB_MERCATOR_CRS = "epsg:3857"
world_polygons = geopandas.read_file(geopandas.datasets.get_path("naturalearth_lowres"))

# Order: bottom-left, bottom-right, top-right, top-left (sw, se, ne, nw)
region_poly = {
    "Caucasus": {"wsen": [35.0000, 30.0009, 60.0000, 50.0000]},
    # Caucasus" : {"wsen": [40.058841, 40.202162, -75.042164, -74.924594]},
    "First Pandemic": {"wsen": [-10, 35, 90, 65]},
    "Second Pandemic": {"wsen": [-10, 35, 90, 65]},
    "Asia": {"wsen": [60, 0, 140, 70]},
    "1.ORI": {"wsen": [-190, -60, 190, 90]},
    "2.MED": {"wsen": [40, 30, 135, 55]},
}

df_columns = ["Region", "Lon", "Lat"]
region_df = pd.DataFrame(columns=df_columns)

for region in region_poly:
    wsen = region_poly[region]["wsen"]

    # add to dataframe
    df = pd.DataFrame(
        [
            [region, wsen[0], wsen[1]],
            [region, wsen[2], wsen[1]],
            [region, wsen[2], wsen[3]],
            [region, wsen[0], wsen[3]],
        ],
        columns=df_columns,
    )
    region_df = region_df.append(df, ignore_index=True)

    # Polygon should be: sw, se, ne, nw
    region_poly[region]["poly"] = shapely.geometry.Polygon(
        [
            (wsen[0], wsen[1]),
            (wsen[2], wsen[1]),
            (wsen[2], wsen[3]),
            (wsen[0], wsen[3]),
        ]
    )
    region_poly[region]["geoseries"] = geopandas.GeoSeries(region_poly[region]["poly"])
    region_poly[region]["geoseries"].crs = CRS
    region_poly[region]["xlim"] = (wsen[0], wsen[2])
    region_poly[region]["ylim"] = (wsen[1], wsen[3])

region_gdf = geopandas.GeoDataFrame(
    region_df, geometry=geopandas.points_from_xy(region_df.Lon, region_df.Lat)
)
region_gdf.set_crs(CRS, inplace=True)

# ------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------

# This code is from the biopython Phylo module


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


# My own conversion function
def convert_timetree_ticks(tree, step):
    """
    Return a dict of axis locations and labels for an input timetree tree.
    """
    # Step 1: Figure out offset to convert year branch length to calendar date
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
    tick_locs = tick_vals - offset
    # New tick labels
    tick_labels = ["%d" % (int(x)) for x in tick_vals]
    return {"tick_locs": tick_locs, "tick_labels": tick_labels}


def color_tree(
    tree, df, attribute, attribute_confidence, threshold_confidence, color_pal="rainbow"
):
    """
    Color branches of a tree using a data frame and addtribute.
    Returns a color dictionary and modifies the tree in place.
    """
    # Create the custom color map
    attr_states = list(dict.fromkeys(df[attribute]))

    # Create the custom color map (pyplot)
    cmap = plt.get_cmap(color_pal, len(attr_states))
    # Convert the color map to a list of RGB values
    cmaplist = [cmap(i) for i in range(cmap.N)]
    # Convert RGB values to hex colors
    attr_hex = [colors.to_hex(col) for col in cmaplist]

    hex_dict = {}

    # Assign states colors based on tip order (Low Conf first as grey)
    for state, hex_col in zip(attr_states, attr_hex):
        hex_dict[state] = hex_col

    # A flag for whether we'll add the low confidence color
    add_low_conf = False

    # Change colors on tree
    for c in tree.find_clades():
        clade_state = df[attribute][c.name]
        clade_color = hex_dict[clade_state]
        # OPTIONAL: Color grey if low confidence
        if df[attribute_confidence][c.name] < threshold_confidence:
            clade_color = "#808080"
            add_low_conf = True

        # Modify tree
        c.color = clade_color
        # Save to dictionary

    if add_low_conf:
        hex_dict["Low Confidence"] = "#808080"

    return hex_dict


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

                # We need the value assigned to the trait by mugration
                if (
                    "Mugration" in attr
                    and "Confidence" not in attr
                    and "Entropy" not in attr
                ):
                    attr = attr.replace("Mugration_", "")
                    # Check again for a type conversion
                    if attr in type_convert:
                        attr_val = type_convert[attr](attr_val)
                    augur_data["nodes"][c.name][attr.lower()] = attr_val

                    # Prepare an empty dict for the confidence values
                    attr_conf = attr + "_Confidence"
                    augur_data["nodes"][c.name][attr_conf.lower()] = {attr_val: "NA"}

                # Get the mugration entropy associated with the trait
                elif "Mugration" in attr and "Entropy" in attr:
                    attr = attr.replace("Mugration_", "")
                    augur_data["nodes"][c.name][attr.lower()] = attr_val

                # Get the mugration confidence associated with the trait
                elif "Mugration" in attr and "Confidence" in attr:
                    attr = attr.replace("Mugration_", "")
                    attr_assoc = attr.replace("_Confidence", "")

                    # If original attr is not yet in dict, set to NA
                    try:
                        attr_mug_val = augur_data["nodes"][c.name][attr_assoc.lower()]

                    except KeyError:
                        attr_mug_val = "NA"
                    attr_val = {attr_mug_val: attr_val}

                    augur_data["nodes"][c.name][attr_conf.lower()] = attr_val

                # Remove the timetree prefix
                elif "timetree" in attr and "Confidence" not in attr:
                    attr = attr.replace("timetree_", "")
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


"""
# Testing
tree_path=(
  "../../../../docs/results/latest/parse_tree/parse_tree.nwk"
)
aln_path = (
  "../../../../docs/results/latest/snippy_multi/snippy-core_chromosome.snps.filter5.aln"
)

tree_div = Phylo.read(tree_path, "newick")
tree_div.ladderize(reverse=False)
tree=tree_div

tree_df_path = "../../../../docs/results/latest/timetree/timetree.tsv"
tree_df = pd.read_csv(tree_df_path, sep='\t')
# Fix the problem with multiple forms of NA in the table
# Consolidate missing data to the NO_DATA_CHAR
tree_df.fillna(NO_DATA_CHAR, inplace=True)
tree_df.set_index("Name", inplace=True)

augur_dict = augur_export(
    tree_path=tree_path,
    aln_path=aln_path,
    tree=tree_div,
    tree_df=tree_df,
    color_keyword_exclude=["color", "coord"],
    type_convert = {
        "Branch_Number" : (lambda x : str(x))
    },
)

print(augur_dict["nodes"]["NODE0"])
"""

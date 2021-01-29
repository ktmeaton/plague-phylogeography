# Config

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
import geopandas
import shapely

# ------------------------------------------------------------------------
# VARIABLES
# ------------------------------------------------------------------------

# Pandas
# pd.set_option("display.max_rows", None, "display.max_columns", None)
pd.set_option("display.max_rows", 10, "display.max_columns", None)

# Branch Support Thresholds (from IQTREE docs)
ALRT_THRESH = 80
UFBOOT_THRESH = 95

# Significant digits for writing newick files
BRANCH_LEN_SIG_DIG = 12

# Data parsing
NO_DATA_CHAR = "NA"

# Mugration Parameters
DATE_COL = "Date"
ATTRIBUTE_LIST = ["Branch_Major", "Branch_Minor", "Country", "Province", "Biovar"]
MUG_CONF_THRESH = 0.95

# Clock models
REF_DATE = 1992.0
REF_LEN = 4653728
CONFIDENCE = 0.95
N_IQD = 3
TIME_MARGINAL = False
SEQ_MARGINAL = False
MAX_ITER = 3
RELAXED_CLOCK = {"slack": 1.0, "coupling": 0}
# RELAXED_CLOCK = False
TC = "skyline"

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

# Plotting Graphics
figsize = (6.4, 4.8)
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
    "Second Pandemic": {"wsen": [-10, 40, 90, 65]},
    "Asia": {"wsen": [60, 0, 140, 70]},
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
            clade_color = "grey"
            add_low_conf = True

        # Modify tree
        c.color = clade_color
        # Save to dictionary

    if add_low_conf:
        hex_dict["Low Confidence"] = "grey"

    return hex_dict

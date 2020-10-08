import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

# -----------------------------------------------------------------------------#
# Argument Parsing                                                             #
# -----------------------------------------------------------------------------#

parser = argparse.ArgumentParser(
    description=("Create a table image from directory."), add_help=True,
)

# Argument groups for the program

parser.add_argument(
    "--indir",
    help="Directory to search.",
    action="store",
    dest="dirPath",
    required=True,
)

parser.add_argument(
    "--outdir",
    help="Directory to store output table.",
    action="store",
    dest="outDir",
    required=True,
)

parser.add_argument(
    "--ext",
    help="Extension to search for.",
    action="store",
    dest="fileExt",
    required=True,
)

# Retrieve user parameters
args = vars(parser.parse_args())
dir_path = args["dirPath"]
out_dir = args["outDir"]
file_ext = args["fileExt"]

# -----------------------------------------------------------------------------#
# Processing                                                                   #
# -----------------------------------------------------------------------------#

sample_dict = {}
for root, _dirs, files in os.walk(dir_path, topdown=False):
    for name in files:
        if name.endswith(file_ext):
            sample_path = os.path.join(root, name)
            sample_dir = os.path.basename(os.path.dirname(sample_path))
            if sample_dir not in sample_dict:
                sample_dict[sample_dir] = []
            sample_dict[sample_dir].append(name)

df = pd.DataFrame()
df["Sample"] = [key for key in sample_dict]
df["Files"] = [", ".join(values) for keys, values in sample_dict.items()]


def render_mpl_table(
    data,
    col_width=3.0,
    row_height=0.625,
    font_size=14,
    header_color="#40466e",
    row_colors=None,
    edge_color="w",
    bbox=None,
    header_columns=0,
    ax=None,
    **kwargs
):
    """
    User: volodymyr, StackOverflow
    https://stackoverflow.com/questions/19726663/
    how-to-save-the-pandas-dataframe-series-data-as-a-figure
    """
    if bbox is None:
        bbox = [0, 0, 1, 1]
    if row_colors is None:
        row_colors = ["#f1f1f2", "w"]
    if ax is None:
        size = (np.array(data.shape[::-1]) + np.array([0, 1])) * np.array(
            [col_width, row_height]
        )
        fig, ax = plt.subplots(figsize=size)
        ax.axis("off")
    mpl_table = ax.table(
        cellText=data.values, bbox=bbox, colLabels=data.columns, **kwargs
    )
    mpl_table.auto_set_font_size(False)
    mpl_table.set_fontsize(font_size)

    for k, cell in mpl_table._cells.items():
        cell.set_edgecolor(edge_color)
        if k[0] == 0 or k[1] < header_columns:
            cell.set_text_props(weight="bold", color="w")
            cell.set_facecolor(header_color)
        else:
            cell.set_facecolor(row_colors[k[0] % len(row_colors)])
    return ax.get_figure(), ax


fig, ax = render_mpl_table(df, header_columns=0, col_width=12.0)
table_name = os.path.join(
    out_dir,
    "table_" + os.path.basename(dir_path) + "_" + file_ext.replace(".", "-") + ".pdf",
)
fig.savefig(table_name)

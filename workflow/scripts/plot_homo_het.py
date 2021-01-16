#!/usr/bin/env python3

# Import modules
import matplotlib.pyplot as plt
from matplotlib import offsetbox, ticker
import pandas as pd
import argparse
import math

# Command-line argument capture
parser = argparse.ArgumentParser(
    description="Plot histograms of homozygous and heterozygous sites.", add_help=True
)

parser.add_argument(
    "--homo",
    help="Path to the homozygosity count file.",
    action="store",
    dest="homoPath",
    required=True,
)

parser.add_argument(
    "--het",
    help="Path to the heterozygosity count file.",
    action="store",
    dest="hetPath",
    required=True,
)

# Retrieve user parameters
args = vars(parser.parse_args())
homo_path = args["homoPath"]
het_path = args["hetPath"]
out_file = "".join(het_path.split(".")[0:-2]) + ".homo_het.jpg"
sample = het_path.split(".")[0].split("/")[-1].replace("_", " ")

# Global plot settings
dpi = 400
figsize = (6.4, 4.8)

# Font
SM_FONT = 5
MED_FONT = 8
LG_FONT = 10

plt.rc("font", size=SM_FONT)  # controls default text sizes
plt.rc("figure", titlesize=LG_FONT)  # fontsize of the figure title
# plt.rc('axes', labelsize=MED_FONT)    # fontsize of the x and y labels

# Create data frames
het_df = pd.read_csv(het_path, names=["depth"])
homo_df = pd.read_csv(homo_path, names=["depth"])

# Get maximum x value (depth)
max_depth = 0
if het_df["depth"].any():
    if max(het_df["depth"]) > max_depth:
        max_depth = max(het_df["depth"])
if homo_df["depth"].any():
    if max(homo_df["depth"]) > max_depth:
        max_depth = max(homo_df["depth"])

if max_depth == 0:
    quit()
# Buffer
max_depth += 1

hist_bins = [x for x in range(0, max_depth, 1)]
xticks_major = [x for x in range(0, max_depth, 5)]
xticks_minor = [x for x in range(0, max_depth, 1)]

# Get maximum y value (count)
homo_y, homo_x, _ = plt.hist(homo_df["depth"], bins=hist_bins,)
het_y, het_x, _ = plt.hist(het_df["depth"], bins=hist_bins,)

max_count = 0
if max(het_y) > max_count:
    max_count = max(het_y)
max_count = max(het_y)
if max(homo_y) > max_count:
    max_count = max(homo_y)
# Buffer
max_count += math.ceil(0.01 * max_count)

fig, (ax1, ax2) = plt.subplots(
    2,
    sharex=False,
    sharey=False,
    gridspec_kw={"hspace": 0.25},
    # figsize=figsize,
    dpi=dpi,
    # constrained_layout=True,
)
# ---------------------------
# Axis 1: Homozygosity
ax1.hist(
    x=homo_df["depth"], bins=hist_bins, color="#1f77b4",
)

ax1.xaxis.set_major_locator(ticker.MultipleLocator(5))
ax1.xaxis.set_minor_locator(ticker.MultipleLocator(1))
ax1.set_xlabel("")
ax1.set_ylabel("Number of Sites")
ax1.set_xlim(0, max_depth)
# Option: Make yaxis the same
ax1.set_ylim(0, max_count)

at = offsetbox.AnchoredText(
    "Homozygous (n={})".format(len(homo_df["depth"])),
    prop=dict(size=MED_FONT),
    frameon=True,
    loc="upper right",
)
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax1.add_artist(at)

# ---------------------------
# Axis 2: Heterozygosity

ax2.hist(x=het_df["depth"], bins=hist_bins, color="#ff7f0e")

ax2.xaxis.set_major_locator(ticker.MultipleLocator(5))
ax2.xaxis.set_minor_locator(ticker.MultipleLocator(1))
ax2.set_xlabel("Depth (X)")
ax2.set_ylabel("Number of Sites")
# Option: Make yaxis the same
ax2.set_ylim(0, max_count)

at = offsetbox.AnchoredText(
    "Heterozygous (n={})".format(len(het_df["depth"])),
    prop=dict(size=MED_FONT),
    frameon=True,
    loc="upper right",
)
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax2.add_artist(at)

# -------------------------------
# Figure title
fig.suptitle(sample)


# Save plot
plt.savefig(out_file, dpi=dpi, bbox_inches="tight")

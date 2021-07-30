#!/usr/bin/env python3

from math import pi
from pathlib import Path
from typing import Union
from Bio import Phylo

import click
import pandas as pd
import numpy as np
from bokeh.io import save, output_file
from bokeh.models import BasicTicker, ColorBar, LinearColorMapper
from bokeh.palettes import brewer
from bokeh.plotting import figure

PathLike = Union[str, Path]
VALUE = "distance"
TOOLS = "hover,save,pan,box_zoom,reset,wheel_zoom"
TITLE = "Pairwise SNP Matrix"
WIDTH = 900
HEIGHT = 900


def load_matrix(fpath: PathLike, delim: str, filter_df) -> pd.DataFrame:
    matrix = []
    with open(fpath) as instream:
        header = next(instream).rstrip()
        ind_names = header.split(delim)[1:]
        if filter_df is not None:
            ind_names = [name for name in ind_names if name in filter_df.index]
            ind_names = [
                filter_df["Strain"][name]
                + "_"
                + str(filter_df["Date"][name]).lstrip("[").rstrip("]")
                + "_"
                + filter_df["Country"]
                for name in ind_names
            ]
            filter_i = [
                i for i in range(0, len(ind_names)) if ind_names[i] in filter_df.index
            ]

        col_names = ind_names

        for row in map(str.rstrip, instream):
            # Filter sample rows
            sample = row.split(delim)[0]
            if sample not in ind_names:
                continue
            # Filter distance column
            dists = row.split(delim)[1:]
            if filter_df:
                dists = [dists[i] for i in filter_i]

            matrix.append([int(d) for d in dists])
    return pd.DataFrame(matrix, index=ind_names, columns=col_names)


@click.command()
@click.help_option("--help", "-h")
@click.option(
    "-i",
    "--matrix",
    help="Distance matrix to plot.",
    type=click.Path(exists=True, dir_okay=False, allow_dash=True),
    default="-",
    show_default=True,
)
@click.option(
    "-o",
    "--output",
    help="Path to save HTML plot to.",
    type=click.Path(dir_okay=False, writable=True),
    default="heatmap.html",
    show_default=True,
)
@click.option(
    "-d",
    "--delim",
    help="Delimiter used in the matrix. [ default: '\\t']",
    default="\t",
    # show_default=True,
)
@click.option(
    "-p",
    "--palette",
    help="ColorBrewer palette to use for heatmap.",
    type=click.Choice(choices=brewer.keys(), case_sensitive=True),
    default="RdBu",
    show_default=True,
)
@click.option(
    "-t", "--title", help="Title for the heatmap.", default=TITLE, show_default=True
)
@click.option("--width", help="Plot width in pixels", default=WIDTH, show_default=True)
@click.option(
    "--height", help="Plot height in pixels", default=HEIGHT, show_default=True
)
@click.option(
    "-a", "--attribute", help="Attribute to filter on.", default=None,
)
@click.option(
    "-s", "--state", help="Attribute state to filter on.", default=None,
)
@click.option(
    "-m",
    "--metadata",
    help="Metadata to filter with.",
    type=click.Path(exists=True, dir_okay=False, allow_dash=True),
    default=None,
)
@click.option(
    "-p",
    "--phylogeny",
    help="Phylogeny to order samples.",
    type=click.Path(exists=True, dir_okay=False, allow_dash=True),
    default=None,
)
def main(
    matrix: str,
    output: str,
    delim: str,
    palette: str,
    title: str,
    height: int,
    width: int,
    attribute: str,
    state: str,
    metadata: str,
    phylogeny: str,
):
    """This script generates an interactive heatmap (HTML) for a distance matrix."""

    # Load metadata
    filter_df = None
    if metadata:
        metadata_df = pd.read_csv(metadata, sep="\t")
        metadata_df.fillna("NA", inplace=True)
        metadata_df.set_index("sample", inplace=True)
        # if attribute and state:
        #    filter_df = metadata_df[metadata_df[attribute] == state]
        # else:
        #    filter_df = metadata_df

    # Load distance matrix
    matrix_df = load_matrix(matrix, delim, filter_df)

    # Filter by phylogeny
    if phylogeny:
        tree = Phylo.read(phylogeny, format="newick")
        tree.ladderize(reverse=True)
        samples = [t.name for t in tree.get_terminals()]
        matrix_df = matrix_df.filter(items=samples, axis=0)
        matrix_df = matrix_df.filter(items=samples, axis=1)

    # Recode sample names
    if metadata:
        new_names = []
        for t in tree.get_terminals():
            name = metadata_df["strain"][t.name]
            new_names.append(name)

        matrix_df.index = new_names
        matrix_df.columns = new_names

        # Now filter on strains
        # matrix_df = matrix_df.filter(items=filter_strains, axis=0)
        # matrix_df = matrix_df.filter(items=filter_strains, axis=1)

    # Lower triangle
    df_lt = matrix_df.where(np.tril(np.ones(matrix_df.shape)).astype(np.bool))
    samples = list(df_lt.index)

    df = (
        # matrix_df.stack()
        df_lt.stack()
        .rename(VALUE)
        .reset_index()
        .rename(columns={"level_0": "sample1", "level_1": "sample2"})
    )

    pal_idx = max(brewer[palette].keys())
    colors = brewer[palette][pal_idx]
    mapper = LinearColorMapper(
        palette=colors, low=df[VALUE].min(), high=df[VALUE].max()
    )

    output_file(output, title=title, mode="inline")

    plot = figure(
        title=title,
        # x_range=samples,
        x_range=list(reversed(samples)),
        y_range=samples,
        # y_range=list(reversed(samples)),
        active_drag="box_zoom",
        active_scroll="wheel_zoom",
        # x_axis_location="above",
        x_axis_location="below",
        plot_width=width,
        plot_height=height,
        tools=TOOLS,
        toolbar_location="below",
        tooltips=[("pair", "@sample1 x @sample2"), (VALUE, f"@{VALUE}")],
    )

    plot.grid.grid_line_color = None
    plot.axis.axis_line_color = None
    plot.axis.major_tick_line_color = None
    # plot.axis.major_label_text_font_size = "7px"
    plot.axis.major_label_text_font_size = "10px"
    plot.axis.major_label_standoff = 0
    plot.xaxis.major_label_orientation = pi / 3

    plot.rect(
        x="sample1",
        y="sample2",
        width=1,
        height=1,
        source=df,
        fill_color={"field": VALUE, "transform": mapper},
        line_color=None,
    )

    color_bar = ColorBar(
        color_mapper=mapper,
        major_label_text_font_size="7px",
        ticker=BasicTicker(desired_num_ticks=len(colors)),
        label_standoff=6,
        border_line_color=None,
        location=(0, 0),
    )
    plot.add_layout(color_bar, "right")

    save(plot)


if __name__ == "__main__":
    main()

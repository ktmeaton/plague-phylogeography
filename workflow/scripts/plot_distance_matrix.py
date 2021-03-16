from math import pi
from pathlib import Path
from typing import Union

import click
import pandas as pd
from bokeh.io import save, output_file
from bokeh.models import BasicTicker, ColorBar, LinearColorMapper
from bokeh.palettes import brewer
from bokeh.plotting import figure

PathLike = Union[str, Path]
VALUE = "distance"
TOOLS = "hover,save,pan,box_zoom,reset,wheel_zoom"
TITLE = "Pairwise distance matrix"
WIDTH = 900
HEIGHT = 900


def load_matrix(fpath: PathLike, delim: str) -> pd.DataFrame:
    matrix = []
    with open(fpath) as instream:
        header = next(instream).rstrip()
        names = header.split(delim)[1:]
        for row in map(str.rstrip, instream):
            matrix.append(map(int, row.split(delim)[1:]))
    return pd.DataFrame(matrix, index=names, columns=names)


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
def main(
    matrix: str,
    output: str,
    delim: str,
    palette: str,
    title: str,
    height: int,
    width: int,
):
    """This script generates an interactive heatmap (HTML) for a distance matrix."""
    matrix_df = load_matrix(matrix, delim)
    samples = list(matrix_df.index)
    df = (
        matrix_df.stack()
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
        x_range=samples,
        y_range=list(reversed(samples)),
        active_drag="box_zoom",
        active_scroll="wheel_zoom",
        x_axis_location="above",
        plot_width=width,
        plot_height=height,
        tools=TOOLS,
        toolbar_location="below",
        tooltips=[("pair", "@sample1 x @sample2"), (VALUE, f"@{VALUE}")],
    )

    plot.grid.grid_line_color = None
    plot.axis.axis_line_color = None
    plot.axis.major_tick_line_color = None
    plot.axis.major_label_text_font_size = "7px"
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

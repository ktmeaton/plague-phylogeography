import plotly.graph_objects as go

data = {
    "missing_data": [0, 1, 2, 3, 4, 5],
    "all_variants": [0, 52, 396, 1786, 3653, 6166],
    "parsimony_variants": [0, 26, 162, 716, 1479, 2460],
}

fig = go.Figure()

fig.add_trace(
    go.Scatter(
        x=data["missing_data"],
        y=data["all_variants"],
        mode="lines+markers",
        name="All Variants",
        marker=dict(
            # color='LightSkyBlue',
            size=20,
            line=dict(color="DarkSlateGrey", width=2,),
        ),
        line=dict(width=5),
    )
)


fig.add_trace(
    go.Scatter(
        x=data["missing_data"],
        y=data["parsimony_variants"],
        mode="lines+markers",
        name="Parsimony Informative Variants",
        marker=dict(
            # color='LightSkyBlue',
            size=20,
            line=dict(color="DarkSlateGrey", width=2,),
        ),
        line=dict(width=5),
    )
)

fig.update_layout(
    template="simple_white",
    width=1080,
    height=720,
    title="<b>Variant Sites Across Missing Data Site Thresholds</b>",
    title_x=0.5,
    xaxis=dict(title="Missing Data Threshold Per Site (%)", tick0=0, dtick=1,),
    yaxis_title="Number of Variant Sites",
)

fig.show()

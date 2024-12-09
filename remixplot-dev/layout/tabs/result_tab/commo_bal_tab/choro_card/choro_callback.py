import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objs as go
from dash import Dash
from dash.dependencies import Output, Input

from scr.utils import Dataprep_utils
from scr.color_schema import Color_schema
from scr.functionalities import df_filter
from scr.map_tools import (
    Fig,
    get_map_center,
    get_contour_coordinates,
    get_center,
)


def plot_choro(
    gen_df_filt: pd.DataFrame,
    size: int,
    utils: Dataprep_utils,
    color_dict: dict,
) -> go.Figure:

    fig = Fig.fig_layout(size)

    node_dict: dict = {}
    for _, row in gen_df_filt.iterrows():

        reg_name: str = row["accNodesModel"]

        if reg_name in utils.mapping.keys():

            val: float = row["value"]
            val_norm: float = norm_value(val, gen_df_filt)

            map_id: str = utils.mapping[reg_name]

            geom_df = utils.geojson[utils.geojson["id"] == map_id]

            lat, lon = get_contour_coordinates(geom_df)
            node_lat, node_lon = get_center(lat, lon, map_id)
            node_dict[reg_name] = [node_lon, node_lat]

            fig.add_trace(
                Fig.fig_country(
                    lat,
                    lon,
                    reg_name,
                    color_dict,
                    min(color_dict.keys(), key=lambda x: abs(x - val_norm)),
                )
            )

    fig = get_map_center(fig, node_dict)

    return fig


def norm_value(val: float, df):
    return (val - np.min(df["value"])) / (np.max(df["value"]) - np.min(df["value"]))


def choro_callback(
    app: Dash, utils: Dataprep_utils, color_schema: Color_schema
) -> go.Figure:

    if utils.commodity_balance_df.empty != True:

        @app.callback(
            Output("choro_plot", "figure"),
            [
                Input("choro_year_dd", "value"),
                Input("choro_good_dd", "value"),
                Input("choro_tech_dd", "value"),
                Input("choro_size_fac", "value"),
            ],
        )
        def update_fig(year: str, good: str, tech: str, size: int) -> go.Figure:

            df: pd.DataFrame = utils.generation_annual_df

            color_dict = {
                round(n * 0.1, 1): color
                for n, color in enumerate(px.colors.sequential.Viridis)
            }

            bal_type_list: list[str] = df.index.get_level_values("balanceType").tolist()
            if "net" in bal_type_list:
                bal_type: str = "net"
            elif "netto" in bal_type_list:
                bal_type: str = "netto"

            gen_df_filt = df_filter(df, year, good, tech, bal_type).reset_index()

            return plot_choro(gen_df_filt, size, utils, color_dict)

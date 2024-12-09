import numpy as np
import pandas as pd
import networkx as nx
import geopandas as gpd
import plotly.graph_objects as go


from remix.framework.tools.gdx import GDXEval

from scr.map_tools import (
    get_map_center,
    get_contour_coordinates,
    get_center,
    Fig,
)


def plot_trans_flow_map(
    graph: GDXEval,
    geo_df: gpd.GeoDataFrame,
    width_fac: float,
    df: pd.DataFrame,
    mapping: dict,
    contours_bool: bool,
    nodes_bool: bool,
    connection_bool: bool,
    zoom: int,
) -> go.Figure:

    fig = Fig.fig_layout(zoom=zoom).add_trace(go.Scattergeo())

    node_dict: dict = {}
    for node in graph.nodes():
        if node in mapping.keys():
            id: str = mapping[node]
            id_df: gpd.GeoDataFrame = geo_df[(geo_df["id"] == id)]

            lat, lon = get_contour_coordinates(id_df)
            node_lat, node_lon = get_center(lat, lon, id)
            node_dict[node] = [node_lon, node_lat]

            # plot contours
            if contours_bool == True:
                node_name: str = list(id_df["NAME_LATN"])[0]
                fig.add_trace(Fig.fig_country(lat, lon, node_name))

            # plot nodes
            if nodes_bool == True:
                fig.add_trace(Fig.fig_nodes(node_lat, node_lon, node_name))

    for edge in list(graph.edges()):
        node_start: str = edge[0]
        node_end: str = edge[1]
        if node_start in node_dict.keys() and node_end in node_dict.keys():

            if width_fac != 1.5:
                df_filt = df[
                    (df["nodesModel_start"] == node_start)
                    & (df["nodesModel_end"] == node_end)
                ]
                if len(list(df_filt["value"])) != 0:

                    val: float = list(df_filt["value"])[0]
                    maxVal: float = np.max(list(df["value"]))
                    ratioVal: float = abs(val / maxVal)

                    # plot connections
                    if connection_bool == True:
                        fig.add_trace(
                            Fig.fig_trans_connections(
                                node_dict, node_start, node_end, ratioVal, width_fac
                            )
                        ).add_trace(
                            Fig.fig_trans_arrows(
                                node_dict,
                                node_start,
                                node_end,
                                ratioVal,
                                val,
                                width_fac,
                            )
                        )

            else:
                # plot connection
                if connection_bool == True:
                    fig.add_trace(Fig.fig_connections(node_dict, node_start, node_end))

    if node_dict != {}:
        fig: go.Figure = get_map_center(fig, node_dict)

    return fig


def plot_network(
    graph: nx.Graph,
    geo_df: gpd.GeoDataFrame,
    mapping_dict: dict,
    contours_bool: bool,
    nodes_bool: bool,
    connection_bool: bool,
    zoom: int,
) -> go.Figure:

    df: pd.DataFrame = pd.DataFrame({})
    return plot_trans_flow_map(
        graph,
        geo_df,
        1.5,
        df,
        mapping_dict,
        contours_bool,
        nodes_bool,
        connection_bool,
        zoom,
    )

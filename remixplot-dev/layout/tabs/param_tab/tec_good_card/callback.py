import pandas as pd
import networkx as nx
import plotly.graph_objs as go
from dash import Dash
from dash.dependencies import Output, Input

from scr.utils import Dataprep_utils
from scr.functionalities import df_filter


def size_by_freq(tg_link, techs, goods, min_size, max_size):
    tg_link_list: list = []
    for tg in tg_link:
        tg_link_list.append(tg[0])
        tg_link_list.append(tg[1])

    size_dict: dict = {}
    for tg in techs + goods:
        size_dict[tg] = tg_link_list.count(tg)

    sized_dict: dict = {}
    for element, count in size_dict.items():
        size = min_size + (count / max(size_dict.values())) * (max_size - min_size)
        sized_dict[element] = size

    sized_dict = {key: int(round(val)) for key, val in sized_dict.items()}

    return sized_dict


def plot_tg(
    graph: nx.Graph, pos: nx.spring_layout, sized_dict: dict, goods: list[str]
) -> go.Figure:

    fig = go.Figure()
    for edge in graph.edges():
        fig.add_trace(
            go.Scatter3d(
                x=[pos[edge[0]][0], pos[edge[1]][0]],
                y=[pos[edge[0]][1], pos[edge[1]][1]],
                z=[pos[edge[0]][2], pos[edge[1]][2]],
                mode="lines",
                line=dict(color="black", width=2),
                opacity=0.8,
                hoverinfo="skip",
            )
        )

    for node in graph.nodes():
        color = "blue" if node in goods else "red"
        fig.add_trace(
            go.Scatter3d(
                x=[pos[node][0]],
                y=[pos[node][1]],
                z=[pos[node][2]],
                mode="markers",
                marker=dict(size=sized_dict[node], color=color),
                hovertemplate=node,
                opacity=0.8,
            )
        )

    fig.update_layout(
        showlegend=False,
        margin=dict(t=0, b=0, r=0, l=0),
        height=600,
        scene=dict(
            xaxis=dict(showgrid=False, showticklabels=False, visible=False),
            yaxis=dict(showgrid=False, showticklabels=False, visible=False),
            zaxis=dict(showgrid=False, showticklabels=False, visible=False),
        ),
    )

    return fig


def tg_callback(app: Dash, utils: Dataprep_utils) -> go.Figure:

    @app.callback(
        Output("tg_plot", "figure"),
        [
            Input("tg_node_dd", "value"),
            Input("tg_year_dd", "value"),
            Input("tg_min_input", "value"),
            Input("tg_max_input", "value"),
        ],
    )
    def update_fig(node: str, year: str, min_size: int, max_size: int) -> go.Figure:

        commo_bal_df: pd.DataFrame = df_filter(
            utils.commodity_balance_annual_df, "gross", year, node
        ).reset_index()

        techs: list[str] = commo_bal_df["techs"].unique().tolist()
        goods: list[str] = commo_bal_df["commodity"].unique().tolist()

        tg_link: list[str] = [
            (row["techs"], row["commodity"]) for _, row in commo_bal_df.iterrows()
        ]

        graph = nx.Graph(tg_link)
        pos = nx.spring_layout(graph, dim=3)

        sized_dict = size_by_freq(
            tg_link, techs, goods, min_size=min_size, max_size=max_size
        )

        return plot_tg(graph, pos, sized_dict, goods)

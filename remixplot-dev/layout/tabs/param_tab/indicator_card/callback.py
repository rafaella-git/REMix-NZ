import pandas as pd
import plotly.graph_objs as go
from dash import Dash
from dash.dependencies import Output, Input

from scr.utils import Dataprep_utils
from scr.functionalities import df_filter, get_links_dict, restructure_df


def get_data(
    utils: Dataprep_utils, node: str, year: str, obj: str, year_a: str
) -> tuple[list[str], dict]:

    del_col: list[str] = [
        "accNodesModel_a",
        "accNodesModel",
        "accYears",
        "indicator",
        "accYears_a",
        "indicator_a",
    ]

    df: pd.DataFrame = restructure_df(
        df_filter(utils.indicator_accounting_comp, node, year)
    )
    filt_df: pd.DataFrame = df_filter(df, obj, year_a)

    nodes: list[str] = (
        filt_df.index.get_level_values("accNodesModel_a").unique().to_list()
    )

    if (
        filt_df.index.get_level_values("indicator").unique().to_list()
        != filt_df.index.get_level_values("indicator_a").unique().to_list()
    ):
        # -- level 1 --
        san1_df: pd.DataFrame = pd.DataFrame(
            {
                "source": [
                    filt_df.index.get_level_values("indicator").unique().to_list()[0]
                    for _ in nodes
                ],
                "target": nodes,
                "value": [df_filter(filt_df, node)["value"].sum() for node in nodes],
            }
        )
        # -- level 2 --
        san2_df: pd.DataFrame = filt_df.reset_index()
        san2_df["source"] = san2_df["accNodesModel_a"]
        san2_df["target"] = san2_df["indicator_a"]
        san2_df = san2_df.drop(columns=del_col, axis=0)
        # -- total --
        links: pd.DataFrame = pd.concat([san1_df, san2_df], axis=0)

    else:
        links: pd.DataFrame = filt_df.reset_index()
        links["source"] = links["accNodesModel_a"]
        links["target"] = links["indicator_a"]
        san2_df = links.drop(columns=del_col, axis=0)

    unique_source_target, links_dict = get_links_dict(links)

    return unique_source_target, links_dict


def plot_indi(unique_source_target: list[str], links_dict: dict) -> go.Figure:

    fig = go.Figure()
    fig.add_trace(
        go.Sankey(
            valueformat=".0f",
            valuesuffix="mio.€",
            # -- nodes --
            node=dict(
                pad=15,
                thickness=20,
                label=unique_source_target,
            ),
            # -- links --
            link=dict(
                source=links_dict["source"],
                target=links_dict["target"],
                value=links_dict["value"],
            ),
        )
    )
    fig.update_layout(height=600)

    return fig


# -- callback -----------------------------------------------------------------
def indicator_callback(app: Dash, utils: Dataprep_utils) -> go.Figure:

    @app.callback(
        Output("indi_plot", "figure"),
        [
            Input("indi_node_dd", "value"),
            Input("indi_year_dd", "value"),
            Input("indi_indi_dd", "value"),
            Input("indi_yeara_dd", "value"),
        ],
    )
    def update_fig(node: str, year: str, obj: str, year_a: str) -> go.Figure:

        return plot_indi(*get_data(utils, node, year, obj, year_a))

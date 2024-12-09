import pandas as pd
import plotly.graph_objs as go
from dash import Dash
from dash.dependencies import Output, Input

from scr.utils import Dataprep_utils
from scr.color_schema import Color_schema
from scr.functionalities import df_filter, get_links_dict


def thresholder(df: pd.DataFrame, th: float) -> pd.DataFrame:

    if df["value"][0] > 0:
        df_filt: pd.DataFrame = df[df["value"] >= th]
    else:
        df_filt: pd.DataFrame = df[df["value"] <= (-1) * th]

    return df_filt


def get_data(
    utils: Dataprep_utils, year: str, node: list[str] | str, th: float
) -> tuple[list[str], dict]:

    cb_an_pos_df: pd.DataFrame = df_filter(
        utils.generation_annual_df, year, node, "positive"
    ).reset_index()
    cb_an_neg_df: pd.DataFrame = df_filter(
        utils.demand_annual_df, year, node, "negative"
    ).reset_index()

    cb_an_pos_df = thresholder(cb_an_pos_df, th)
    cb_an_neg_df = thresholder(cb_an_neg_df, th)

    pos_df = pd.DataFrame(
        {
            "source": cb_an_pos_df["techs"].tolist(),
            "target": cb_an_pos_df["commodity"].tolist(),
            "value": cb_an_pos_df["value"].tolist(),
        }
    )
    neg_df = pd.DataFrame(
        {
            "source": cb_an_neg_df["commodity"].tolist(),
            "target": cb_an_neg_df["techs"].tolist(),
            "value": cb_an_neg_df["value"].abs().tolist(),
        }
    )
    links: pd.DataFrame = pd.concat([pos_df, neg_df], axis=0)

    unique_source_target, links_dict = get_links_dict(links)

    return unique_source_target, links_dict, pos_df, neg_df


def sankey_callback(
    app: Dash, utils: Dataprep_utils, color_schema: Color_schema
) -> go.Figure:

    if (
        utils.generation_annual_df.empty != True
        and utils.demand_annual_df.empty != True
    ):

        @app.callback(
            Output("sank_plot", "figure"),
            [
                Input("sank_year_dd", "value"),
                Input("sank_node_dd", "value"),
                Input("sank_thresh_input", "value"),
            ],
        )
        def update_fig(year: str, node: str | list[str], threshold: int) -> go.Figure:

            unique_source_target, links_dict, pos_df, neg_df = get_data(
                utils, year, node, threshold
            )

            color_list: list[str] = [
                (
                    color_schema.tech_color_dict[st]
                    if st not in set(list(pos_df["target"]) + list(neg_df["source"]))
                    else "blue"
                )
                for st in unique_source_target
            ]

            fig = go.Figure().update_layout(margin=dict(t=35))
            fig.add_trace(
                go.Sankey(
                    valueformat=".0f",
                    valuesuffix="GW",
                    # -- nodes --
                    node=dict(
                        pad=15,
                        thickness=20,
                        label=unique_source_target,
                        color=color_list,
                    ),
                    # -- links --
                    link=dict(
                        source=links_dict["source"],
                        target=links_dict["target"],
                        value=links_dict["value"],
                    ),
                )
            )
            fig.update_layout(height=800)

            return fig

    else:
        pass

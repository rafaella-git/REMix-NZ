import pandas as pd
import plotly.graph_objs as go
from dash import Dash
from dash.dependencies import Output, Input
from plotly.subplots import make_subplots

from scr.utils import Dataprep_utils
from scr.color_schema import Color_schema
from scr.functionalities import df_filter
from layout.tabs.result_tab.functs import transform_to_heatmap
from remix.framework.tools.gdx import fill_missing_values_in_timeseries


def timeseries_callback(
    app: Dash, utils: Dataprep_utils, color_schema: Color_schema
) -> go.Figure:

    if utils.commodity_balance_df.empty != True:

        @app.callback(
            Output("commo_bal_ts", "figure"),
            [
                Input("commo_bal_ts_tech_dd", "value"),
                Input("commo_bal_ts_good_dd", "value"),
                Input("commo_bal_ts_year_dd", "value"),
                Input("commo_bal_ts_node_dd", "value"),
            ],
        )
        def update_fig(tech: str, good: str, year: str, node: str) -> go.Figure:

            tech_list: list[str] = [tech] if isinstance(tech, str) else tech

            fig: go.Figure = make_subplots(rows=2, cols=1, vertical_spacing=0.1)

            for i, tec in enumerate(tech_list):

                df_fill: pd.DataFrame = fill_missing_values_in_timeseries(
                    df_filter(utils.commodity_balance_df, year, node, good, tec),
                    utils.timemodel_df,
                ).reset_index()

                # -- line plot
                fig.add_trace(
                    go.Scatter(
                        x=df_fill["timeModel"],
                        y=df_fill["value"],
                        mode="lines",
                        name=tec,
                        opacity=0.5,
                        line_color=color_schema.tech_color_dict[tec],
                    ),
                    row=1,
                    col=1,
                )
                # -- heatmap
                if i == 0:
                    df_hm: pd.DataFrame = transform_to_heatmap(df_fill.copy())
                    fig.add_trace(
                        go.Heatmap(
                            z=df_hm.values,
                            x=df_hm.columns,
                            y=df_hm.index,
                            colorscale="Viridis",
                            colorbar=dict(
                                len=0.5,
                                y=0.25,
                                yanchor="middle",
                                thickness=20,
                                title=f"{tec} / GW",
                            ),
                        ),
                        row=2,
                        col=1,
                    )

            fig.update_xaxes(title_text="Time / h", row=1, col=1)
            fig.update_yaxes(title_text="Value / GW", row=1, col=1)
            fig.update_xaxes(
                title_text="Days", row=2, col=1, tickmode="linear", dtick=7
            )
            fig.update_yaxes(
                title_text="Hours", row=2, col=1, tickmode="linear", dtick=3
            )
            fig.update_layout(height=800)

            return fig

    else:
        pass

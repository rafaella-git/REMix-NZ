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


def storage_level_callback(
    app: Dash, utils: Dataprep_utils, color_schema: Color_schema
) -> go.Figure:

    if utils.storage_level_df.empty != True:

        @app.callback(
            Output("stor_lev_ts", "figure"),
            [
                Input("stor_lev_ts_good_dd", "value"),
                Input("stor_lev_ts_year_dd", "value"),
                Input("stor_lev_ts_node_dd", "value"),
            ],
        )
        def update_fig(good: str, year: str, node: str) -> go.Figure:

            df = fill_missing_values_in_timeseries(
                df_filter(utils.storage_level_df, good, year, node),
                utils.timemodel_df,
            ).reset_index()
            df_hm: pd.DataFrame = transform_to_heatmap(df.copy())

            tech = df["techs"].to_list()[0]

            # figure
            fig = make_subplots(rows=2, cols=1, vertical_spacing=0.1)
            fig.update_layout(margin=dict(t=35))

            # plot 1
            fig.add_trace(
                go.Scatter(
                    x=df["timeModel"],
                    y=df["value"],
                    mode="lines",
                    line_color=color_schema.tech_color_dict[tech],
                ),
                row=1,
                col=1,
            )
            fig.update_xaxes(title_text="Time / h", row=1, col=1)
            fig.update_yaxes(title_text="Value / GW", row=1, col=1)

            # plot 2
            fig.add_trace(
                go.Heatmap(
                    z=df_hm.values,
                    x=df_hm.columns,
                    y=df_hm.index,
                    colorscale="Viridis",
                    colorbar=dict(
                        len=0.5,
                        y=0.25,
                        title="Value / GW",
                    ),
                ),
                row=2,
                col=1,
            )
            fig.update_xaxes(
                title_text="Days", tickmode="linear", dtick=7, row=2, col=1
            )
            fig.update_yaxes(
                title_text="Hours", tickmode="linear", dtick=3, row=2, col=1
            )

            fig.update_layout(height=800)

            return fig

    else:
        pass

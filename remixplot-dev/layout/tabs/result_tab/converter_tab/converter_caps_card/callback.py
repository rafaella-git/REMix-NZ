import pandas as pd
import plotly.graph_objs as go
from dash import Dash
from dash.dependencies import Output, Input

from scr.utils import Dataprep_utils
from scr.color_schema import Color_schema
from scr.functionalities import df_filter, get_tec_col
from layout.tabs.result_tab.converter_tab.converter_functs import group_techs, plot_tech


def converter_caps_callback(
    app: Dash, utils: Dataprep_utils, color_schema: Color_schema
) -> go.Figure:

    @app.callback(
        Output("conv_caps_plot", "figure"),
        [
            Input("conv_caps_node_dd", "value"),
            Input("conv_caps_good_dd", "value"),
            Input("conv_caps_captype_dd", "value"),
            Input("conv_caps_group_switch", "on"),
            Input("conv_caps_stack_switch", "on"),
        ],
    )
    def update_fig(
        node: str, good: str, cap_type: str, group: bool, stack: bool
    ) -> go.Figure:

        barmode: str = "stack" if stack == True else "group"

        converter_caps_df: pd.DataFrame = df_filter(
            utils.converter_caps_df, good, node, cap_type
        ).reset_index()

        if group == True:
            converter_caps_df = group_techs(converter_caps_df)

        tec_col: str = get_tec_col(converter_caps_df)

        fig = go.Figure().update_layout(margin=dict(t=35))

        for tec in converter_caps_df[tec_col].unique():
            df_tech = converter_caps_df[converter_caps_df[tec_col] == tec]
            fig.add_trace(plot_tech(df_tech, tec, color_schema))

        fig.update_xaxes(title_text="years")
        fig.update_yaxes(title_text="value / GW")
        fig.update_layout(barmode=barmode, height=600)
        return fig

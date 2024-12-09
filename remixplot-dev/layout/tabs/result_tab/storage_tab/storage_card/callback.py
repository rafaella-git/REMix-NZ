import pandas as pd
import plotly.graph_objs as go
from dash import Dash
from dash.dependencies import Output, Input

from scr.utils import Dataprep_utils
from scr.functionalities import df_filter, get_tec_col
from layout.tabs.result_tab.converter_tab.converter_functs import plot_tech


def storage_unit_callback(app: Dash, utils: Dataprep_utils, color_schema) -> go.Figure:

    if utils.storage_units_df.empty != True:

        @app.callback(
            Output("stor_unit_plot", "figure"),
            [
                Input("stor_unit_node_dd", "value"),
                Input("stor_unit_captype_dd", "value"),
                Input("conv_unit_stack_switch", "on"),
            ],
        )
        def update_dispatch_fig(node: str, cap_type: str, stack: bool) -> go.Figure:

            barmode: str = "stack" if stack == True else "group"

            storage_units_df_filt: pd.DataFrame = df_filter(
                utils.storage_units_df, node, cap_type
            ).reset_index()
            tec_col: str = get_tec_col(storage_units_df_filt)

            fig = go.Figure().update_layout(margin=dict(t=35))
            for tec in storage_units_df_filt[tec_col].unique():
                df_tech = storage_units_df_filt[storage_units_df_filt[tec_col] == tec]
                fig.add_trace(plot_tech(df_tech, tec, color_schema))

            fig.update_xaxes(title_text="years")
            fig.update_yaxes(title_text="value / GW")
            fig.update_layout(barmode=barmode, legend_title="Technologies", height=800)
            return fig

    else:
        pass

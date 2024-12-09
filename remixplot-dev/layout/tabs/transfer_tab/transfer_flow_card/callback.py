import pandas as pd
import plotly.graph_objs as go
from dash import Dash
from dash.dependencies import Output, Input


from scr.utils import Dataprep_utils
from layout.tabs.result_tab.functs import df_filter
from scr.transfer_functs import plot_trans_flow_map


def transfer_flow_callback(app: Dash, utils: Dataprep_utils) -> go.Figure:

    if utils.trans_flow_df.empty != True:

        @app.callback(
            Output("trans_flow_fig", "figure"),
            [
                Input("trans_flow_year_dd", "value"),
                Input("trans_flow_good_dd", "value"),
                Input("trans_flow_bal_type_dd", "value"),
                Input("trans_flow_cont_switch", "on"),
                Input("trans_flow_links_switch", "on"),
                Input("trans_flow_nodes_switch", "on"),
                Input("trans_flow_size_fac", "value"),
                Input("trans_flow_zoom_fac", "value"),
            ],
        )
        def update_network_fig(
            year: str,
            good: str,
            bal_type: str,
            cont_switch: bool,
            links_switch: bool,
            nodes_switch: bool,
            width_fac: int,
            zoom: int,
        ) -> go.Figure:

            df_filt = df_filter(utils.trans_flow_df, year, good, bal_type)

            return plot_trans_flow_map(
                utils.graph,
                utils.geojson,
                width_fac,
                df_filt,
                utils.mapping,
                contours_bool=cont_switch,
                nodes_bool=nodes_switch,
                connection_bool=links_switch,
                zoom=zoom,
            )

    else:
        pass

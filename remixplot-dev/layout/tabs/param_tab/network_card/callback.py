import plotly.graph_objs as go
from dash import Dash
from dash.dependencies import Output, Input


from scr.utils import Dataprep_utils
from scr.transfer_functs import plot_network


def network_callback(app: Dash, utils: Dataprep_utils) -> go.Figure:

    @app.callback(
        Output("netw_fig", "figure"),
        [
            Input("nw_cont_switch", "on"),
            Input("nw_links_switch", "on"),
            Input("nw_nodes_switch", "on"),
            Input("nw_zoom_fac", "value"),
        ],
    )
    def update_network_fig(
        cont_switch: bool, links_switch: bool, nodes_switch: bool, zoom: int
    ) -> go.Figure:

        return plot_network(
            utils.graph,
            utils.geojson,
            utils.mapping,
            contours_bool=cont_switch,
            nodes_bool=nodes_switch,
            connection_bool=links_switch,
            zoom=zoom,
        )

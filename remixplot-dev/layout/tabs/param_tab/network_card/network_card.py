import dash_daq as daq
import dash_bootstrap_components as dbc

from dash import dcc

from scr.color_schema import Color_schema
from layout.tabs.util_headlines import card_headline, dropdown_headline


def network_card(color_schema: Color_schema) -> dbc.Card:

    card = dbc.Card(
        dbc.CardBody(
            [
                card_headline("Network", color_schema),
                dbc.Row(
                    [
                        # Toggles
                        dbc.Col(
                            [
                                dropdown_headline("Contours", color_schema),
                                daq.BooleanSwitch(
                                    id="nw_cont_switch", on=True, className="mb-4"
                                ),
                                #
                                dropdown_headline("Links", color_schema),
                                daq.BooleanSwitch(
                                    id="nw_links_switch", on=True, className="mb-4"
                                ),
                                #
                                dropdown_headline("Nodes", color_schema),
                                daq.BooleanSwitch(
                                    id="nw_nodes_switch", on=True, className="mb-4"
                                ),
                                dropdown_headline("Zoom", color_schema),
                                dcc.Slider(
                                    id="nw_zoom_fac",
                                    min=2,
                                    max=14,
                                    step=2,
                                    value=6,
                                ),
                            ],
                            width=2,
                        ),
                        # Figure
                        dbc.Col(
                            dcc.Graph(id="netw_fig", figure={}),
                            width=10,
                            className="g-0",
                        ),
                    ]
                ),
            ]
        ),
        className="mt-3",
    )
    return card

import pandas as pd
import dash_daq as daq
import dash_bootstrap_components as dbc
from dash import dcc, html

from scr.utils import Dataprep_utils
from scr.color_schema import Color_schema
from scr.functionalities import get_idx_names
from layout.tabs.util_headlines import card_headline, dropdown_headline


def get_lists(utils: Dataprep_utils) -> None:

    df: pd.DataFrame = utils.generation_df

    nodes: list[str] = get_idx_names(df, df.index.names[1])
    years: list[str] = get_idx_names(df, df.index.names[2])
    goods: list[str] = get_idx_names(df, df.index.names[4])

    return years, goods, nodes


def dispatch_card(utils: Dataprep_utils, color_schema: Color_schema) -> dbc.Card:

    years, goods, nodes = get_lists(utils)

    card = dbc.Card(
        dbc.CardBody(
            [
                card_headline("Dispatch", color_schema),
                dbc.Row(
                    [
                        dbc.Col(
                            [
                                dropdown_headline("Commodity", color_schema),
                                dcc.Dropdown(
                                    id="dis_good_dd",
                                    multi=False,
                                    value=goods[0],
                                    options=[{"label": x, "value": x} for x in goods],
                                    className="mb-4",
                                ),
                                #
                                dropdown_headline("Node", color_schema),
                                dcc.Dropdown(
                                    id="dis_node_dd",
                                    multi=False,
                                    value="global" if "global" in nodes else nodes[0],
                                    options=[{"label": x, "value": x} for x in nodes],
                                    className="mb-4",
                                ),
                                #
                                dropdown_headline("Year", color_schema),
                                dcc.Dropdown(
                                    id="dis_year_dd",
                                    multi=False,
                                    value=years[0],
                                    options=[{"label": x, "value": x} for x in years],
                                    className="mb-4",
                                ),
                                #
                                dropdown_headline("Group technologies", color_schema),
                                daq.BooleanSwitch(
                                    id="dis_group_switch",
                                    on=False,
                                    className="mb-4",
                                ),
                            ],
                            width=2,
                        ),
                        dbc.Col(
                            [
                                dcc.Graph(id="dis_plot", figure={}),
                            ],
                            width=10,
                            className="g-0",
                        ),
                    ]
                ),
                # Heatmap
                dbc.Row(
                    [
                        dbc.Col(
                            [
                                dropdown_headline("Technology", color_schema),
                                dcc.Dropdown(
                                    id="dis_tech_dd",
                                    multi=False,
                                    value=None,
                                    options=[],
                                ),
                            ],
                            width=2,
                        ),
                        dbc.Col(
                            [
                                dcc.Store(id="dis_data"),
                                dcc.Graph(id="dis_hm_plot", figure={}),
                            ],
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

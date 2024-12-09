import pandas as pd
import dash_bootstrap_components as dbc

from dash import dcc

from scr.utils import Dataprep_utils
from scr.color_schema import Color_schema
from scr.functionalities import get_idx_names
from layout.tabs.util_headlines import card_headline, dropdown_headline


def get_lists(df: pd.DataFrame) -> None:

    nodes: list[str] = get_idx_names(df, df.index.names[0])
    years: list[str] = get_idx_names(df, df.index.names[1])
    indis: list[str] = get_idx_names(df, df.index.names[2])
    years_a: list[str] = get_idx_names(df, df.index.names[4])

    return nodes, years, indis, years_a


def indicator_card(utils: Dataprep_utils, color_schema: Color_schema) -> dbc.Card:

    nodes, years, indis, years_a = get_lists(utils.indicator_accounting_comp)

    card = dbc.Card(
        dbc.CardBody(
            [
                card_headline("Indicator", color_schema),
                dbc.Row(
                    [
                        # Dropdown
                        dbc.Col(
                            [
                                dropdown_headline("Node", color_schema),
                                dcc.Dropdown(
                                    id="indi_node_dd",
                                    multi=False,
                                    value="global",
                                    options=[{"label": x, "value": x} for x in nodes],
                                    className="mb-4",
                                ),
                                #
                                dropdown_headline("Year", color_schema),
                                dcc.Dropdown(
                                    id="indi_year_dd",
                                    multi=False,
                                    value="horizon",
                                    options=[{"label": x, "value": x} for x in years],
                                    className="mb-4",
                                ),
                                #
                                dropdown_headline("Indicator", color_schema),
                                dcc.Dropdown(
                                    id="indi_indi_dd",
                                    multi=False,
                                    value="SystemCost",
                                    options=[{"label": x, "value": x} for x in indis],
                                    className="mb-4",
                                ),
                                #
                                dropdown_headline("Year_a", color_schema),
                                dcc.Dropdown(
                                    id="indi_yeara_dd",
                                    multi=False,
                                    value=years_a[0],
                                    options=[{"label": x, "value": x} for x in years_a],
                                    className="mb-4",
                                ),
                            ],
                            width=2,
                        ),
                        # Figure
                        dbc.Col(
                            dcc.Graph(id="indi_plot", figure={}),
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

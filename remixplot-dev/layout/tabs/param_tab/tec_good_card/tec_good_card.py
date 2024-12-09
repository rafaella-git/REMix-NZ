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

    return nodes, years


def tec_good_card(utils: Dataprep_utils, color_schema: Color_schema) -> dbc.Card:

    nodes, years = get_lists(utils.commodity_balance_annual_df)

    card = dbc.Card(
        dbc.CardBody(
            [
                card_headline("Technologies + Commodities", color_schema),
                dbc.Row(
                    [
                        # Dropdown
                        dbc.Col(
                            [
                                dropdown_headline("Node", color_schema),
                                dcc.Dropdown(
                                    id="tg_node_dd",
                                    multi=False,
                                    value="global",
                                    options=[{"label": x, "value": x} for x in nodes],
                                    className="mb-4",
                                ),
                                #
                                dropdown_headline("Year", color_schema),
                                dcc.Dropdown(
                                    id="tg_year_dd",
                                    multi=False,
                                    value=years[-1],
                                    options=[{"label": x, "value": x} for x in years],
                                    className="mb-4",
                                ),
                                #
                                dropdown_headline("Minium size", color_schema),
                                dcc.Input(
                                    id="tg_min_input",
                                    type="number",
                                    value=2,
                                    min=2,
                                    max=10,
                                    step=1,
                                ),
                                #
                                dropdown_headline("Maximum size", color_schema),
                                dcc.Input(
                                    id="tg_max_input",
                                    type="number",
                                    value=20,
                                    min=11,
                                    max=30,
                                    step=1,
                                ),
                            ],
                            width=2,
                        ),
                        # Figure
                        dbc.Col(
                            dcc.Graph(id="tg_plot", figure={}),
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

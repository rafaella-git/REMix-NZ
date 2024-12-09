import pandas as pd
import dash_daq as daq
import dash_bootstrap_components as dbc
from dash import dcc

from scr.utils import Dataprep_utils
from scr.color_schema import Color_schema
from layout.tabs.util_headlines import card_headline, dropdown_headline


def get_lists(df: pd.DataFrame) -> None:

    if "accYears" in df.columns:
        years: list[str] = list(df["accYears"].unique())
    elif "years" in df.columns:
        years: list[str] = list(df["years"].unique())

    goods: list[str] = list(df["commodity"].unique())
    techs: list[str] = list(df["techs"].unique())
    bal_types: list[str] = list(df["balanceType"].unique())

    return years, goods, techs, bal_types


def transfer_flow_card(utils: Dataprep_utils, color_schema: Color_schema) -> dbc.Card:

    if utils.trans_flow_df.empty != True:

        years, goods, techs, bal_types = get_lists(utils.trans_flow_df)

        card = dbc.Card(
            dbc.CardBody(
                [
                    card_headline("Transfer flow", color_schema),
                    dbc.Row(
                        [
                            dbc.Col(
                                [
                                    #
                                    dropdown_headline("Year", color_schema),
                                    dcc.Dropdown(
                                        id="trans_flow_year_dd",
                                        multi=False,
                                        value=years[0],
                                        options=[
                                            {"label": x, "value": x} for x in years
                                        ],
                                        className="mb-4",
                                    ),
                                    #
                                    dropdown_headline("Commodity", color_schema),
                                    dcc.Dropdown(
                                        id="trans_flow_good_dd",
                                        multi=False,
                                        value=goods[0],
                                        options=[
                                            {"label": x, "value": x} for x in goods
                                        ],
                                        className="mb-4",
                                    ),
                                    #
                                    dropdown_headline("Balance Type", color_schema),
                                    dcc.Dropdown(
                                        id="trans_flow_bal_type_dd",
                                        multi=False,
                                        value=bal_types[0],
                                        options=[
                                            {"label": x, "value": x} for x in bal_types
                                        ],
                                        className="mb-4",
                                    ),
                                ],
                                width=2,
                            ),
                            # Figure
                            dbc.Col(
                                dcc.Graph(id="trans_flow_fig", figure={}),
                                width=10,
                                className="g-0",
                            ),
                        ],
                        className="mb-4",
                    ),
                    dbc.Row(
                        [
                            dbc.Col(
                                [
                                    dropdown_headline("Contours", color_schema),
                                    daq.BooleanSwitch(
                                        id="trans_flow_cont_switch",
                                        on=True,
                                        className="me-1",
                                    ),
                                ],
                                width=2,
                            ),
                            dbc.Col(
                                [
                                    dropdown_headline("Links", color_schema),
                                    daq.BooleanSwitch(
                                        id="trans_flow_links_switch",
                                        on=True,
                                        className="me-1",
                                    ),
                                ],
                                width=2,
                            ),
                            dbc.Col(
                                [
                                    dropdown_headline("Nodes", color_schema),
                                    daq.BooleanSwitch(
                                        id="trans_flow_nodes_switch",
                                        on=True,
                                        className="me-1",
                                    ),
                                ],
                                width=2,
                            ),
                            dbc.Col(
                                [
                                    dropdown_headline("Arrow size", color_schema),
                                    dcc.Slider(
                                        id="trans_flow_size_fac",
                                        min=2,
                                        max=8,
                                        step=1,
                                        value=5,
                                    ),
                                ],
                                width=2,
                            ),
                            dbc.Col(
                                [
                                    dropdown_headline("Zoom", color_schema),
                                    dcc.Slider(
                                        id="trans_flow_zoom_fac",
                                        min=2,
                                        max=14,
                                        step=2,
                                        value=6,
                                    ),
                                ],
                                width=2,
                            ),
                        ]
                    ),
                ]
            ),
            className="mt-3",
        )
        return card

    else:
        pass

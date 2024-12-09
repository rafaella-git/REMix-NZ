import pandas as pd
import dash_daq as daq
import dash_bootstrap_components as dbc
from dash import dcc

from scr.utils import Dataprep_utils
from scr.color_schema import Color_schema
from scr.functionalities import get_idx_names
from layout.tabs.util_headlines import card_headline, dropdown_headline


def get_lists(df: pd.DataFrame) -> None:

    goods: list[str] = get_idx_names(df, df.index.names[3])
    nodes: list[str] = get_idx_names(df, df.index.names[0])

    return goods, nodes


def pro_con_card(utils: Dataprep_utils, color_schema: Color_schema) -> dbc.Card:

    if (
        utils.generation_annual_df.empty != True
        and utils.demand_annual_df.empty != True
    ):

        goods, nodes = get_lists(utils.generation_annual_df)

        card = dbc.Card(
            dbc.CardBody(
                [
                    card_headline("Energy Producer/Consumer (net)", color_schema),
                    dbc.Row(
                        [
                            # Dropdown
                            dbc.Col(
                                [
                                    #
                                    dropdown_headline("Commodity", color_schema),
                                    dcc.Dropdown(
                                        id="pro_con_good_dd",
                                        multi=False,
                                        value=goods[0],
                                        options=[
                                            {"label": x, "value": x} for x in goods
                                        ],
                                        className="mb-4",
                                    ),
                                    #
                                    dropdown_headline("Node", color_schema),
                                    dcc.Dropdown(
                                        id="pro_con_node_dd",
                                        multi=False,
                                        value=(
                                            "global" if "global" in nodes else nodes[0]
                                        ),
                                        options=[
                                            {"label": x, "value": x} for x in nodes
                                        ],
                                        className="mb-4",
                                    ),
                                    #
                                    dropdown_headline(
                                        "Group technologies", color_schema
                                    ),
                                    daq.BooleanSwitch(
                                        id="pro_con_group_switch",
                                        on=False,
                                        className="mb-4",
                                    ),
                                    #
                                    dropdown_headline("Stacking", color_schema),
                                    daq.BooleanSwitch(
                                        id="pro_con_stack_switch",
                                        on=False,
                                        className="mb-4",
                                    ),
                                ],
                                width=2,
                            ),
                            # Figure
                            dbc.Col(
                                dcc.Graph(id="pro_con_plot", figure={}),
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

import pandas as pd
import dash_daq as daq
import dash_bootstrap_components as dbc
from dash import dcc

from scr.utils import Dataprep_utils
from scr.color_schema import Color_schema
from scr.functionalities import get_idx_names
from layout.tabs.util_headlines import card_headline, dropdown_headline


def get_lists(utils: Dataprep_utils) -> tuple[list[str]]:

    df: pd.DataFrame = utils.storage_units_df

    nodes: list[str] = get_idx_names(df, df.index.names[0])
    cap_types: list[str] = get_idx_names(df, df.index.names[3])

    return nodes, cap_types


def storage_unit_card(utils: Dataprep_utils, color_schema: Color_schema) -> dbc.Card:

    if utils.storage_units_df.empty != True:

        nodes, cap_types = get_lists(utils)

        card = dbc.Card(
            dbc.CardBody(
                [
                    card_headline("Storage units", color_schema),
                    dbc.Row(
                        [
                            # Dropdown
                            dbc.Col(
                                [
                                    #
                                    dropdown_headline("Node", color_schema),
                                    dcc.Dropdown(
                                        id="stor_unit_node_dd",
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
                                    dropdown_headline("Capacity type", color_schema),
                                    dcc.Dropdown(
                                        id="stor_unit_captype_dd",
                                        multi=False,
                                        value=cap_types[0],
                                        options=[
                                            {"label": x, "value": x} for x in cap_types
                                        ],
                                        className="mb-4",
                                    ),
                                    #
                                    dropdown_headline("Stacking", color_schema),
                                    daq.BooleanSwitch(
                                        id="stor_unit_stack_switch",
                                        on=False,
                                        className="mb-4",
                                    ),
                                ],
                                width=2,
                            ),
                            # Figure
                            dbc.Col(
                                dcc.Graph(id="stor_unit_plot", figure={}),
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

    else:
        pass

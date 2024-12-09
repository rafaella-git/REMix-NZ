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


def sankey_card(utils: Dataprep_utils, color_schema: Color_schema) -> dbc.Card:

    if (
        utils.generation_annual_df.empty != True
        and utils.demand_annual_df.empty != True
    ):

        nodes, years = get_lists(utils.generation_annual_df)

        card = dbc.Card(
            dbc.CardBody(
                [
                    card_headline("Commodity flow", color_schema),
                    dbc.Row(
                        [
                            # Dropdown
                            dbc.Col(
                                [
                                    #
                                    dropdown_headline("Year", color_schema),
                                    dcc.Dropdown(
                                        id="sank_year_dd",
                                        multi=False,
                                        value=years[0],
                                        options=[
                                            {"label": x, "value": x} for x in years
                                        ],
                                        className="mb-4",
                                    ),
                                    #
                                    dropdown_headline("Node", color_schema),
                                    dcc.Dropdown(
                                        id="sank_node_dd",
                                        multi=True,
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
                                        "Numerical threshold", color_schema
                                    ),
                                    dcc.Input(
                                        id="sank_thresh_input",
                                        type="number",
                                        value=0,
                                        min=0,
                                        max=1_000_000,
                                        step=1,
                                    ),
                                ],
                                width=2,
                            ),
                            # Figure
                            dbc.Col(
                                dcc.Graph(id="sank_plot", figure={}),
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

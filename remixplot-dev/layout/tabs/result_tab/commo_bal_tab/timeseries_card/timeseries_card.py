import dash_bootstrap_components as dbc
from dash import dcc

from scr.utils import Dataprep_utils
from scr.functionalities import get_idx_names
from layout.tabs.util_headlines import card_headline, dropdown_headline


def get_lists(utils: Dataprep_utils) -> tuple[list[str]]:

    nodes: list[str] = get_idx_names(
        utils.commodity_balance_df, utils.commodity_balance_df.index.names[1]
    )
    years: list[str] = get_idx_names(
        utils.commodity_balance_df, utils.commodity_balance_df.index.names[2]
    )
    techs: list[str] = get_idx_names(
        utils.commodity_balance_df, utils.commodity_balance_df.index.names[3]
    )
    goods: list[str] = get_idx_names(
        utils.commodity_balance_df, utils.commodity_balance_df.index.names[4]
    )

    return techs, goods, years, nodes


def timeseries_card(utils: Dataprep_utils, color_schema) -> dbc.Card:

    if utils.commodity_balance_df.empty != True:

        techs, goods, years, nodes = get_lists(utils)

        card = dbc.Card(
            dbc.CardBody(
                [
                    card_headline("Commodity balance time series", color_schema),
                    dbc.Row(
                        [
                            # Dropdown
                            dbc.Col(
                                [
                                    dropdown_headline("Technology", color_schema),
                                    dcc.Dropdown(
                                        id="commo_bal_ts_tech_dd",
                                        multi=True,
                                        value=techs[0],
                                        options=[
                                            {"label": x, "value": x} for x in techs
                                        ],
                                        className="mb-4",
                                    ),
                                    #
                                    dropdown_headline("Commodity", color_schema),
                                    dcc.Dropdown(
                                        id="commo_bal_ts_good_dd",
                                        multi=False,
                                        value=goods[0],
                                        options=[
                                            {"label": x, "value": x} for x in goods
                                        ],
                                        className="mb-4",
                                    ),
                                    #
                                    dropdown_headline("Year", color_schema),
                                    dcc.Dropdown(
                                        id="commo_bal_ts_year_dd",
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
                                        id="commo_bal_ts_node_dd",
                                        multi=False,
                                        value=nodes[0],
                                        options=[
                                            {"label": x, "value": x} for x in nodes
                                        ],
                                        className="mb-4",
                                    ),
                                ],
                                width=2,
                            ),
                            # Figure
                            dbc.Col(
                                dcc.Graph(id="commo_bal_ts", figure={}),
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

import pandas as pd
import dash_bootstrap_components as dbc
from dash import dcc

from scr.utils import Dataprep_utils
from scr.color_schema import Color_schema
from scr.functionalities import get_idx_names
from layout.tabs.util_headlines import card_headline, dropdown_headline


def get_lists(df: pd.DataFrame) -> None:

    goods: list[str] = get_idx_names(df, df.index.names[3])
    years: list[str] = get_idx_names(df, df.index.names[1])
    techs: list[str] = get_idx_names(df, df.index.names[2])

    return years, goods, techs


def choro_card(utils: Dataprep_utils, color_schema: Color_schema) -> dbc.Card:

    if (
        utils.generation_annual_df.empty != True
        and utils.demand_annual_df.empty != True
    ):

        years, goods, techs = get_lists(utils.generation_annual_df)

        card = dbc.Card(
            dbc.CardBody(
                [
                    card_headline("Choropleth plot", color_schema),
                    dbc.Row(
                        [
                            # Dropdown
                            dbc.Col(
                                [
                                    dropdown_headline("Year", color_schema),
                                    dcc.Dropdown(
                                        id="choro_year_dd",
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
                                        id="choro_good_dd",
                                        multi=False,
                                        value=goods[0],
                                        options=[
                                            {"label": x, "value": x} for x in goods
                                        ],
                                        className="mb-4",
                                    ),
                                    #
                                    dropdown_headline("Technology", color_schema),
                                    dcc.Dropdown(
                                        id="choro_tech_dd",
                                        multi=False,
                                        value=techs[0],
                                        options=[
                                            {"label": x, "value": x} for x in techs
                                        ],
                                        className="mb-4",
                                    ),
                                    #
                                    dropdown_headline("Size", color_schema),
                                    dcc.Slider(
                                        id="choro_size_fac",
                                        min=2,
                                        max=14,
                                        step=2,
                                        value=6,
                                    ),
                                ]
                            ),
                            # Figure
                            dbc.Col(
                                dcc.Graph(id="choro_plot", figure={}),
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

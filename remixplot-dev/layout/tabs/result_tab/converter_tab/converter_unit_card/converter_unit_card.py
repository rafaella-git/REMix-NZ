import dash_daq as daq
import dash_bootstrap_components as dbc
from dash import dcc

from scr.utils import Dataprep_utils
from scr.color_schema import Color_schema
from layout.tabs.util_headlines import card_headline, dropdown_headline

from scr.functionalities import get_idx_names


def get_lists(utils: Dataprep_utils) -> tuple[list[str]]:

    df = utils.converter_units_df

    nodes: list[str] = get_idx_names(df, df.index.names[0])
    years: list[str] = get_idx_names(df, df.index.names[1])
    cap_types: list[str] = get_idx_names(df, df.index.names[4])

    return nodes, years, cap_types


def converter_unit_card(utils: Dataprep_utils, color_schema: Color_schema) -> dbc.Card:

    nodes, years, cap_types = get_lists(utils)

    card = dbc.Card(
        dbc.CardBody(
            [
                card_headline("Converter units", color_schema),
                dbc.Row(
                    [
                        # Dropdown
                        dbc.Col(
                            [
                                #
                                dropdown_headline("Node", color_schema),
                                dcc.Dropdown(
                                    id="conv_unit_node_dd",
                                    multi=False,
                                    value="global" if "global" in nodes else nodes[0],
                                    options=[{"label": x, "value": x} for x in nodes],
                                    className="mb-4",
                                ),
                                #
                                dropdown_headline("Capacity type", color_schema),
                                dcc.Dropdown(
                                    id="conv_unit_captype_dd",
                                    multi=False,
                                    value=cap_types[0],
                                    options=[
                                        {"label": x, "value": x} for x in cap_types
                                    ],
                                    className="mb-4",
                                ),
                                #
                                dropdown_headline("Group technologies", color_schema),
                                daq.BooleanSwitch(
                                    id="conv_unit_group_switch",
                                    on=False,
                                    className="mb-4",
                                ),
                                #
                                dropdown_headline("Stacking", color_schema),
                                daq.BooleanSwitch(
                                    id="conv_unit_stack_switch",
                                    on=False,
                                    className="mb-4",
                                ),
                            ],
                            width=2,
                        ),
                        # Figure
                        dbc.Col(
                            dcc.Graph(id="conv_unit_plot", figure={}),
                            width=10,
                            className="g-0",
                        ),
                    ],
                ),
                dbc.Row(
                    [
                        # Dropdown
                        dbc.Col(
                            [
                                #
                                dropdown_headline("Year", color_schema),
                                dcc.Dropdown(
                                    id="conv_unit_year_dd",
                                    multi=True,
                                    value="horizon" if "horizon" in nodes else years[0],
                                    options=[{"label": x, "value": x} for x in years],
                                    className="mb-4",
                                ),
                                #
                                dropdown_headline("Percentage Others", color_schema),
                                dcc.Input(
                                    id="conv_unit_pie_th",
                                    type="number",
                                    value=5,
                                    min=1,
                                    max=100,
                                    step=1,
                                ),
                            ],
                            width=2,
                        ),
                        # Figure
                        dbc.Col(
                            [
                                dcc.Store(id="conv_unit_data"),
                                dcc.Graph(id="conv_unit_pie_plot", figure={}),
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

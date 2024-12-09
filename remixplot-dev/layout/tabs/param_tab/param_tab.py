import dash_bootstrap_components as dbc
from dash import Dash

from scr.utils import Dataprep_utils
from scr.color_schema import Color_schema

from layout.tabs.param_tab.network_card.network_card import network_card
from layout.tabs.param_tab.network_card.callback import network_callback
from layout.tabs.param_tab.tec_good_card.tec_good_card import tec_good_card

from layout.tabs.param_tab.indicator_card.indicator_card import indicator_card
from layout.tabs.param_tab.indicator_card.callback import indicator_callback
from layout.tabs.param_tab.tec_good_card.callback import tg_callback


def param_cards(utils: Dataprep_utils, color_schema: Color_schema) -> None:

    return dbc.Tab(
        [
            network_card(color_schema),
            indicator_card(utils, color_schema),
            tec_good_card(utils, color_schema),
        ],
        label="Parameterization",
        activeTabClassName="fw-bold",
        label_style={"color": color_schema.green},
    )


def param_callbacks(app: Dash, utils: Dataprep_utils) -> None:

    network_callback(app, utils)
    indicator_callback(app, utils)
    tg_callback(app, utils)

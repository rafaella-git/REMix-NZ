import dash_bootstrap_components as dbc
from dash import Dash

from scr.utils import Dataprep_utils
from scr.color_schema import Color_schema

from layout.tabs.transfer_tab.transfer_flow_card.transfer_flow_card import (
    transfer_flow_card,
)
from layout.tabs.transfer_tab.transfer_flow_card.callback import transfer_flow_callback


def transfer_cards(utils: Dataprep_utils, color_schema: Color_schema) -> None:

    return dbc.Tab(
        [
            # transfer cards
            transfer_flow_card(utils, color_schema)
        ],
        label="Transfer",
        activeTabClassName="fw-bold",
        label_style={"color": color_schema.green},
    )


def transfer_callbacks(app: Dash, utils: Dataprep_utils) -> None:
    transfer_flow_callback(app, utils)

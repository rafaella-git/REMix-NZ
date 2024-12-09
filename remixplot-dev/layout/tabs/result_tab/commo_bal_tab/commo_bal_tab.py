import dash_bootstrap_components as dbc
from dash import Dash

from scr.utils import Dataprep_utils
from scr.color_schema import Color_schema
from layout.tabs.result_tab.commo_bal_tab.dispatch_card.dispatch_card import (
    dispatch_card,
)
from layout.tabs.result_tab.commo_bal_tab.pro_con_card.pro_con_card import pro_con_card
from layout.tabs.result_tab.commo_bal_tab.sankey_card.sankey_card import sankey_card
from layout.tabs.result_tab.commo_bal_tab.choro_card.choro_card import choro_card
from layout.tabs.result_tab.commo_bal_tab.dispatch_card.callback import (
    dispatch_callback,
)
from layout.tabs.result_tab.commo_bal_tab.pro_con_card.callback import pro_con_callback
from layout.tabs.result_tab.commo_bal_tab.sankey_card.callback import sankey_callback
from layout.tabs.result_tab.commo_bal_tab.choro_card.choro_callback import (
    choro_callback,
)


def commo_bal_tab(utils: Dataprep_utils, color_schema: Color_schema) -> dbc.Tab:
    return dbc.Tab(
        [
            dispatch_card(utils, color_schema),
            pro_con_card(utils, color_schema),
            sankey_card(utils, color_schema),
            choro_card(utils, color_schema),
        ],
        label="Commo. Balance",
        activeTabClassName="fw-bold",
        label_style={"color": color_schema.blue},
    )


def commo_bal_callbacks(
    app: Dash, utils: Dataprep_utils, color_schema: Color_schema
) -> None:

    dispatch_callback(app, utils, color_schema)
    pro_con_callback(app, utils, color_schema)
    sankey_callback(app, utils, color_schema)
    choro_callback(app, utils, color_schema)

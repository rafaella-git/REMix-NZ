import dash_bootstrap_components as dbc
from dash import Dash

from scr.utils import Dataprep_utils
from scr.color_schema import Color_schema

from layout.tabs.result_tab.storage_tab.storage_card.storage_card import (
    storage_unit_card,
)
from layout.tabs.result_tab.storage_tab.storage_level_card.storage_level_card import (
    storage_level_card,
)


def storage_tab(utils: Dataprep_utils, color_schema: Color_schema) -> dbc.Tab:

    return dbc.Tab(
        [
            storage_unit_card(utils, color_schema),
            storage_level_card(utils, color_schema),
        ],
        label="Storage",
        activeTabClassName="fw-bold",
        label_style={"color": color_schema.blue},
    )


from layout.tabs.result_tab.storage_tab.storage_card.callback import (
    storage_unit_callback,
)
from layout.tabs.result_tab.storage_tab.storage_level_card.callback import (
    storage_level_callback,
)


def storage_callbacks(
    app: Dash, utils: Dataprep_utils, color_schema: Color_schema
) -> None:
    storage_unit_callback(app, utils, color_schema)
    storage_level_callback(app, utils, color_schema)

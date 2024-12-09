import dash_bootstrap_components as dbc
from dash import Dash

from scr.utils import Dataprep_utils
from scr.color_schema import Color_schema
from layout.tabs.result_tab.converter_tab.converter_unit_card.converter_unit_card import (
    converter_unit_card,
)
from layout.tabs.result_tab.converter_tab.converter_caps_card.converter_caps_card import (
    converter_caps_card,
)
from layout.tabs.result_tab.converter_tab.converter_unit_card.callback import (
    converter_unit_callback,
)
from layout.tabs.result_tab.converter_tab.converter_caps_card.callback import (
    converter_caps_callback,
)


def converter_tab(utils: Dataprep_utils, color_schema: Color_schema) -> dbc.Tab:

    return dbc.Tab(
        [
            converter_unit_card(utils, color_schema),
            converter_caps_card(utils, color_schema),
        ],
        label="Converter",
        activeTabClassName="fw-bold",
        label_style={"color": color_schema.blue},
    )


def converter_callbacks(
    app: Dash, utils: Dataprep_utils, color_schema: Color_schema
) -> None:
    converter_unit_callback(app, utils, color_schema)
    converter_caps_callback(app, utils, color_schema)

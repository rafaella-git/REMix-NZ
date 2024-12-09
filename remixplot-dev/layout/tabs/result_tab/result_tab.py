import dash_bootstrap_components as dbc
from dash import Dash

from scr.utils import Dataprep_utils
from scr.color_schema import Color_schema


from layout.tabs.result_tab.commo_bal_tab.commo_bal_tab import (
    commo_bal_tab,
    commo_bal_callbacks,
)
from layout.tabs.result_tab.converter_tab.converter_tab import (
    converter_tab,
    converter_callbacks,
)
from layout.tabs.result_tab.storage_tab.storage_tab import (
    storage_tab,
    storage_callbacks,
)


def result_tabs(utils: Dataprep_utils, color_schema: Color_schema) -> dbc.Tabs:

    return dbc.Tab(
        dbc.Tabs(
            [
                commo_bal_tab(utils, color_schema),
                converter_tab(utils, color_schema),
                storage_tab(utils, color_schema),
            ]
        ),
        label="Results",
        activeTabClassName="fw-bold pt-2",
        label_style={"color": color_schema.green},
    )


def result_callbacks(
    app: Dash, utils: Dataprep_utils, color_schema: Color_schema
) -> None:

    commo_bal_callbacks(app, utils, color_schema)
    converter_callbacks(app, utils, color_schema)
    storage_callbacks(app, utils, color_schema)

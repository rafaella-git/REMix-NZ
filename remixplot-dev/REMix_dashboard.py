import yaml
import dash
import dash_bootstrap_components as dbc

from dash import Dash, html, dcc
from dash.dependencies import Input, Output, State


from layout.headline.headline import headline
from layout.tabs.param_tab.param_tab import param_callbacks, param_cards
from layout.tabs.result_tab.result_tab import (
    result_callbacks,
    result_tabs,
)
from layout.tabs.transfer_tab.transfer_tab import transfer_callbacks, transfer_cards
from scr.color_schema import Color_schema
from scr.utils import Dataprep_utils


def scen_card():
    return dbc.Card(
        dbc.CardBody(),
        className="mt-3",
    )


def main():

    # -- load config ----------------------------------------------------------
    with open("c:/Local/REMix/remixplot-dev/config.yaml") as file:
        cfg = yaml.load(file, Loader=yaml.Loader)

    # -- create app -----------------------------------------------------------
    app = Dash(__name__, external_stylesheets=[dbc.themes.LUMEN])

    # -- load data in ---------------------------------------------------------
    utils: Dataprep_utils = Dataprep_utils(
        gdx_path=cfg["gdx_path"],
        geojson_path=cfg["geojson_path"],
        mapping_path=cfg["mapping_path"],
    )

    color_schema: Color_schema = Color_schema(utils)

    # -- app layout -----------------------------------------------------------
    app.layout = dbc.Container(
        [
            headline(color_schema),
            # -- Base ------------------------------
            dbc.Tabs(
                [
                    param_cards(utils, color_schema),
                    transfer_cards(utils, color_schema),
                    result_tabs(utils, color_schema),
                ]
            ),
        ]
    )

    # -- callbacks ------------------------------------------------------------
    param_callbacks(app, utils)
    transfer_callbacks(app, utils)
    result_callbacks(app, utils, color_schema)

    # -- run app --------------------------------------------------------------
    app.run_server(debug=cfg["debug"], host=cfg["host"], port=cfg["port"])


if __name__ == "__main__":
    main()

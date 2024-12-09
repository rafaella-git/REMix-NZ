import dash_bootstrap_components as dbc
from dash import html
from PIL import Image

from scr.color_schema import Color_schema

remix_logo_path: str = "./remixplot-dev/layout/headline/logos/REMix_logo_text.PNG" #
remix_img = Image.open(remix_logo_path)

dlr_logo_path: str = "./remixplot-dev/layout/headline/logos/DLR_Signet_black.png"
dlr_img = Image.open(dlr_logo_path)

height: int = 100  # px


def headline(color_schema: Color_schema) -> dbc.Row:

    headline = dbc.Row(
        [
            dbc.Col(
                html.Img(src=dlr_img, style={"height": f"{height}px", "width": "auto"}),
                width=2,
                align="left",
            ),
            dbc.Col(
                html.H1(
                    "REMix Dashboard",
                    className="text-center fs-1 fw-bold",
                    style={"color": "blue"},  # color_schema.blue},
                ),
                width=8,
                align="center",
            ),
            dbc.Col(
                html.Img(
                    src=remix_img, style={"height": f"{height}px", "width": "auto"}
                ),
                width=2,
                align="right",
            ),
        ],
        className="mb-3",
    )
    return headline

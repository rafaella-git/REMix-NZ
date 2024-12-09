from dash import html
import dash_bootstrap_components as dbc


def dropdown_headline(hl_text: str, color_schema) -> html.H5:
    return html.H5(hl_text, style={"color": color_schema.blue})


def card_headline(hl_text: str, color_schema) -> dbc.Row:
    return dbc.Row(
        dbc.Col(
            html.H2(
                hl_text,
                className="text-center fw-bold",
                style={"color": color_schema.green},
            ),
            width=12,
            align="center",
        )
    )

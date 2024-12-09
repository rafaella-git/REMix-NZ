import pandas as pd
import plotly.graph_objs as go
from dash import Dash
from plotly.subplots import make_subplots
from dash.dependencies import Output, Input, State

from scr.utils import Dataprep_utils
from scr.color_schema import Color_schema
from scr.functionalities import df_filter, get_tec_col
from layout.tabs.result_tab.converter_tab.converter_functs import group_techs, plot_tech


def group_others(df: pd.DataFrame, per_th: int) -> pd.DataFrame:

    th_df: pd.DataFrame = df.copy()
    th_df["techs"] = th_df.apply(
        lambda row: (
            row["techs"]
            if row["value"] >= df["value"].sum() * per_th / 100
            else "Others"
        ),
        axis=1,
    )
    th_df = th_df.drop(
        columns=["accNodesModel", "accYears", "vintage", "capType"], axis=0
    )
    th_df = th_df.groupby("techs").sum().reset_index()

    return th_df


def converter_unit_callback(
    app: Dash, utils: Dataprep_utils, color_schema: Color_schema
) -> go.Figure:

    @app.callback(
        [
            # Outputs
            Output("conv_unit_plot", "figure"),
            Output("conv_unit_data", "data"),
        ],
        [
            Input("conv_unit_node_dd", "value"),
            Input("conv_unit_captype_dd", "value"),
            Input("conv_unit_group_switch", "on"),
            Input("conv_unit_stack_switch", "on"),
        ],
    )
    def update_conv_unit_fig(node: str, cap_type: str, group: bool, stack: bool):

        barmode: str = "stack" if stack == True else "group"

        converter_units_df: pd.DataFrame = df_filter(
            utils.converter_units_df, node, cap_type
        ).reset_index()

        if group == True:
            converter_units_df = group_techs(converter_units_df)

        tec_col: str = get_tec_col(converter_units_df)

        fig = go.Figure().update_layout(margin=dict(t=35))

        for tec in converter_units_df[tec_col].unique():
            df_tech = converter_units_df[converter_units_df[tec_col] == tec]
            fig.add_trace(plot_tech(df_tech, tec, color_schema))

        fig.update_xaxes(title_text="years")
        fig.update_yaxes(title_text="value / GW")
        fig.update_layout(barmode=barmode, legend_title="Technologies", height=600)

        return fig, [node, cap_type]

    @app.callback(
        Output("conv_unit_pie_plot", "figure"),
        [
            # Inputs
            Input("conv_unit_year_dd", "value"),
            Input("conv_unit_pie_th", "value"),
        ],
        State("conv_unit_data", "data"),
    )
    def update_conv_unit_pie_plot(
        years: str | list[str], per_th: int, data
    ) -> go.Figure:

        converter_units_df = df_filter(utils.converter_units_df, *data).reset_index()

        fig = make_subplots(
            rows=1, cols=len(years), specs=[[{"type": "pie"} for _ in years]]
        ).update_layout(margin=dict(t=35))

        for i, year in enumerate(years):
            col = i + 1

            converter_units_th_df = group_others(
                converter_units_df[converter_units_df["accYears"] == year], per_th
            )

            fig.add_trace(
                go.Pie(
                    labels=converter_units_th_df["techs"],
                    values=converter_units_th_df["value"],
                    # name=year,
                    hole=0.2,
                    title=year,
                ),
                row=1,
                col=col,
            )

        return fig

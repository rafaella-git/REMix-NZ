import pandas as pd
import plotly.graph_objs as go
from dash import Dash
from dash.dependencies import Output, Input
from plotly.subplots import make_subplots

from scr.utils import Dataprep_utils
from scr.color_schema import Color_schema
from scr.functionalities import df_filter, get_tec_col


def group_techs(df: pd.DataFrame) -> pd.DataFrame:

    group_df: pd.DataFrame = df.copy()
    group_df["tech_group"] = [tec.split("_")[0] for tec in df["techs"]]
    group_df: pd.DataFrame = group_df.groupby(
        ["accYears", "tech_group"], as_index=False
    )["value"].sum()

    return group_df


def plot_tech(
    df_tech: pd.DataFrame, tec: str, color_schema: Color_schema
) -> go.Scatter:

    return go.Bar(
        x=df_tech["accYears"],
        y=df_tech["value"],
        name=tec,
        marker_color=color_schema.tech_color_dict[tec],
    )


def pro_con_callback(
    app: Dash, utils: Dataprep_utils, color_schema: Color_schema
) -> go.Figure:

    if utils.commodity_balance_df.empty != True:

        @app.callback(
            Output("pro_con_plot", "figure"),
            [
                Input("pro_con_good_dd", "value"),
                Input("pro_con_node_dd", "value"),
                Input("pro_con_group_switch", "on"),
                Input("pro_con_stack_switch", "on"),
            ],
        )
        def update_pro_con_fig(
            good: str, node: str, group: bool, stack: bool
        ) -> go.Figure:

            barmode: str = "stack" if stack == True else "group"

            # -- data -------------------------------------------------------------------
            bal_type_list: list[str] = (
                utils.generation_annual_df.index.get_level_values(
                    "balanceType"
                ).tolist()
            )
            if "net" in bal_type_list:
                bal_type: str = "net"
            elif "netto" in bal_type_list:
                bal_type: str = "netto"

            generation_df: pd.DataFrame = df_filter(
                utils.generation_annual_df, good, node, bal_type
            ).reset_index()
            demand_df: pd.DataFrame = df_filter(
                utils.demand_annual_df, good, node, bal_type
            ).reset_index()

            if group == True:
                generation_df: pd.DataFrame = group_techs(generation_df)
                demand_df: pd.DataFrame = group_techs(demand_df)

            # -- figure -----------------------------------------------------------------
            fig = make_subplots(rows=1, cols=2, subplot_titles=("Producer", "Consumer"))
            fig.update_layout(margin=dict(t=35))

            tec_col: str = get_tec_col(generation_df)

            for tec in generation_df[tec_col].unique():
                df_tech = generation_df[generation_df[tec_col] == tec]
                fig.add_trace(
                    plot_tech(df_tech, tec, color_schema),
                    row=1,
                    col=1,
                )
            fig.update_xaxes(title_text="years", row=1, col=1)
            fig.update_yaxes(title_text="value / GW", row=1, col=1)

            for tec in demand_df[tec_col].unique():
                df_tech = demand_df[demand_df[tec_col] == tec]
                fig.add_trace(
                    plot_tech(df_tech, tec, color_schema),
                    row=1,
                    col=2,
                )
            fig.update_xaxes(title_text="years", row=1, col=2)

            fig.update_layout(barmode=barmode, legend_title="Technologies", height=800)
            return fig

    else:
        pass

import pandas as pd
import plotly.graph_objs as go

from dash import Dash
from dash.dependencies import Output, Input, State

from scr.utils import Dataprep_utils
from scr.color_schema import Color_schema
from scr.functionalities import df_filter, get_tec_col
from layout.tabs.result_tab.functs import transform_to_heatmap
from remix.framework.tools.gdx import fill_missing_values_in_timeseries


def group_techs(df: pd.DataFrame) -> pd.DataFrame:

    group_df: pd.DataFrame = df.copy()
    group_df["tech_group"] = [tec.split("_")[0] for tec in df["techs"]]
    group_df: pd.DataFrame = group_df.groupby(
        ["timeModel", "tech_group"], as_index=False
    )["value"].sum()

    return group_df


def plot_tech(
    df_filt: pd.DataFrame, tec: str, stackgroup: str, color_schema: Color_schema
) -> go.Scatter:

    return go.Scatter(
        x=df_filt["timeModel"],
        y=df_filt["value"],
        stackgroup=stackgroup,
        name=tec,
        mode="none",
        fillcolor=color_schema.tech_color_dict[tec],
    )


def dispatch_callback(
    app: Dash, utils: Dataprep_utils, color_schema: Color_schema
) -> go.Figure:

    @app.callback(
        [
            Output("dis_plot", "figure"),
            Output("dis_tech_dd", "options"),
            Output("dis_data", "data"),
        ],
        [
            Input("dis_good_dd", "value"),
            Input("dis_node_dd", "value"),
            Input("dis_year_dd", "value"),
            Input("dis_group_switch", "on"),
        ],
    )
    def update_tech_dropdown(good: str, node: str, year: str, group: bool):

        generation_df: pd.DataFrame = fill_missing_values_in_timeseries(
            df_filter(utils.generation_df, year, node, good),
            utils.timemodel_df,
        ).reset_index()
        demand_df: pd.DataFrame = fill_missing_values_in_timeseries(
            df_filter(utils.demand_df, year, node, good),
            utils.timemodel_df,
        ).reset_index()

        if group == True:
            generation_df: pd.DataFrame = group_techs(generation_df)
            demand_df: pd.DataFrame = group_techs(demand_df)

        tec_col_gen: str = get_tec_col(generation_df)
        tec_col_dem: str = get_tec_col(demand_df)
        all_tecs: list[str] = list(
            set(
                list(generation_df[tec_col_gen].unique())
                + list(demand_df[tec_col_dem].unique())
            )
        )

        fig = go.Figure().update_layout(margin=dict(t=35))

        tec_col_gen: str = get_tec_col(generation_df)
        for tec in list(generation_df[tec_col_gen].unique()):
            df_filt: pd.DataFrame = generation_df[generation_df[tec_col_gen] == tec]
            fig.add_trace(plot_tech(df_filt, tec, "positive", color_schema))

        tec_col_dem: str = get_tec_col(demand_df)
        for tec in list(demand_df[tec_col_dem].unique()):
            df_filt: pd.DataFrame = demand_df[demand_df[tec_col_dem] == tec]
            fig.add_trace(plot_tech(df_filt, tec, "negative", color_schema))

        fig.update_xaxes(title_text="time / h")
        fig.update_yaxes(title_text="value / GW")
        fig.update_layout(legend_title="Technologies", height=500)

        return fig, all_tecs, [year, node, good]

    @app.callback(
        Output("dis_hm_plot", "figure"),
        Input("dis_tech_dd", "value"),
        State("dis_data", "data"),
    )
    def update_fig(tech: str, data: tuple) -> go.Figure:

        df_fill: pd.DataFrame = fill_missing_values_in_timeseries(
            df_filter(utils.commodity_balance_df, *data, tech),
            utils.timemodel_df,
        ).reset_index()
        df_hm: pd.DataFrame = transform_to_heatmap(df_fill.copy())

        fig = go.Figure().update_layout(margin=dict(t=35))
        fig.add_trace(
            go.Heatmap(
                z=df_hm.values,
                x=df_hm.columns,
                y=df_hm.index,
                colorscale="Viridis",
                colorbar=dict(thickness=20, title=f"tech. / GW"),
            )
        )
        fig.update_xaxes(title_text="Days", tickmode="linear", dtick=7)
        fig.update_yaxes(title_text="Hours", tickmode="linear", dtick=3)
        fig.update_layout(height=500)

        return fig

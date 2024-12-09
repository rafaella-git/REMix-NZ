import pandas as pd
import plotly.graph_objs as go

from scr.color_schema import Color_schema


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

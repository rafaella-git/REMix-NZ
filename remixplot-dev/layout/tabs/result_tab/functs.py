import pandas as pd


def transform_to_heatmap(df: pd.DataFrame) -> pd.DataFrame:

    df["day"] = ((df["timeModel"] - 1) // 24) + 1
    df["hour"] = ((df["timeModel"] - 1) % 24) + 1

    df_pivot: pd.DataFrame = df.pivot(index="hour", columns="day", values="value")

    all_days: pd.DataFrame = pd.DataFrame(index=range(1, 25), columns=range(1, 366))
    df_hm: pd.DataFrame = all_days.combine_first(df_pivot)

    return df_hm


# weg
def df_filter(df: pd.DataFrame, *args) -> pd.DataFrame:

    for arg in args:
        col: str = df[df.eq(arg)].dropna(axis=1, how="all").columns[0]
        df: pd.DataFrame = df[df[col] == arg]

    return df

import pandas as pd


def restructure_df(df: pd.DataFrame) -> pd.DataFrame:

    idxs: list[str] = df.index.names
    df_res: pd.DataFrame = df.reset_index()

    for col in df_res.columns:
        df_res[col] = df_res[col].to_list()

    df2: pd.DataFrame = df_res.set_index(idxs)

    return df2


def get_links_dict(links: pd.DataFrame) -> tuple[list[str], dict]:

    unique_source_target: list[str] = list(
        pd.unique(links[["source", "target"]].values.ravel("K"))
    )

    mapping_dict: dict = {k: v for v, k in enumerate(unique_source_target)}
    links["source"] = links["source"].map(mapping_dict)
    links["target"] = links["target"].map(mapping_dict)
    links_dict: dict = links.to_dict(orient="list")

    return unique_source_target, links_dict


def get_tec_col(df: pd.DataFrame) -> str:

    if "tech_group" in df.columns:
        tec_col: str = "tech_group"
    elif "techs" in df.columns:
        tec_col: str = "techs"

    return tec_col


def get_col_names(df: pd.DataFrame, *args) -> tuple[list[str], list[str]]:

    col_names: list[str] = []
    idx_names: list[str] = []
    for arg in args:
        key: str = arg[0] if isinstance(arg, list) else arg
        col: str = next(
            name for name, level in zip(df.index.names, df.index.levels) if key in level
        )
        col_names.append(col)
        idx_names.append(arg)

    return col_names, idx_names


def get_idx_names(df: pd.DataFrame, idx: str) -> list[str]:

    idx_names: list[str] = df.index.get_level_values(idx).unique().tolist()

    return idx_names


def df_filter_index_slicing(df: pd.DataFrame, **kwargs) -> pd.DataFrame:

    names: list[str] = df.index.names
    lst: list[slice] = [
        kwargs[name] if name in kwargs else slice(None) for name in names
    ]
    idx_lst: list[slice] = [[l] if not isinstance(l, slice | list) else l for l in lst]

    df_filt: pd.DataFrame = df.loc[tuple(idx_lst)]
    return df_filt


def df_filter(df: pd.DataFrame, *args) -> pd.DataFrame:

    for i, arg in enumerate(args):
        if i == 0:
            col_names, idx_names = get_col_names(df, arg)
            df_filt: pd.DataFrame = df_filter_index_slicing(
                df,
                **{key: val for key, val in zip(col_names, idx_names)},
            )
        else:
            col_names, idx_names = get_col_names(df_filt, arg)
            df_filt: pd.DataFrame = df_filter_index_slicing(
                df_filt,
                **{key: val for key, val in zip(col_names, idx_names)},
            )

    return df_filt

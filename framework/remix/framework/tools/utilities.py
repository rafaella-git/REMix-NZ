# collection of functions to evaluate REMix results
import os
from collections.abc import Mapping
from copy import deepcopy
from io import StringIO
from pathlib import Path
from types import FunctionType
from typing import Union

import numpy as np
import pandas as pd

PROFILEIDS = ["profile", "timeseries"]


def merge_dfs(df_list: list):
    """Merge list of pandas.DataFrames to one common DataFrame.

    Parameters
    ----------
    df_list : list of pandas.DataFrames
        list of pandas.DataFrames.

    Returns
    -------
    pandas.DataFrame
        single pandas.DataFrame.
    """
    if len(df_list) == 1:
        return df_list[0]

    elif len(df_list) > 0:
        all_columns = []
        for df in df_list:
            all_columns += df.columns.tolist()
        all_columns = list(set(all_columns))

        df_list_out = [df for df in df_list if (not df.empty)]

        if len(df_list_out) == 1:
            df_out = df_list_out[0]

        elif len(df_list_out) == 0:
            new_index = df_list[0].index
            for df in df_list[1:]:
                new_index = new_index.union(df.index)

            new_index.names = df_list[0].index.names
            return pd.DataFrame(index=new_index)

        else:
            df_out = pd.concat(df_list_out, axis=0)

        for missing in set(all_columns) - set(df_out.columns):
            df[missing] = np.nan

        df_out.index.names = df_list[0].index.names
        df_out = df_out.groupby(level=list(range(df_out.index.nlevels))).agg("last")
        df_out = df_out.sort_index()

    else:
        df_out = pd.DataFrame()

    return df_out


def merge_dicts(dict1, dict2, update_only=False):
    """Merge the contents of two dicts.

    Parameters
    ----------
    dict1 : dict
        Dictionary to be updated with dict2.

    dict2 : dict
        Dictionary to merge into dict1.

    update_only : bool
        Only update values of key existing in dict1.

    Returns
    -------
    dict
        Merged dictionary.

    Example
    -------
    Run the function with updating and appending:

    >>> from remix.framework.tools.utilities import merge_dicts
    >>> dict1 = {"a": {"b": 2}, "c": 1}
    >>> dict2 = {"a": {"b": 3, "c": 5}, "d": {"e": 0}}
    >>> result = merge_dicts(dict1, dict2)
    >>> result["a"]["b"]
    3
    >>> result["a"]["c"]
    5
    >>> result["d"]
    {'e': 0}

    Allowing updating only will throw an error:

    >>> result = merge_dicts(dict1, dict2, update_only=True)
    Traceback (most recent call last):
      ...
    AttributeError: The key "c" of dict2 is not in the current level keys of dict1. If you call this method with update_only=True you cannot add new keys to dict1. Available keys of dict1 at this level are "b".
    """
    result = deepcopy(dict1)

    for key, value in dict2.items():
        if update_only and (key not in result):
            msg = (
                f'The key "{key}" of dict2 is not in the current level keys of dict1. If you call this method '
                "with update_only=True you cannot add new keys to dict1. Available keys of dict1 at this level "
                'are "' + '", "'.join(result.keys()) + '".'
            )
            raise AttributeError(msg)
        if isinstance(value, Mapping):
            result[key] = merge_dicts(
                result.get(key, {}), value, update_only=update_only
            )
        else:
            result[key] = deepcopy(dict2[key])

    return result


def concat_ts(df_list):
    """Concatenate profiles from list of pandas.DataFrames.

    Parameters
    ----------
    df_list : list
        List of pandas.DataFrames.

    Returns
    -------
    pandas.DataFrame
        Single pandas.DataFrame.
    """
    df_out = pd.DataFrame()

    if len(df_list) > 0:
        idx = pd.MultiIndex.from_tuples(np.concatenate([i.index for i in df_list]))
        col = df_list[0].columns
        data = np.concatenate([i.to_numpy() for i in df_list])
        df_out = pd.DataFrame(data=data, index=idx, columns=col).dropna().astype(float)

    return df_out


def read_dat(file: str):
    """Read :code:`.dat` files into pandas.DataFrames.

    Parameters
    ----------
    file : str
        Name of the :code:`.dat` file.

    Returns
    -------
    pandas.DataFrame : pandas.DataFrame
        Same content as the input code:`.dat` file.
    """
    if os.stat(file).st_size == 0:
        return None

    skiplines = []
    with open(file) as f:
        skiplines.extend(ln for ln, line in enumerate(f) if line.startswith("*"))
    df = pd.read_csv(file, sep=r"\s+", skiprows=skiplines)

    if isinstance(df.index, pd.RangeIndex):
        df = pd.read_csv(file, sep=r"\s+", skiprows=skiplines, header=None, dtype=str)
        df = df.set_index(list(df.columns))
        if type(df.index) != pd.MultiIndex:
            df.index = df.index.astype(str)

    if isinstance(df.index, pd.MultiIndex):
        drop_index_levels = list(range(1, len(df.index.levels), 2))
        df.index = df.index.droplevel(drop_index_levels)
        if isinstance(df.index, pd.MultiIndex):
            df.index.names = [None] * len(df.index.levels)

            for i in range(len(df.index.levels)):
                idx = df.index.levels[i].astype(str)
                df.index = df.index.set_levels(idx, level=i)

    elif len(df.index) == 0 and len(df.columns) == 3:
        # special case of a map with one item
        df = pd.read_csv(file, sep=r"\s+", skiprows=skiplines, header=None)
        df = df.set_index([0, 2])
        df.index.names = [None, None]
        df = pd.DataFrame(index=df.index)

    if (len(df.index) == 0) and (len(df.columns) == 1):
        # special case of a list with one item
        df = pd.read_csv(file, header=None)
        df[0] = df[0].astype(str)
        df = df.set_index(0)

    if (len(df.index) == 0) and (len(df.columns) == 2):
        # special case of a list with one item and comments
        df = pd.read_csv(file, header=None)
        df[0] = df[0].astype(str).str.split('"')
        df[0] = df[0][0][0].strip()
        df = df.set_index(0)

    return df


def to_dat(df: pd.DataFrame, fname: str, float_format=None):
    """Transform pandas.DataFrames to :code:`.dat` files.

    Parameters
    ----------
    df : pandas.DataFrame
        pandas.DataFrame.
    fname : str
        File name; path to file the pandas.DataFrame is supposed to be written to.
    """
    Path(fname).parent.mkdir(parents=True, exist_ok=True)
    df_out = df.copy()

    if isinstance(df_out.index, pd.MultiIndex):
        idx_gms = pd.Index(
            [" . ".join(map(str, df_out.iloc[i].name)) for i in range(len(df_out))]
        )
        df_out.set_index(idx_gms, inplace=True)

    contents = ""
    # if precision is not None:
    #     float_string = '{:.' + str(precision) + 'f}'
    #     float_format = float_string.format
    # else:
    #     float_format = None
    if not df.empty:
        contents = df_out.to_string(
            index_names=False, sparsify=False, float_format=float_format
        )
    else:
        if len(df.index) > 0 and len(df.columns) == 0:
            setlist = df.index.tolist()
        elif len(df.index) == 0 and len(df.columns) > 0:
            setlist = df.columns.tolist()
        else:
            return

        if isinstance(setlist[0], tuple):
            setlist = [" . ".join(i) for i in setlist]
        for element in setlist:
            contents += str(element) + "\n"

    with open(fname, "w") as outData:
        outData.write(contents)


# monkey patching intended here?
pd.DataFrame.to_dat = to_dat


def get_limited_df(df: pd.DataFrame, file: str, _logger: FunctionType = lambda x: None):
    """Limit a DataFrame to 80,000 lines.

    Parameters
    ----------
    df : pandas.DataFrame
        Original DataFrame from raw input data.
    file : str
        Path of the input file
    _logger : FunctionType, optional
        A logging function that takes strings as input, by default lambda x: None.

    Returns
    -------
    pandas.DataFrame
        Rounded data conform with the 80,000 lines limit of GAMS.
    """
    for digits in [None, 4, 3, 2, 1, 0]:
        textbuf = StringIO()
        if digits is None:
            df_new = df
            df_new.to_csv(textbuf)
        else:
            df_new = df.round(digits)
            df_new.to_csv(textbuf)
        lmax = max([len(l) for l in textbuf.getvalue().splitlines()])
        if lmax < 80e3:
            if _logger is not None:
                _logger(
                    "Maximum line length {} for file {} - continuing".format(lmax, file)
                )
            return df_new
        else:
            if _logger is not None:
                _logger(
                    "Maximum line length {} for file {} - automatic rounding to fewer digits".format(
                        lmax, file
                    )
                )
    return df_new


def read_remix_csv(file, schema=None):
    """Read input data with no column names in the csv inputs.

        This function expects that the index columns are empty.

    Parameters
    ----------
    file : str
        Path to the file.

    Returns
    -------
    pandas.DataFrame
        pandas.DataFrame with the data stored from the .csv file.
    """
    if os.stat(file).st_size == 0:
        return None
    df = pd.read_csv(file)
    if len(df.columns) == 1:
        df = pd.read_csv(file, header=None)
        df[0] = df[0].astype(str)
        df = df.set_index(0)
        if all(" . " in s for s in list(df.index)):
            # this seems to be a map
            df.index = pd.MultiIndex.from_tuples(
                [tuple(s.split(" . ")) for s in df.index]
            )

    else:
        if schema is None:
            index_cols = [col for col in df.columns if "Unnamed" in col]
            index_length = len(index_cols)
            # Turn the indeces into strings (?) performance penalty?
            df = df.set_index(index_cols)
            df.index.names = [None] * index_length

        elif isinstance(schema, dict):
            index_columns = [fk["fields"][0] for fk in schema["foreignKeys"]]
            try:
                df = df.set_index(index_columns)
            except KeyError:
                index_columns = [
                    ix for ix in index_columns if ix.lower() not in ["timedata"]
                ]
                df = df.set_index(index_columns)

        else:
            msg = "The provided schema must be a dictionary."
            raise TypeError(msg)

        for i in range(len(df.index.levels)):
            idx = df.index.levels[i].astype(str)
            df.index = df.index.set_levels(idx, level=i)
    return df


def read_file(file):
    """Read input data for a DataFrame from a .csv or .dat file.

    Parameters
    ----------
    file : str
        Path to the file.

    Returns
    -------
    pandas.DataFrame
        pandas.DataFrame with the data stored from the file.
    """
    if file.endswith(".csv"):
        df = read_remix_csv(file)
    elif file.endswith(".dat"):
        df = read_dat(file)
    else:
        raise TypeError("File extension not valid")
    return df


def build_timeseries_range(as_list: bool = False):
    """Build a list for the time steps of the REMix model.

    Parameters
    ----------
    as_list : bool
        Get the time series as list, by default False.

    Returns
    -------
    list
        Timesteps from 1 through 8760.
    """
    return (
        [[f"t{str(n).zfill(4)}"] for n in range(1, 8761)]
        if as_list
        else [f"t{str(n).zfill(4)}" for n in range(1, 8761)]
    )


def is_empty(value: Union[pd.DataFrame, list]):
    """Empty list OR DataFrame.

    Parameters
    ----------
    value: pandas.DataFrame, list
        List or pandas.DataFrame.

    Returns
    -------
    bool
        True if the list or pandas.DataFrame is empty.
    """
    if isinstance(value, pd.DataFrame):
        if value.empty and isinstance(value.index, pd.MultiIndex):
            return all(len(i) == 0 for i in value.index.values)
        return value.empty
    if isinstance(value, list):
        return len(value) == 0


def is_timeseries(name: str) -> bool:
    """Check if the name corresponds to a time series.

    Parameters
    ----------
    name : str
        Name of the parameter to be checked.

    Returns
    -------
    bool
        True if the name corresponds to a time series.
    """
    return any(identifier in name.lower() for identifier in PROFILEIDS)


def is_map(name: str) -> bool:
    """Check if the name corresponds to a map.

    Parameters
    ----------
    name : str
        Name of the parameter to be checked.

    Returns
    -------
    bool
        True if the name corresponds to a map.
    """
    return "map_" in name.lower()

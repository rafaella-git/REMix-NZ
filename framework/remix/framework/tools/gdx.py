import os
from multiprocessing import Pool
import re
from typing import Optional, Union

import pandas as pd


class GAMSVersionError(Exception):
    """Exeption for newer GAMS Versions."""

    pass


try:
    import gamstransfer as gt
except ModuleNotFoundError:
    gt = None
    from remix.framework.api.run import _find_gams_pythonapi

    version = _find_gams_pythonapi()

    if 37 <= version < 42:
        import gamstransfer as gt
    elif version < 37:
        raise GAMSVersionError(
            f'The REMix Python API is not compatible with GAMS version {version} or anything below 37, \
            you can still run your models using "gams path/to/remix/framework/model/run_remix.gms" \
            but you might get unexpected problems.'
        )
    else:
        try:
            import gams.transfer as gt
        except ModuleNotFoundError:
            raise GAMSVersionError(
                f"It seems that you are using a GAMS version {version} which is higher than 41, \
                please install the API as instructed in \
                https://www.gams.com/{version}/docs/API_PY_GETTING_STARTED.html#PY_PIP_INSTALL_BDIST"
            )


symbols_without_time = [
    "indicator_accounting",
    "indicator_accounting_comp",
    "indicator_accounting_detailed",
    "commodity_balance_annual",
    "transfer_flows_annual",
    "storage_flows_annual",
    "converter_caps",
    "transfer_caps",
    "storage_caps",
]


def read_gdx_symbol(gdx_path: Union[str, dict], gdx_symbol: str, sel_dict=None, fill_ts=True) -> pd.DataFrame:
    """Read in gdx files into pandas.DataFrames with GAMS Transfer.

    Parameters
    ----------
    gdx_path : str
        The path to the gdx file holding the result data.
    gdx_symbol : str
        String of gdx_symbol in REMix model output, e.g. "converter_caps", "commodity_balance", ...
    sel_dict : dict
        Optional dictionary to use for data selection.
    fill_ts: bool
        Try to automatically fill gaps in the time index created by GAMS when squeezing zeros.

    Returns
    -------
    pandas.DataFrame
        DataFrame with the result data from the gdx_symbol.
    """
    cont_meta = gt.Container()
    cont_meta.read(gdx_path, records=False)
    time_symbols = [t for t in ["timeModel", "timeModelToCalc"] if t in cont_meta.listSymbols()]

    cont = gt.Container()
    cont.read(gdx_path, [gdx_symbol] + time_symbols)
    df = cont.data[gdx_symbol].records

    if df is None:
        return pd.DataFrame()

    new_cols = {col: re.sub(r"\_(\d+)$", "", col) for col in df.columns if "_" in col}
    if len(new_cols.values()) != len(set(new_cols.values())):
        unique_cols = {}
        counter = {}
        for k, v in new_cols.items():
            if v in unique_cols.values():
                if v in counter:
                    counter[v] += 1
                else:
                    counter[v] = 1
                counter = counter[v]
                unique_cols[k] = f"{v}_{counter}"
            else:
                unique_cols[k] = v
    else:
        unique_cols = new_cols

    df = df.rename(columns=unique_cols)
    df = df.set_index(df.columns.to_list()[:-1])

    if isinstance(sel_dict, dict):
        df = select_from_dataframe(df, sel_dict)

    if fill_ts and ("timeModel" in df.index.names) and isinstance(df.index, pd.MultiIndex):
        if "timeModelToCalc" in cont.listSymbols():
            timetocalc_records = cont.data["timeModelToCalc"].records
            timetocalc_records.columns = [re.sub("_\d+$", "", c) for c in timetocalc_records.columns]
            df = fill_missing_values_in_timeseries(df, timetocalc_records.set_index("timeModel"))
        elif "timeModel" in cont.listSymbols():
            timemodel_records = cont.data["timeModel"].records
            timemodel_records.columns = [re.sub("_\d+$", "", c) for c in timemodel_records.columns]
            df = fill_missing_values_in_timeseries(df, timemodel_records.set_index("timeModel"))
        else:
            df = fill_missing_values_in_timeseries(df, pd.DataFrame(index=df.index.get_level_values("timeModel").unique()))

    return df


def select_from_dataframe(df: pd.DataFrame, sel_dict: dict) -> pd.DataFrame:
    """Filter a dataframe via dictionary for each index (key) using a string or list of allowed entries (value).

    Args:
        df (pd.DataFrame): Dataframe to select from.
        sel_dict (dict): Dictionary with index value pairs to select by.

    Returns:
        pd.DataFrame: Filtered dataframe
    """
    for idx, sel in sel_dict.items():
        if idx in df.index.names:
            _indexer = [slice(None)] * len(df.index.levels)
            idx_pos = df.index.names.index(idx)
            _indexer[idx_pos] = sel
            indexer = tuple(_indexer)
            df = df.loc[indexer, slice(None)]
    return df


def write_gdx_symbol_as_csv(gdx: str, symbol: str, directory: str):
    """GDX symbol translation function

    Args:
        gdx (str): Name of gdx file
        symbol (str): Symbol in the gdx file
        directory (str): Main directory to save gdx file
    """
    df = read_gdx_symbol(gdx, symbol)
    gams_43 = all(
        title in df.reset_index().columns for title in ["uni", "element_text"]
    )
    gams_older = all(title in df.reset_index().columns for title in ["uni", "element"])
    if gams_older or gams_43:
        output = df.index.to_list()
        if len(output) > 0:
            with open(
                os.path.join(directory, f"set_{symbol}.csv"), "w", encoding="UTF-8"
            ) as out_data:
                for element in output:
                    out_data.write(str(element) + "\n")
    elif "element" in df.columns:
        df = df.drop(columns=["element"])
        df.to_csv(os.path.join(directory, f"{symbol}.csv"), index=True)
    else:
        if not df.empty:
            df.to_csv(os.path.join(directory, f"{symbol}.csv"))


def fill_missing_values_in_timeseries(df: pd.DataFrame, timemodel: pd.DataFrame) -> pd.DataFrame:
    """Fill missing values (0.0) in timeseries.

    GAMS truncates the result data where a value is equal to 0. In time series plotting, this is disadvantagous,
    because depending on the commodity time series have different lenght.

    Parameters
    ----------
    df : pandas.DataFrame
        Pandas DataFrame containing time series result data, e.g. commodity_balance, or transfer_flows.
    timemodel : pandas.DataFrame
        Pandas DataFrame with either the timeModel or preferrably the timeModelToCalc.

    Returns
    -------
    pandas.DataFrame
        The timeseries with missing values added as 0.0.

    Raises
    ------
    ValueError
        If the provided input dataframe does not have a time dimension.
    """
    if "timeModel" not in df.index.names:
        msg = "The data provided does not have a timeModel dimension."
        raise ValueError(msg)

    lvl_pre = list(range(len(df.index.names))[: df.index.names.index("timeModel")])
    lvl_time = [df.index.names.index("timeModel")]
    lvl_post = list(range(len(df.index.names))[df.index.names.index("timeModel") :][1:])
    lvl_reorder = lvl_time + lvl_pre + lvl_post

    mi = pd.MultiIndex.from_product(
        [timemodel.index.to_list(), df.index.droplevel(lvl_time).unique().to_list()]
    )
    mi_time = mi.get_level_values(0)
    mi_sub = pd.MultiIndex.from_tuples(mi.get_level_values(1))

    idx_new = pd.MultiIndex.from_arrays(
        [mi_time] + [mi_sub.get_level_values(i) for i in range(len(mi_sub.levels))]
    )
    idx_new = idx_new.reorder_levels(lvl_reorder)
    idx_new.names = df.index.names
    df = df.reindex(idx_new).fillna(0).sort_index()

    return df


class GDXEval:
    def __init__(
        self,
        gdx: Union[str, dict],
        fill_ts: bool = True,
        lazy: bool = True,
        parallel: bool = True,
        processes: Optional[int] = None,
        scen_names: list[str] = ["scenario"],
    ) -> None:
        """Initialise attributes for evaluation of gdx files

        Parameters
        ----------
        gdx : str | dict
            Path to gdx file or dictionary with gdx paths as values.
        fill_ts : bool, optional
            Specify if timeseries with missing intervals should be filled in
            with zeros during loading. Defaults to True.
        lazy : bool, optional
            Specify if loading of symbols should be done lazily. Defaults to True.
        parallel : bool, optional
            Specify if multiprocessing should be used when loading from multiple
            gdx files. Defaults to True.
        processes : int, optional
            Number of parallel processes to use for parallel reading. Defaults
            to None (=CPU count).
        scen_names : list, optional
            List of scenario names for multi-index keys. Defaults to ["scenario"].
        """
        self._fill_ts = fill_ts
        self._lazy = lazy
        self._parallel = parallel
        self._processes = processes
        self._gdx_dict = isinstance(gdx, dict)
        self._gdx = self.filter_scens(gdx)
        self._symbols = self.get_metadata(self._gdx)
        self._scen_names = scen_names
        self._dfs = {}
        if not self._lazy:
            c = gt.Container()
            c.read(gdx)
            for k, p in c.data.items():
                self._dfs[k] = p.records

    def __repr__(self) -> str:
        """Prints the list of available symbols.

        Returns
        -------
        str
            message listing all available symbols
        """
        return "Available symbols:\n" + "\n".join(
            [f"{m}" for m in sorted(self._symbols)]
        )

    def __getitem__(self, key: str) -> pd.DataFrame:
        """Returns the cached DataFrame or loads it from the gdx file in case of
        lazy loading.

        Parameters
        ----------
        key : str
            Symbol name from the gdx file.

        Returns
        -------
        pd.DataFrame
            DataFrame containing the results.
        """
        if key not in self._dfs.keys():
            self._dfs[key] = self.get_symbol(key, scen_names=self._scen_names)
        if isinstance(self._dfs[key].index, pd.MultiIndex):
            if ("timeModel" in self._dfs[key].index.names) & (
                isinstance(self._dfs[key].index.levels[0].values[0], str)
            ):
                self._dfs[key].index = self._dfs[key].index.set_levels(
                    self._dfs[key].index.levels[0].str.strip("tm").astype(int),
                    level="timeModel",
                )

                self._dfs[key] = self._dfs[key].sort_index()

        return self._dfs[key]

    def filter_scens(self, gdx: Union[str, dict]) -> Union[str, dict]:
        """Checks if the path to the gdx file exists or filters the gdx paths in the case of a dictionary.

        Parameters
        ----------
        gdx : str | dict
            Path to gdx file or dictionary with gdx paths as values.

        Returns
        -------
        str | dict
            Verified path or dictionary with valid gdx paths.
        """
        if self._gdx_dict:
            drop_scens = [scen for scen, pth in gdx.items() if not os.path.isfile(pth)]
            for k in drop_scens:
                print(f"File {gdx[k]} does not exist, skipping.")
                gdx.pop(k)
        elif not os.path.isfile(gdx):
            raise IOError(f"File {gdx} does not exist, aborting.")
        return gdx

    def get_metadata(self, gdx: Union[str, dict]) -> list:
        """Open the gdx file to extract symbol list from the metadata.

        Parameters
        ----------
        gdx : str | dict
            Path to gdx file or dictionary with gdx paths as values.

        Returns
        -------
        list
            List of available symbols in the gdx file.
        """
        c = gt.Container()
        if self._gdx_dict:
            c.read(list(gdx.values())[0], records=False)
        else:
            c.read(gdx, records=False)
        return list(c.data)

    def get_symbol(
        self, key: str, scen_names: list[str], sel_dict: Optional[dict] = None
    ) -> pd.DataFrame:
        """Provide the DataFrame based on the symbol name.
        An additional index is prepended in case of a gdx dictionary.

        Parameters
        ----------
        key : str
            Symbol name from the gdx file.
        scen_names : list[str]
            names of scenarios
        sel_dict : dict
            Key/value pairs corresponding to index names and values to select by.

        Returns
        -------
        pd.DataFrame
            Dataframe containing the results.
        """
        if self._gdx_dict:
            # Build tuple of paths and symbols to be loaded
            arg_list = tuple((pth, key, sel_dict) for pth in list(self._gdx.values()))

            # Parallel read-in of partial gdx files
            if self._parallel:
                with Pool(processes=self._processes) as p:
                    df_list = p.starmap(read_gdx_symbol, arg_list)
            else:
                df_list = [
                    read_gdx_symbol(pth, key, sel_dict)
                    for pth, key, sel_dict in arg_list
                ]

            # Concat dataframes with additional scenario index
            df_list = [d for d in df_list if not d.empty]
            df = pd.concat(df_list, keys=self._gdx.keys())
            idxn = list(df.index.names)
            for i, j in enumerate(scen_names):
                idxn[i] = j
            df.index.names = idxn
        else:
            df = read_gdx_symbol(self._gdx, key)

        return df

    def export(
        self,
        directory: str = ".",
        keys: list = [],
        format: str = "csv",
        sel_dict: Optional[dict] = None,
        stacked: bool = False,
    ):
        """Export gdx files to the given format.

        Parameters
        ----------
        directory : str, optional
            Path to the directory. Defaults to "." .
        keys : list, optional
            list of symbols to export. Defaults to [].
        format : str, optional
            Format of the output file. Defaults to "csv".
        sel_dict : dict
            Key/value pairs corresponding to index names and values to select by.
        stacked : bool
            If true, values will be aggregated.
        """
        if not os.path.exists(directory):
            os.mkdir(directory)
        if len(keys) == 0:
            keys = self._symbols
        if self._gdx_dict:
            if not stacked:
                for name, gdx in self._gdx.items():
                    for k in keys:
                        outdir = os.path.join(directory, name)
                        if not os.path.exists(outdir):
                            os.mkdir(outdir)
                        write_gdx_symbol_as_csv(gdx, k, outdir)
            else:
                raise NotImplementedError("This feature is in development")
                # for k in keys:
                #     symbol = self.get_symbol(k,sel_dict=sel_dict)
                #     symbol.to_csv(os.path.join(directory, f"{k}.csv"))
        else:
            if format == "csv":
                for k in keys:
                    write_gdx_symbol_as_csv(self._gdx, k, directory)
            else:
                raise IOError(f"Format {format} is invalid, aborting.")

    def to_xlsx(self, filename: str = "results.xlsx", sel_dict: Optional[dict] = None):
        """Export gdx symbols to xlsx file.

        Parameters
        ----------
        filename : str, optional
            Name of *.xlsx file to be created. Defaults to "results.xlsx".
        """
        with pd.ExcelWriter(filename, engine="openpyxl") as writer:
            for i in symbols_without_time:
                if i in self._symbols:
                    df = self[i]
                    if sel_dict is not None:
                        df = select_from_dataframe(df, sel_dict)
                    df.to_excel(writer, sheet_name=i, merge_cells=False)

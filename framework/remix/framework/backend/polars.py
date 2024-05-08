import re
from io import StringIO
from os.path import join
from pathlib import Path

from remix.framework.api.dataset import DirectoryDataset
from remix.framework.schema.templates import sm
from remix.framework.tools.utilities import build_timeseries_range

try:
    import polars as pl
except ImportError as e:
    pl = None


class PolarsDataset(DirectoryDataset):
    """Polars based REMix directory dataset.

    Parameters
    ----------
    dataframes: dict, optional
        DataFrames contained in the object, by default None.
    files : list, optional
        List of files contained in the object, by default empty list.
    name : str, optional
        Name of the scenario, by default ".".
    parent_files : list, optional
        Files of parent scenario, by default empty list.
    path : str, optional
        Path where the item's input files are located, by default None.
    roundts : int, optional
        If 1, the output functions will round time-series values, by default 0.
    version : str, optional
        REMix model version, by default "latest".
    """

    def __init__(self, **kwargs):
        if pl is None:
            raise ImportError("Please install the polars library to use this feature. Tip: use pip install remix.framework[polars]")
        super().__init__(**kwargs)

    def _get_file_schema(self, path=None):
        return self._schemas[Path(path).stem.lower()]

    def _set_dataframe_to_list(self, dataframe):
        return dataframe["set"].to_list()

    def _merge_dfs(self, df_list, file_schema=None):
        """Function to merge list of polars DataFrames to one common DataFrame.

        Parameters
        ----------
        df_list : list
            list of polars DataFrames.
        file_schema : dict
            Schema.

        Returns
        -------
        polars.DataFrame
            single polars DataFrame.
        """
        all_columns = [col["name"] for col in file_schema["fields"]]
        index_columns = [fk["fields"][0] for fk in file_schema["foreignKeys"]]
        value_columns = [col for col in all_columns if col not in index_columns]
        if len(df_list) > 0:
            df_out = pl.concat(df_list)
            value_columns = [pl.col(vc).last() for vc in value_columns if vc in df_out.columns]
            if "Value" in df_out.columns:
                df_out = df_out.pivot(values="Value", index=index_columns[:-1], columns=index_columns[-1], aggregate_function="first")
            if len(value_columns) != 0:
                if "t0001" not in df_out.columns:
                    df_out.group_by(by=index_columns).agg(value_columns)
                else:
                    index_columns = [ix for ix in index_columns if ix.lower() not in ["timedata"]]
                    value_columns = [pl.col(vc).last() for vc in build_timeseries_range()]
                    df_out.group_by(by=index_columns).agg(value_columns)
                df_out = df_out.fill_nan(0)
        else:
            df_out = pl.DataFrame()

        return df_out

    def _get_limited_df(self, df, file, _logger=None, file_schema=None):
        """Limit a DataFrame to 80,000 lines.

        Parameters
        ----------
        df : polars.DataFrame
            Original DataFrame from raw input data.
        file : str
            Path to the input file.
        _logger : function
            Logging function. Defaults to None.
        _logger : dict
            Dictionary with the shema associated with the input file. Defaults to None.

        Returns:
            polars.DataFrame
                Rounded data conform with the 80,000 lines limit of GAMS.
        """
        # all_columns = [col["name"] for col in file_schema["fields"]]
        index_columns = [fk["fields"][0] for fk in file_schema["foreignKeys"] if fk["fields"][0] .lower() not in ["timedata"]]
        for digits in [None, 4, 3, 2, 1, 0]:
            if digits is None:
                df_new = df
                textbuf = StringIO(df_new.write_csv())
            else:
                value_columns = [pl.col(vc).round(decimals=digits).last() for vc in build_timeseries_range()]
                df_new = df.groupby(index_columns).agg(value_columns)
                textbuf = StringIO(df_new.write_csv())
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

    def write_map_file(self, data, parent, name, extension, **kwargs):
        """Write a map file.

        Parameters
        ----------
        data : polars.Dataframe, list
            Data to write.
        parent : str
            Path to the parent directory.
        name : str
            Name of the map file
        extension : str
            Datatype extension, e.g. .csv, .dat.
        map_format : str
            Format of the map files to be used, use csv for a better shareable file. Defaults to dat.
        """
        filename = name + extension
        map_format = kwargs.get("map_format", "dat")
        if map_format == "dat":
            data = data.with_columns(pl.lit(" . ").alias("middle"))
            data = data.select(pl.concat_str([data.columns[0], data.columns[2], data.columns[1]]))
            data.columns = [""]
        data.write_csv(join(parent, filename), include_header=True)

    def write_parameter_file(self, data, parent, name, extension, **kwargs):
        """Write a parameter file.

        Parameters
        ----------
        data : polars.Dataframe, list
            Data to write.
        parent : str
            Path to the parent directory.
        name : str
            Name of the parameter file
        extension : str
            Datatype extension, e.g. .csv, .dat.
        """
        filename = name + extension
        self.logger(f"Inheriting data for file {name}")
        data.write_csv(join(parent, filename))

    def write_profile_file(self, data, parent, name, extension, **kwargs):
        """Write a profile file.

        Parameters
        ----------
        data : polars.Dataframe, list
            Data to write.
        parent : str
            Path to the parent directory.
        name : str
            Name of the profile file
        extension : str
            Datatype extension, e.g. .csv, .dat.
        """
        profile_format = kwargs.get("profile_format", "wide")
        if profile_format == "tall":
            index_columns = [c for c in data.columns if not re.match("t\d\d\d\d", c)]
            value_cols = [c for c in data.columns if re.match("t\d\d\d\d", c)]
            data = data.melt(id_vars=index_columns, value_vars=value_cols, value_name="Value", variable_name="timeData")
        self.write_parameter_file(data, parent, name, extension, **kwargs)

    def read_file(self, _input: str):
        """This interface was added to have flexible input file reading."""
        if self._schemas is None:
            self._schemas = {k.lower(): v for k, v in sm.get_schemas().items()}
        if "map_" in _input:
            dataframe = pl.read_csv(_input, has_header=True)
            file_schema = self._schemas[Path(_input).stem.lower()]
            index_columns = [fk["fields"][0] for fk in file_schema["foreignKeys"] if fk["fields"][0]]
            if not (set([id.lower() for id in index_columns]) == set([id.lower() for id in dataframe.columns])):
                dataframe.columns = ["map"]
                dataframe = dataframe.with_columns([
                    pl.col('map'),
                    *[pl.col('map').map_elements(lambda s, i=i: s.split(' . ')[i]).alias(col_name)
                        for i, col_name in enumerate(index_columns)]
                ])
                dataframe = dataframe.select(index_columns)
        else:
            if "set_" in _input:
                has_header = False
                columns = ["set"]
            else:
                has_header = True
                columns = None
            dataframe = pl.read_csv(_input, has_header=has_header, new_columns=columns)
        return dataframe

    def write(self, output_path=None, inspect=False, map_format="dat", profile_format="wide", **kwargs):
        """Main inheritance method. Converts DataFrames into .csv files.

        Parameters
        ----------
        output_path : str, optional
            Directory where the scenario will be written to, by default None.
        inspect : bool, optional
            If true, this function will return the DataFrames, by default False.
        map_format : str, optional
            The shape of the map file, using "dat" will produce a GAMS-compatible map, by default "dat".
        profile_format : str, optional
            If wide, the time indeces will be used as columns, if tall, there will be only one value column, by default
            "wide".
        support_sets : bool, optional
            If true, structural validation sets will be part of the output, by default False.

        Returns
        -------
        polars.dataframe
            In case inspect is True, else returns NoneType.
        """
        dataframes = super().write(output_path=output_path, inspect=inspect, map_format=map_format, profile_format=profile_format, **kwargs)
        return dataframes

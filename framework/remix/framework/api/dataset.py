# -*- coding: utf-8 -*-
"""REMix datasets.

This module contains class definitons of objects used to handle REMix instances
input data.

Todo:
    * Implement caching to harmonize with inheritance.
"""
from __future__ import annotations

import glob
import os
from os.path import exists
from os.path import join
from pathlib import Path
from types import FunctionType
from typing import Union

from pandas import DataFrame

from remix.framework.tools.utilities import build_timeseries_range
from remix.framework.tools.utilities import get_limited_df
from remix.framework.tools.utilities import is_map
from remix.framework.tools.utilities import is_timeseries
from remix.framework.tools.utilities import merge_dfs
from remix.framework.tools.utilities import read_file
from remix.framework.tools.utilities import read_remix_csv
from remix.framework.tools.utilities import to_dat


class Dataset:
    """Dataset interface of the REMix scenario dataset.

    This class is the interface for REMix dataset implementations.

    Implements the write method that will be shared across implementations.
    The back end consists of Pandas DataFrames.

    Parameters
    ----------
    dataframes: dict, optional
        Dataframes contained in the object, by default None.
    files : list, optional
        List of files contained in the object, by default empty list.
    name : str, optional
        Name of the scenario, by default ".".
    parent_files : list, optional
        Files of parent scenario, by default empty list.
    path : str, optional
        Path where the items input files are located, by default None.
    roundts : int, optional
        If 1, the output functions will round time series values, by default 0.
    version : str, optional
        REMix model version, by default "latest".
    """

    def __init__(self, **kwargs) -> None:

        self.dataframes = kwargs.get("dataframes", None)
        self.files = kwargs.get("files", [])
        self.name = kwargs.get("name", ".")
        self.parent_files = kwargs.get("parent_files", [])
        self.path = kwargs.get("path", None)
        self.version = kwargs.get("version", "latest")

        self.logger = lambda x: None
        self.ready = False

    def write(self, output_path=None, inspect=False, fileformat="csv", map_format="dat", profile_format="wide", support_sets=False, **kwargs):
        """Main inheritance method. Converts dataframes into .csv or .dat files.

        Parameters
        ----------
        output_path : str, optional
            Directory where the scenario will be written to, by default None.
        inspect : bool, optional
            If true, this function will return the dataframes, by default False.
        fileformat : str, optional
            Format of the inherited files, by default "csv".
        map_format : str, optional
            The shape of the map file, using "dat" will produce a GAMS-compatible map, by default "dat".
        profile_format : str, optional
            If "wide", the time indeces will be used as columns, if "tall", there will be only one value column,
            by default "wide".
        support_sets : bool, optional
            If true, structural validation sets will be part of the output, by default False.

        Returns
        -------
        pandas.dataframe
            In case inspect is True, else returns NoneType.

        Example
        -------
        m = Instance.from_path("./data")
        m.write(fileformat="dat", output_path="./data2", float_format="{:.4g}".format)
        """
        dataframes = kwargs.get("dataframes", None)
        if dataframes is None:
            dataframes = self.dataframes
        sets = dataframes["sets"]
        parameters = dataframes["parameters"]
        timeseries = dataframes["timeseries"]
        maps = dataframes["maps"]
        extension = f".{fileformat}"
        if not os.path.isdir(output_path):
            os.makedirs(output_path)
        if support_sets:
            # FIXME: We could use the schemas not to hard code this
            sets["set_timedata"] = build_timeseries_range()
            sets["set_profiletypes"] = ["lower", "fixed", "upper"]
            sets["set_accnodesdata"] = list(set(dataframes["sets"].get("set_nodesdata", []) + dataframes["sets"].get("set_accnodes", []) + ["global"]))
            sets["set_linksdata"] = list(set(dataframes["sets"].get("set_linksdata", []) + ["global"]))
            sets["set_years"] = list(set(dataframes["sets"].get("set_years", []) + ["horizon"]))

        for ifile, setlist in sets.items():
            self.write_set_file(
                data=setlist, parent=output_path, name=ifile, extension=extension, **kwargs
            )
        for ifile, dataframe in parameters.items():
            self.write_parameter_file(
                data=dataframe, parent=output_path, name=ifile, extension=extension, **kwargs
            )
        for ifile, dataframe in timeseries.items():
            self.write_profile_file(
                data=dataframe, parent=output_path, name=ifile, extension=extension, profile_format=profile_format, **kwargs
            )

        for ifile, dataframe in maps.items():
            self.write_map_file(
                data=dataframe, parent=output_path, name=ifile, extension=extension, map_format=map_format, **kwargs
            )

        dataframes = {
            "sets": sets,
            "maps": maps,
            "parameters": parameters,
            "timeseries": timeseries,
        }
        return dataframes if inspect else None

    def write_map_file(self, data, parent, name, extension, **kwargs):
        """Write a map file.

        Parameters
        ----------
        data : pandas.Dataframe, list
            Data to write.
        parent : str
            Path to the parent directory.
        name : str
            Name of the map file
        extension : str
            Datatype extension, e.g. .csv, .dat.
        """
        self.write_parameter_file(data, parent, name, extension, **kwargs)

    def write_parameter_file(self, data, parent, name, extension, **kwargs):
        """Write a parameter file.

        Parameters
        ----------
        data : pandas.Dataframe, list
            Data to write.
        parent : str
            Path to the parent directory.
        name : str
            Name of the parameter file.
        extension : str
            Datatype extension, e.g. .csv, .dat.
        """
        filename = name + extension
        self.logger(f"Inheriting data for file {name}")
        if isinstance(data, list):
            data = self._merge_dfs(data)

        if len(data.index) != 0:
            map_format = kwargs.get("map_format", None)
            float_format = kwargs.get("float_format", None)
            if (map_format == "dat") | (extension == ".dat"):
                to_dat(data, join(parent, filename), float_format=float_format)
            else:
                no_sets = kwargs.get("no_sets", False)
                if no_sets:
                    data.index.names = len(data.index.names) * [""]
                    data.to_csv(join(parent, filename), float_format=float_format)
                else:
                    data.to_csv(join(parent, filename), float_format=float_format)

    def write_set_file(self, data, parent, name, extension, **kwargs):
        """Write a set file.

        Parameters
        ----------
        data : pandas.Dataframe, list
            Data to write.
        parent : str
            Path to the parent directory.
        name : str
            Name of the set file
        extension : str
            Datatype extension, e.g. .csv, .dat.
        """
        filename = name + extension
        if len(data) > 0:
            with open(join(parent, filename), "w", encoding="UTF-8") as out_data:
                self.logger(f"Inheriting data for file {name}")
                for element in data:
                    out_data.write(str(element) + "\n")

    def write_profile_file(self, data, parent, name, extension, **kwargs):
        """Write a profile file.

        Parameters
        ----------
        data : pandas.Dataframe, list
            Data to write.
        parent : str
            Path to the parent directory.
        name : str
            Name of the profile file
        extension : str
            Datatype extension, e.g. .csv, .dat.
        """
        filename = name + extension
        if isinstance(data, list):
            data = self._merge_dfs(data)

        roundts = kwargs.get("roundts", 0)
        if roundts:
            data = get_limited_df(data, name, _logger=self.logger)
        if not data.empty:
            self.logger(f"Inheriting data for file {name}")
            float_format = kwargs.get("float_format", None)
            if extension == ".dat":
                to_dat(data, join(parent, filename), float_format=float_format)
            else:
                profile_format = kwargs.get("profile_format", "wide")
                if profile_format == "tall":
                    old_index = list(data.index.names)
                    data = DataFrame({"Value": data.stack()})
                    if any(old_index):
                        old_index.append("timeData")
                        data.index.names = old_index
                no_sets = kwargs.get("no_sets", False)
                if no_sets:
                    data.index.names = len(data.index.names) * [""]
                    data.to_csv(join(parent, filename), float_format=float_format)
                else:
                    data.to_csv(join(parent, filename), float_format=float_format)

    @classmethod
    def from_dataframes(cls, dataframes: dict, version="latest", name="default", **kwargs):
        """Factory method to create an instance of AbstractDataset from a dictionary of DataFrames.

        Parameters
        ----------
        dataframes : dict
             Dictionary of input DataFrames.
        version : str, optional
            REMix model version, by default "latest".
        name : str, optional
            name of the scenario, by default "default".

        Returns
        -------
        classinstance
            Instance of Dataset or any childclass
        """
        instance = cls(version=version, name=name, dataframes=dataframes, **kwargs)
        instance.ready = True
        return instance


class DirectoryDataset(Dataset):
    """REMix scenario directory dataset.

    This class's main function is to contain a REMix scenario dataset stored in a file directory.

    It is mainly designed to be used with its :py:meth:`remix.framework.api.dataset.DirectoryDataset.from_path` method
    which parses a file path containing .dat or .csv files into a dictionary of DataFrames. csv files read by this
    class have the legacy form without index names.

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
        Path where the items input files are located, by default None.
    roundts : int, optional
        If 1, the output functions will round time series values, by default 0.
    version : str, optional
        REMix model version, by default "latest".
    """
    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)

        self._schemas = None
        self.errors = kwargs.get("errors", [])
        self.own_files = []

    def _find_files(self) -> None:
        """Find all the files with the extension .dat. If none are found it will look for the csv extension instead."""
        files = list(Path(self.path).glob("*.dat"))
        if len(files) == 0:
            files = list(Path(self.path).glob("*.csv"))
        files = [str(p) for p in files]
        self.own_files.extend(files)
        self.files.extend(files)
        self.files = list(set(self.files))

    def _arrange_files(self) -> None:
        """If multiple subscenarios are present, it ensures all file names are lowercase and creates a dictionary with
        file names, which are the same but belong to different subscenarios.
        """
        self.arranged_files = {}
        for ifile in self.files:
            file_name = str(Path(ifile).name).lower()
            if file_name not in self.arranged_files:
                self.arranged_files[file_name] = []
            self.arranged_files[file_name].append(ifile)

    def _split_parameters_sets(self) -> None:
        """Separate the different types of input files: sets, parameters, profiles and maps.

        This only checks for file names.
        """
        self.sets = {
            k: v for k, v in self.arranged_files.items() if k.lower().startswith("set_")
        }
        parameters = {
            k: v
            for k, v in self.arranged_files.items()
            if not (k.lower().startswith("set_"))
        }
        parameters = {k: v for k, v in parameters.items() if not k.lower() == "version"}
        self.parameters = {
            k: v for k, v in parameters.items() if not (is_timeseries(k) or is_map(k))
        }
        self.timeseries = {k: v for k, v in parameters.items() if is_timeseries(k)}
        self.maps = {k: v for k, v in parameters.items() if is_map(k)}

    def _search_duplicates(self) -> None:
        """Checks for the possibility of an input file to be in the dataset twice.

        Mostly related to capitalized titles being considered different files in most systems.

        Raises
        ------
        RuntimeError
            If some file was found to exist twice, it will stop the program.
        """
        duplicates = [
            str(ifile).lower()
            for ifile in self.files
            if len([f for f in self.files if str(f).lower() == str(ifile).lower()]) > 1
        ]
        if len(duplicates) > 0:
            raise RuntimeError("Duplicate files found in the input")

    def _inherit_dataframes(self):
        """Run the inheritance on all dataframes"""
        sets = self._inherit_sets()
        parameters = self._inherit_parameters()
        timeseries = self._inherit_timeseries()
        maps = self._inherit_maps()
        self.dataframes = {
            "sets": sets,
            "parameters": parameters,
            "timeseries": timeseries,
            "maps": maps,
        }

    def _inherit_sets(self):
        """Extract the data from the set input files.

        These are mostly csv with one entry and no header.
        """
        files = {}
        for name, path_list in self.sets.items():
            if len(path_list) > 1:
                path_list.sort(key=lambda x: len(Path(x).parents))
                self.logger(f"{name} was inherited from {path_list[-1]}")
            dataframe = self.read_file(
                path_list[-1]
            )
            setlist = self._set_dataframe_to_list(dataframe)
            _name = name.replace(".dat", "").replace(".csv", "").lower()
            # FIXME: Remove global and horizon from "working" dataframe?
            if _name in ["set_accnodesdata", "set_linksdata", "set_years"]:
                setlist = [set_element for set_element in setlist if set_element not in ["global", "horizon"]]
            files[_name] = setlist

        return files

    def _set_dataframe_to_list(self, dataframe):
        """
        Get the index of a dataframe as list.

        Parameters
        ----------
        dataframe : pandas.DataFrame
            Pandas DataFrame to extract the index of.

        Returns
        -------
        list
            DataFrame index converted to a list.
        """
        if dataframe is None:
            return []
        else:
            return dataframe.index.tolist()

    def _inherit_maps(self):
        files = {}
        for name, path_list in self.maps.items():
            if len(path_list) > 1:
                path_list.sort(key=lambda x: len(Path(x).parents))
                self.logger(f"{name} was inherited from {path_list[-1]}")
            dataframe = self.read_file(
                path_list[-1]
            )
            files[name.replace(".dat", "").replace(".csv", "").lower()] = dataframe
        return files

    def _inherit_parameters(self):
        """Extract the data from the input files.

        The input files can be maps, which usually are point delimited files with two columns or they can be
        parameters, which are tabular data in .dat or .csv format.

        Parameters
        ----------
        data : dict
            Dictonary of names and paths to the data.
        """
        files = {}
        for name, path_list in self.parameters.items():
            path_list.sort(key=lambda x: len(Path(x).parents))
            df_list = [
                self.read_file(path) for path in path_list
            ]
            df_list = [frame for frame in df_list if frame is not None]
            file_schema = self._get_file_schema(path_list[0])
            dataframe = self._merge_dfs(df_list, file_schema)
            files[name.replace(".dat", "").replace(".csv", "").lower()] = dataframe
        return files

    def _inherit_timeseries(self):
        """Extract the data from the profile input files.

        Profiles have a special format very similar to parameters but with a time index.
        """
        files = {}
        for name, path_list in self.timeseries.items():
            path_list.sort(key=lambda x: len(Path(x).parents))
            df_list = [
                self.read_file(path) for path in path_list
            ]
            df_list = [frame for frame in df_list if frame is not None]
            file_schema = self._get_file_schema(path_list[0])
            dataframe = self._merge_dfs(df_list, file_schema)
            files[name.replace(".dat", "").replace(".csv", "").lower()] = dataframe
        return files

    def _get_file_schema(self, path=None):
        """Placeholder method to get file schema.

        This method is not requried in this class, but allows children classes to override it in order to implement
        different behavior.

        Parameters
        ----------
        path : Nonetype, optional
            Path to the file schema, by default None.

        Returns
        -------
        Nonetype
            This placeholder always returns None.
        """
        return None

    def _merge_dfs(self, df_list, file_schema=None):
        """Placeholder method for merging dataframes.

        This method pipes the default :code:`merge_dfs` functionality from
        :py:func:`remix.framework.tools.utilities.merge_dfs`. It does not need the name of the file schema but the
        method can be overriden by children classes.

        Parameters
        ----------
        df_list : list
            List of dataframes to be merged.
        file_schema : str, optional
           Name of the file schema, by default None.

        Returns
        -------
        pandas.DataFrame
            Merged DataFrame.
        """
        dataframe = merge_dfs(df_list)
        if set(dataframe.columns) == {"Value"}:
            dataframe = dataframe.unstack()
            dataframe.columns = dataframe.columns.get_level_values(1)
        return dataframe

    def _get_limited_df(self, df, file, _logger=None):
        """Placeholder method to get a limited DataFrame.

        This method pipes the default :code:`get_limited_df` functionality from
        :py:func:`remix.framework.tools.utilities.get_limited_df`. It does not need the name of the file schema but the
        method can be overriden by children classes.

        Parameters
        ----------
        df : pandas.DataFrame
            Original DataFrame from raw input data.
        file : str
            Path of the input file.
        _logger : FunctionType, optional
            A logging function that takes strings as input, by default lambda x: None.
        file_schema : dict, optional
            Schema for the current file, only a placeholder in this method, by default None.

        Returns
        -------
        pandas.DataFrame
            Limited DataFrame.
        """
        return get_limited_df(df, file, _logger=_logger)

    def inherit(self):
        """Build dataframes based on the given input files."""
        if isinstance(self.path, str):
            self.path = Path(self.path)
        self._arrange_files()
        self._split_parameters_sets()
        self._search_duplicates()
        self._inherit_dataframes()
        self.ready = True

    def read_file(self, input: str) -> DataFrame:
        """This interface was added to have flexible input file reading."""
        try:
            dataframe = read_file(input)
        except ValueError as e:
            # FIXME: Very problematic, use something that is not a try except clause.
            if not ".csv" in input:
                raise  e
            if self._schemas is None:
                from remix.framework.schema.templates import sm
                self._schemas = {k.lower(): v for k, v in sm.get_schemas().items()}
            file_schema = self._schemas[Path(input).stem.lower()]
            dataframe = read_remix_csv(input, schema=file_schema)
        return dataframe

    @classmethod
    def from_path(
        cls,
        path: Union[str, Path],
        parent_files: list = None,
        version: str = "latest",
        name: str = ".",
        roundts: int = 1,
        logger: FunctionType = lambda x: None,
        inherit: bool = False,
        **kwargs
    ):
        """Factory method to create an instance of DirectoryDataset from a file path.

        Parameters
        ----------
        path : str, Path
            String or path of the input files to be parsed.
        parent_files : list, optional
            Files of a parent scenario.
        version : str, optional
            REMix model version, by default "latest".
        name : str, optional
            name of the scenario, by default ".".
        roundts : int, optional
            Round the timeseries to match the GAMS maximum file-length restriction, by default 1.
        logger : FunctionType, optional
            A logging function that takes strings as input, by default lambda x: None.
        inherit : bool, optional
            Inherit the scenario, by default False.

        Returns
        -------
        classinstance
            Instance of the desired class.
        """
        if parent_files is None:
            parent_files = []

        instance = cls(
            path=path,
            parent_files=parent_files,
            version=version,
            name=name,
            logger=logger,
            roundts=roundts,
            **kwargs
        )
        instance.files.extend(parent_files)
        instance._find_files()
        if inherit:
            instance.inherit()
        return instance

class DataPackage:
    """REMix data package. Container of one or more scenarios.

    Parameters
    ----------
    path : str, optional
        Main directory of the package, by default None.
    version : str, optional
        REMix version, by default "latest".
    dataset_calss : Class, optional
        Dataset managing class, by default DirectoryDataset.
    """

    def __init__(self, path=None, version="latest", **kwargs):
        self.path = path
        self.version = version
        self.logger = kwargs.get("logger", lambda x: None)
        self.scenarios = kwargs.get("scenarios", {})
        if isinstance(self.path, str):
            self.path = Path(self.path)
        self.dataset_class = kwargs.get("dataset_class", DirectoryDataset)
        self.find_scenarios()

    def get_scenario(self, scn_name, inherit=True):
        """Get a scenario of the DataPackage instance.

        Parameters
        ----------
        scn_name : str
            Name of the scenario.
        inherit : bool, optional
            inherit data if True, by default True

        Returns
        -------
        DirectoryDataset
            The scenario with its data.

        Raises
        ------
        KeyError
            If scn_name not in the available scenarios.
        """
        if scn_name not in self.scenarios:
            msg = (
                "Could not find specified scenario. Available scenario names are: \"" +
                "\", \"".join(self.scenarios.keys()) + "\"."
            )
            raise KeyError(msg)
        if inherit:
            self.scenarios[scn_name].inherit()
        return self.scenarios[scn_name]

    def find_scenarios(self):
        """Find all the scenarios and its subscenarios in the data folder."""
        current = self.dataset_class.from_path(
            self.path,
            version=self.version,
            logger=self.logger,
            inherit=False,
        )
        self.scenarios["."] = current
        subscenarios = [
            Path(p)
            for p in glob.glob(f"{self.path.as_posix()}/**/", recursive=True)
            if Path(p) != self.path
        ]
        parent = current
        for sub in subscenarios:
            name = sub.relative_to(self.path).as_posix()
            parent = self.scenarios[sub.parent.relative_to(self.path).as_posix()]
            self.scenarios[name] = self.dataset_class.from_path(
                sub,
                parent_files=parent.files,
                version=self.version,
                name=name,
                logger=self.logger,
                inherit=False,
            )

    def inherit_scenario(self, scenario, output_path=None, fileformat="dat", roundts=0):
        """Inherit one scenario of the input files.

        Parameters
        ----------
        scenario : str
            Main target directory.
        output_path : str, optional
            Path to the base scenario, by default None.
        fileformat : str, optional
            File format for the data.

        Returns
        -------
        dict
            Dictionary containing the scenario's data in pandas.DataFrames.
        """
        if not self.scenarios[scenario].ready:
            self.scenarios[scenario].inherit()
        if output_path is not None:
            dataframes = self.scenarios[scenario].write(
                output_path=output_path, inspect=True, fileformat=fileformat, roundts=roundts
            )
        else:
            dataframes = self.scenarios[scenario].dataframes
        return dataframes

    def inherit_all(self, output_path, fileformat="dat", roundts=0):
        """Inherit all scenarios in the package to individual directories.

        Parameters
        ----------
        output_path : str
            Path to the base scenario.
        fileformat : str, optional
            File format for the data, by default "dat".
        """
        for scenario in self.scenarios:
            output_name = scenario.replace("/", "_")
            if output_name == ".":
                if "base" not in self.scenarios:
                    output_name = "base"
                else:
                    output_name = "scenBase"
            scenario_output = join(output_path, f"{output_name}")
            if not exists(scenario_output):
                Path(scenario_output).mkdir(parents=True, exist_ok=True)
            self.inherit_scenario(
                scenario=scenario, output_path=scenario_output, fileformat=fileformat, roundts=roundts
            )

    def export(self, output_path, fileformat="csv", map_format="dat", profile_format="wide", support_sets=False, **kwargs):
        """This function writes the stored scenarios as a folder structure in the output_directory.

        Parameters
        ----------
        output_path : str
            Path to the base scenario.
        fileformat : str, optional
            Format for the files, by default "csv".
        map_format : str, optional
            Formal for the map files, by default "dat".
        profile_format : str, optional
            "tall" or "wide" time series format, by default "wide".
        support_sets : bool, optional
            If true, timedata and profile types will be part of the model, by default False.
        """
        # export base scenario first
        if not self.scenarios["."].ready:
            self.scenarios["."].inherit()
        self.scenarios["."].write(
            output_path=output_path,
            inspect=False,
            fileformat=fileformat,
            map_format=map_format,
            profile_format=profile_format,
            support_sets=support_sets,
            **kwargs
        )

        for k in self.scenarios:
            if k != ".":
                scenario_folder = Path(output_path).joinpath(k)
                scenario_folder.mkdir(exist_ok=True, parents=True)
                tmp_subscenario = self.dataset_class.from_path(
                    Path(self.path).joinpath(k),
                    parent_files=[],
                    version=self.version,
                    name=k,
                    logger=self.logger,
                    inherit=True,
                    **kwargs
                )
                tmp_subscenario.write(
                    output_path=scenario_folder,
                    inspect=False,
                    fileformat=fileformat,
                    map_format=map_format,
                    profile_format=profile_format,
                    infer=False,
                    support_sets=False,
                    **kwargs
                )

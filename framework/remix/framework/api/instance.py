"""

"""
import os

import numpy as np
from pandas import DataFrame
from pandas import MultiIndex
from pandas.errors import IntCastingNaNError

from remix.framework.api.dataset import DirectoryDataset
from remix.framework.api.run import default_run_remix_kwargs
from remix.framework.api.run import run_remix
from remix.framework.schema.templates import default_dataframes
from remix.framework.schema.templates import distribute_dictionary
from remix.framework.schema.templates import extract_sets
from remix.framework.schema.templates import sm
from remix.framework.tools.gdx import GDXEval
from remix.framework.tools.utilities import is_empty
from remix.framework.tools.utilities import merge_dfs
from remix.framework.tools.utilities import merge_dicts

NON_FILE_SETS = [
    "set_profiletypes",
]  # Sets that are in the parameters but are no input files


class Instance(DirectoryDataset):
    """
    This class sets up the data structure necessary to build a full energy system model in REMix.

    The data is set up in pandas.DataFrames and can be written to different file formats (dat or csv). Sets are
    auto-generated based on the input in the parameters, profiles and maps. It is possible to provide default values to
    missing values in the data.
    """

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        as_lists = kwargs.get("as_lists", False)
        self.custom_defaults = kwargs.get("custom_defaults", None)
        self.schemas = {k.lower(): v for k, v in sm.get_schemas().items()}
        self.assign_attributes(
            descriptor=default_dataframes(
                as_lists=as_lists,
                index_names=kwargs.get("index_names", True),
                column_names=kwargs.get("column_names", True),
                schemas=self.schemas,
            )
        )

        datafields = {
            "sets": "set",
            "parameters": "parameter",
            "maps": "map",
            "timeseries": "profile",
        }

        if self.dataframes is not None:
            for k, v in self.dataframes.items():
                for df_key, df in v.items():
                    obj = getattr(self, datafields[k])
                    obj.add(df, df_key.replace("set_", "").replace("map_", ""))

        self.update_dataframes()
        self.datadir = kwargs.get("datadir", "./data")
        self.config = kwargs.get("config", None)
        self.ready = True
        self.result = {}

    def _inherit_dataframes(self):
        """Build and update the data of an instance."""
        sets = self._inherit_sets()
        parameters = self._inherit_parameters()
        timeseries = self._inherit_timeseries()
        maps = self._inherit_maps()
        descriptor = {
            "set": sets,
            "parameter": parameters,
            "profile": timeseries,
            "map": maps,
        }
        for k, v in descriptor.items():
            for sk in v.keys():
                obj = getattr(self, k)
                obj.add(descriptor[k][sk], sk.replace(f"{k}_", ""))

        self.update_dataframes()

    def run(self, resultfile="remix", read_result=False, **kwargs):
        """Method to run a REMix.Instance.

        This method pipes all necessary data of the Instance to the :py:func:`remix.framework.api.run.run_remix`
        function. If :code:`read_results` is True, the result file is read after the optimization and results are
        available in the :code:`.result` attribute of the Instance object.

        A list with the most important parameters for the :code:`run` method, can be obtained from by typing
        `remix run --help` in command line interface. A full list of all important remix parameters can be found
        :ref:`here <cli_label>`.

        Parameters
        ----------
        resultfile : str, optional
            Name of the result file, by default "remix".
        read_result : bool, optional
            Read the results after the optimization, by default False.

        Returns
        -------
        int
            Return value of the subprocess running the GAMS code.
        """
        status = run_remix(
            datadir=self.datadir,
            resultfile=resultfile,
            **kwargs,
        )

        if read_result:
            resultdir = kwargs.get("resultdir", default_run_remix_kwargs["resultdir"])
            self.result = GDXEval(f"{resultdir}/{resultfile}.gdx")
        return status

    def write(self, output_path=None, inspect=False, fileformat="csv", infer=True, reassign=False, **kwargs):
        """Main inheritance method. Converts dataframes into .csv or .dat files.

        Other parameters available are documented in the method of the parent class:
        :py:meth:`remix.framework.api.dataset.Dataset.write`

        Parameters
        ----------
        output_path : str, optional
            Directory the scenario will be written to, by default None.
        inspect : bool, optional
            If True, this function will return the DataFrames, by default False.
        fileformat : str, optional
            Format of the inherited files, by default "csv".
        infer : bool, optional
            Autogenerate sets, by default True.
        reassign : bool, optional
            When writing a scenario assign the output_path as datadir property, by default False.

        Returns
        -------
        pandas.DataFrame
            In case inspect is True, else returns NoneType.
        """
        if infer:
            self.infer_set_data()
        self.update_dataframes()
        if output_path is None:
            output_path = self.datadir
        else:
            if reassign:
                self.datadir = output_path
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        dataframes = self.build_dataframes_to_write()
        dataframes = super().write(output_path, inspect=inspect, fileformat=fileformat, dataframes=dataframes, **kwargs)

        return dataframes

    def infer_set_data(self):
        """Value-inference functionality for sets."""
        serialized_params = self.parameter.serialize()
        serialized_profiles = self.profile.serialize()
        serialized_maps = self.map.serialize()

        non_empty_parameters = {k.lower(): v for k, v in serialized_params.items() if not is_empty(v)}
        non_empty_profiles = {k.lower(): v for k, v in serialized_profiles.items() if not is_empty(v)}
        non_empty_maps = {k.lower(): v for k, v in serialized_maps.items() if not is_empty(v)}

        infer_sets = {}
        for element in [non_empty_parameters, non_empty_profiles]:
            for parameter in element:
                # ASSUMING IS ORDERED, THIS IS UGLY, MAYBE CONSIDER HAVING NAMES IN SETS
                if parameter in self.schemas:
                    parameter_sets = extract_sets(self.schemas[parameter])
                    parameter_sets = [p for p in parameter_sets if p not in ["set_timedata"]]
                    if isinstance(element[parameter], list):
                        element[parameter] = self._merge_dfs(element[parameter])
                    else:
                        element[parameter] = element[parameter]
                    if isinstance(element[parameter].index, MultiIndex):
                        for i in range(len(element[parameter].index.levels)):
                            level = element[parameter].index.get_level_values(i).unique()
                            current_level = [l for l in level if l not in ["global", "horizon"]]
                            if parameter_sets[i] not in infer_sets:
                                infer_sets[parameter_sets[i]] = set(current_level)
                            else:
                                infer_sets[parameter_sets[i]].update(current_level)
                    else:
                        current_level = [l for l in element[parameter].index.tolist() if l not in ["global", "horizon"]]
                        if parameter_sets[0] not in infer_sets:
                            infer_sets[parameter_sets[0]] = set(current_level)
                        else:
                            infer_sets[parameter_sets[0]].update(current_level)
        for _map in non_empty_maps:
            map_name = f"map_{_map}"
            if map_name in self.schemas:
                map_sets = extract_sets(self.schemas[map_name])
                if isinstance(non_empty_maps[_map], list):
                    for list_elem in non_empty_maps[_map]:
                        infer_sets = self._infer_sets_from_maps_dfs(infer_sets, list_elem, map_sets)
                else:
                    infer_sets = self._infer_sets_from_maps_dfs(infer_sets, non_empty_maps[_map], map_sets)

        infer_sets = {
            k: list(v) for k, v in infer_sets.items() if k.lower() not in ["set_profiletypes", "set_accnodesdata"]
        }

        # write automatically generated sets into .set attribute
        for _set, data in infer_sets.items():
            if _set.lower() not in NON_FILE_SETS:
                self.set.add(data, _set[4:])

    def _infer_sets_from_maps_dfs(self, infer_sets, df, map_sets):
        """Generate sets from map data.

        Parameters
        ----------
        infer_sets : dict
            Dictionary with set names and set data to be extended from the map information.
        df : pandas.DataFrame
            pandas DataFrame to infer from.
        map_sets : list
            List of sets of a specific map.

        Returns
        -------
        dict
            Updated infer_sets.
        """
        if isinstance(df.index, MultiIndex):
            for i, level in enumerate(df.index.levels):
                if map_sets[i] not in infer_sets:
                    infer_sets[map_sets[i]] = set(level)
                else:
                    infer_sets[map_sets[i]].update(level)
        return infer_sets

    def update_dataframes(self):
        """Update DataFrames using the current state of containers."""
        sets = self.set.serialize()
        sets = {f"set_{k}": v for k, v in sets.items()}
        parameters = self.parameter.serialize()
        timeseries = self.profile.serialize()
        maps = self.map.serialize()
        maps = {f"map_{k}": v for k, v in maps.items()}
        dataframes = {
            "sets": sets,
            "maps": maps,
            "parameters": parameters,
            "timeseries": timeseries,
        }
        self.dataframes = dataframes

    def build_dataframes_to_write(self):
        """Build the DataFrames before writing.

        This method allows to fill missing values with default values and checking for datatypes.

        Returns
        -------
        dict
            Dictionary of the sets, maps, parameters and time series DataFrames.
        """
        sets = self.set.serialize()
        sets = {f"set_{k}": v for k, v in sets.items()}
        parameters = self.parameter.get_dataframes_with_defaults()
        timeseries = self.profile.serialize()
        maps = self.map.serialize()
        maps = {f"map_{k}": v for k, v in maps.items()}
        dataframes = {
            "sets": sets,
            "maps": maps,
            "parameters": parameters,
            "timeseries": timeseries,
        }
        return dataframes

    def assign_attributes(self, descriptor=None):
        """This method is the main updating method.

        If it is run without an argument, it will assign its DataFrames to its Container objects. If it is run with a
        descriptor, the descriptor can be used to initialize the Containers. It is recommended that this last option is
        not done outside of the from_dataframes method.

        Parameters
        ----------
        descriptor : dict
            Descriptor containing REMix dataframes in deep form, by default None.
        """
        if descriptor is None:
            short_sets = {k.replace("set_", ""): v for k, v in self.dataframes["sets"].items()}
            self.set = container_factory(short_sets, "Sets")()
            self.parameter = container_factory(self.dataframes["parameters"], "Parameters")(
                self.schemas, custom_defaults=self.custom_defaults
            )
            self.profile = container_factory(self.dataframes["timeseries"], "Profiles")()
            short_maps = {k.replace("map_", ""): v for k, v in self.dataframes["maps"].items()}
        else:
            short_sets = {k.replace("set_", ""): v for k, v in descriptor["sets"].items()}
            self.set = container_factory(short_sets, "Sets")()
            self.parameter = container_factory(descriptor["parameters"], "Parameters")(
                self.schemas, custom_defaults=self.custom_defaults
            )
            self.profile = container_factory(descriptor["timeseries"], "Profiles")()
            short_maps = {k.replace("map_", ""): v for k, v in descriptor["maps"].items()}

        self.map = container_factory(short_maps, "Maps")()

    @classmethod
    def from_path(
        cls,
        path,
        parent_files=None,
        version="latest",
        name=".",
        roundts=1,
        logger=lambda x: None,
        inherit=True,
        custom_defaults=None,
        index_names=True,
        **kwargs,
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
            Round the timeseries to match the GAMS maximum file length restriction, by default 1.
        logger : FunctionType, optional
            A logging function that takes strings as input, by default lambda x: None.
        inherit : bool, optional
            Inherit the scenario, by default False.
        as_lists : bool, optional
            Read the DataFrames nested in lists if True, by default False.
        custom_defaults : dict, optional
            Dictionary with custom default values for some data fields, by default None.
        index_names : bool, optional
            Read with index names if True, by default False.

        Returns
        -------
        classinstance
            Instance of the desired class.

        Example
        -------
        m = Instance.from_path(path="./data")
        """
        kwargs["index_names"] = index_names
        instance = super().from_path(
            path=path,
            parent_files=parent_files,
            version=version,
            name=name,
            roundts=roundts,
            logger=logger,
            inherit=inherit,
            custom_defaults=custom_defaults,
            **kwargs,
        )
        if not inherit:
            instance.ready = False
        instance.datadir = path
        return instance

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
        defaults_deep = default_dataframes(version=version, flat=False)
        if dataframes.keys() == defaults_deep.keys():
            for element in dataframes:
                missing = [ky for ky in dataframes[element] if ky not in defaults_deep[element]]
                if len(missing) > 0:
                    mssn = "[" + " ".join(missing).strip(" ") + "]"
                    raise KeyError(f"the keys: {mssn} at {element} are not valid")
        else:
            defaults_shallow = default_dataframes(version=version, flat=True)
            missing = [ky for ky in dataframes if ky not in defaults_shallow]
            if len(missing) > 0:
                mssn = "[" + " ".join(missing).strip(" ") + "]"
                raise KeyError(f"the keys: {mssn} are not valid")
            dataframes = distribute_dictionary(dataframes)

        return super().from_dataframes(dataframes, version, name, **kwargs)


def container_factory(descriptor, name="Container"):
    """Dynamic instantiation of classes is very tricky. This factory function takes care of this.

    Parameters
    ----------
    descriptor : Mapping
        A mapping describing the properties of the instantiated class.
    name : str, optional
        Name of the container class, by default "Container".
    custom_defaults : dict, optional
        Dictionary with default values for some datafields, by default None.

    Returns
    -------
    object
        container_factory Class or Sets instance.
    """

    class GenericDC:
        """Generic DataContainer class."""

        def __init__(self):
            self.name = name

        def __repr__(self):
            return f"{self.name} with the following items: {self._attrs}"

        def __str__(self):
            return f"{self.name} with the following items: {self._attrs}"

        def _check_label_on_add(self, label):
            """Check if the label of the data field specified is valid.

            Parameters
            ----------
            label : str
                Label for the data field.

            Raises
            ------
            AttributeError
                If the label is available.
            """
            if label not in self._attrs:
                msg = (
                    f"The Instance {self.name} does not have any data with the label '{label}'. Please check, if you "
                    "specified the correct label."
                )
                raise AttributeError(msg)

        def _make_attrs(self):
            """Extract all datafields of a class object."""
            _attrs = [a for a in dir(self) if not a.startswith("_") and not callable(getattr(self, a))]
            self._attrs = [a for a in _attrs if a not in ["name", "custom_defaults"]]
            self._org_cols_params = {a: getattr(self, a) for a in self._attrs}

        def add(self, new_data, label):
            """Add new data to the existing data while checking the validity of the new data beforehand.

            Parameters
            ----------
            new_data : pd.DataFrame, list
                New data to be added to the data container.
            label : str
                Name of the data to be added.
            """
            self._check_label_on_add(label)
            self._check_datatype_on_add(new_data)

            # Enforce dataframes have only Strings as index
            if isinstance(new_data, DataFrame):
                new_data = new_data.copy()
                if len(new_data.index) > 0:
                    if isinstance(new_data.index, MultiIndex):
                        for i, j in enumerate(new_data.index.levels):
                            new_data.index = new_data.index.set_levels(j.astype(str), level=i)
                    else:
                        new_data.index = new_data.index.astype(str)

            elif isinstance(new_data, list):
                new_data = [str(i) for i in new_data]

            current_data = getattr(self, f"_{label}")
            self._add(current_data, new_data, label)

        def serialize(self):
            """Convert the container into a dictionary.

            Returns
            -------
            dict
                Dictionary containing data fields only.
            """
            return {s: self.__dict__[f"_{s}"] for s in self._attrs}

    class ListDC(GenericDC):
        """Generic DataContainer class saving data in lists."""

        def __init__(self):
            super().__init__()

            for key, value in descriptor.items():
                if isinstance(value, list):
                    setattr(self, "_" + key.lower(), value)
                else:
                    wrong_type = type(value)
                    raise TypeError(f"The element {k} has type {wrong_type} which is not valid, use a list instead.")

            self._make_attrs()

        def _add(self, current_data, new_data, label):
            """Concatenate the new data given with the existing data.

            Parameters
            ----------
            current_data : list
                Existing data of the model.
            new_data : list
                New data to be added to the model.
            label : str
                Name of the data field.
            """
            current_data.extend(new_data)
            setattr(self, f"_{label}", sorted(list(set(current_data))))

        def _check_datatype_on_add(self, data):
            """Check if data types of new data are correct.

            Parameters
            ----------
            data : list
                Data to be added to the model.

            Raises
            ------
            TypeError
                If the data is not a list.
            """
            if not isinstance(data, list):
                wrong_type = type(data)
                raise TypeError(f"Datatype {wrong_type} cannot be added, please use a list.")

    class DataFrameDC(GenericDC):
        """Generic container class storing data in pandas.DataFrames"""

        def __init__(self):
            super().__init__()

            for key, value in descriptor.items():
                if isinstance(value, (DataFrame, list)):
                    setattr(self, "_" + key.lower(), value)
                else:
                    wrong_type = type(value)
                    raise TypeError(
                        f"The element {k} has type {wrong_type} which is not valid, try list or pandas.DataFrame"
                    )

            self._make_attrs()

        def _add(self, current_data, new_data, label):
            """Concatenate the new data given with the existing data.

            Parameters
            ----------
            current_data : pd.DataFrame
                Existing data of the model.
            new_data : pd.DataFrame
                New data to be added to the model.
            label : str
                Name of the data field.

            Raises
            ------
            ValueError
                In case illegal column names are provided in the new data.
            """
            if isinstance(current_data, DataFrame) and isinstance(new_data, DataFrame):
                legal = self._org_cols_params[label].columns
                if not set(new_data.columns).issubset(legal) and len(legal) > 0:
                    illegal = list(set(new_data.columns) - set(legal))
                    msg = (
                        'The following column names were provided, but are not available: "'
                        + '", "'.join(illegal)
                        + f'". The available column names for {label} are: "'
                        + '", "'.join(legal)
                        + '".'
                    )
                    raise ValueError(msg)
                setattr(self, f"_{label}", merge_dfs([current_data, new_data]))
            else:
                # FIXME: this does nothing, does it?
                current_data.append(new_data)

        def _check_datatype_on_add(self, data):
            """Check if data types of new data are correct.

            Parameters
            ----------
            data : pandas.DataFrame
                Data to be added to the model.

            Raises
            ------
            TypeError
                If the data is not a pandas.DataFrame instance.
            """
            if not isinstance(data, (DataFrame, list)):
                wrong_type = type(data)
                raise TypeError(f"Datatype {wrong_type} cannot be added, please use a pandas.DataFrame.")

    class DataFrameDCWithDefaults(DataFrameDC):
        """Generic DataContainer class storing data in pandas.DataFrames with default values."""

        def __init__(self, schemas, custom_defaults=None):
            super().__init__()

            self.schemas = schemas
            self.custom_defaults = custom_defaults
            self._set_default_values_and_dtypes_from_schemas()

        def _set_default_values_and_dtypes_from_schemas(self):
            """Read schemas as set the defaults for all columns of all parameters."""
            required_schemas = {schema: self.schemas[schema] for schema in self.schemas if schema in self._attrs}
            self.defaults = {}
            self.datatypes = {}

            for schema_name in required_schemas:
                index_names = [fk["fields"][0] for fk in required_schemas[schema_name]["foreignKeys"]]
                self.defaults[schema_name] = {
                    f.get("name"): float(f.get("constraints", {}).get("default", np.nan))
                    for f in required_schemas[schema_name]["fields"]
                    if f.get("name") not in index_names
                }
                dtype_lookup = {
                    "number": float,
                    "integer": int,
                    "boolean": int,
                }
                self.datatypes[schema_name] = {
                    f.get("name"): dtype_lookup[f.get("type", "number")]
                    for f in required_schemas[schema_name]["fields"]
                    if f.get("name") not in index_names
                }
            if self.custom_defaults is not None:
                self.defaults = merge_dicts(self.defaults, self.custom_defaults, update_only=True)

        def get_dataframes_with_defaults(self):
            """Get dictionary of DataFrames with missing values replaced by the default values.

            Returns
            -------
            dict
                Dictionary with the DataFrames and missing values filled with defaults.
            """
            dataframes = {}
            for label in self._attrs:
                if isinstance(self.__dict__[f"_{label}"], list):
                    df_out = self.__dict__[f"_{label}"][0]
                    for df in self.__dict__[f"_{label}"][1:]:
                        df_out = self._merge_dfs([df_out, df])

                    dataframes[label] = df_out.fillna(self.defaults[label])
                else:
                    dataframes[label] = self.__dict__[f"_{label}"].infer_objects(copy=False).fillna(self.defaults[label])

                if not dataframes[label].empty:
                    dataframes[label] = dataframes[label].loc[:, ~(np.isnan(dataframes[label]).all())]
                    dataframes[label] = dataframes[label].fillna(0)

                    try:
                        dataframes[label] = dataframes[label].astype(
                            {k: v for k, v in self.datatypes[label].items() if k in dataframes[label].columns}
                        )
                    except IntCastingNaNError:
                        msg = (
                            f"Could not convert the data in a column of {label} to integer. Make sure, that your data "
                            "do not contain any non-convertable datatypes such as inf or NaN."
                        )
                        raise IntCastingNaNError(msg)

            return dataframes

    if name == "Sets":
        c = ListDC
    elif name == "Parameters":
        c = DataFrameDCWithDefaults
    else:
        c = DataFrameDC

    for k, _ in descriptor.items():
        setattr(c, k.lower(), __managed_attribute(k))
    c.__name__ = name.title()
    return c


def __managed_attribute(name: str):
    storage_name = f"_{name.lower()}"

    @property
    def prop(self):
        return getattr(self, storage_name)

    @prop.setter
    def prop(self, value):
        setattr(self, storage_name, value)

    return prop

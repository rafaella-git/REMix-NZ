from pandas import DataFrame

from remix.framework.schema.manager import Manager
from remix.framework.tools.utilities import build_timeseries_range
from remix.framework.tools.utilities import is_map
from remix.framework.tools.utilities import is_timeseries
from remix.framework.tools.utilities import read_remix_csv

IGNORED_SETS = ["set_profiletypes", "set_soc"]
INCLUDED_SETS = ["set_yearssel"]

sm = Manager()


def extract_sets(schema: dict) -> list:
    """Extract the set(Foreign Keys) elements from the schema.

    Parameters
    ----------
    schema : dict
        Schema of the profile.

    Returns
    -------
    list
        List with foreign key elements from the schema.
    """
    this_sets = [f["reference"]["resource"] for f in schema["foreignKeys"]]
    return this_sets


def default_dataframes(
    version: str = "latest",
    index_names: bool = False,
    column_names: bool = False,
    as_lists: bool = False,
    flat: bool = False,
    schemas: dict = None,
) -> dict:
    """Build a dictionary with default DataFrame values.

    This can be used as a reference to build model instances.

    Parameters
    ----------
    version : str, optional
        Model version to be fetched from the schemas, by default "latest".
    index_names : bool, optional
        If True, the index names will appear in the DataFrames, by default False.
    column_names : bool, optional
        If True, the empty DataFrames will be given default columns, by default False.
    flat : bool, optional
        If True, the dictionary won't have its upper level, by default False.
    schemas : dict, optional
        If given, the schema will be used as a reference to build the dataframe, by default None.

    Returns
    -------
    dict
        Dictionary with the requested dataframes.
    """
    if schemas is None:
        schemas = sm.get_schemas(version)
        schemas = {k.lower(): v for k, v in schemas.items()}

    _parameters = {k: v for k, v in schemas.items() if not k.lower().startswith("set_")}
    _parameters = {k: v for k, v in _parameters.items() if k.lower() != "version"}

    parameters = {}
    timeseries = {}
    maps = {}
    for k, v in _parameters.items():
        if not is_timeseries(k) and not is_map(k):
            this_df = empty_dataframe(v, index_names, column_names)
            parameters[k] = [this_df] if as_lists else this_df
        elif is_timeseries(k):
            this_df = empty_profile(v, index_names, column_names)
            timeseries[k] = [this_df] if as_lists else this_df
        elif is_map(k):
            this_df = empty_dataframe(v, index_names, column_names)
            maps[k] = [this_df] if as_lists else this_df
    sets = {k: [] for k in schemas if k.lower().startswith("set_")}
    if flat:
        dictionary = {}
        for element in [parameters, sets, maps, timeseries]:
            dictionary.update(element)
    else:
        dictionary = {
            "sets": sets,
            "parameters": parameters,
            "timeseries": timeseries,
            "maps": maps,
        }
    return dictionary


def distribute_dictionary(dictionary: dict) -> dict:
    """Helper function to transform from a flat format dictionary to a deep one.

    Parameters
    ----------
    dictionary : dict
        Dictionary with flat REMix parameters.

    Returns
    -------
    dict
        Deep dictionary.
    """
    _parameters = {
        k: v for k, v in dictionary.items() if not k.lower().startswith("set_")
    }
    _parameters = {k: v for k, v in _parameters.items() if not k.lower() == "version"}
    parameters = {
        k: v for k, v in _parameters.items() if not (is_timeseries(k) or is_map(k))
    }
    timeseries = {k: v for k, v in _parameters.items() if is_timeseries(k)}
    maps = {k: v for k, v in _parameters.items() if is_map(k)}
    sets = {k: v for k, v in dictionary.items() if k.lower().startswith("set_")}
    new = {
        "sets": sets,
        "parameters": parameters,
        "timeseries": timeseries,
        "maps": maps,
    }
    return new


def empty_dataframe(
    schema: dict, index_names: bool = False, column_names: bool = False
):
    """Helper function to generate an empty DataFrame based on schema.

    Parameters
    ----------
    schema : dict
        Schema of the tabular data.
    index_names : bool, optional
        If True, the index names will be added, by default False.
    column_names : bool, optional
        If True, the column names will be added, by default False.

    Returns
    -------
    pandas.DataFrame
        Empty DataFrame.
    """

    if column_names:
        all_columns = [col["name"] for col in schema["fields"]]
        index_columns = [fk["fields"][0] for fk in schema["foreignKeys"]]
        out_columns = [col for col in all_columns if col not in index_columns]

        complete_columns = index_columns.copy()
        complete_columns.extend(out_columns)
        dataframe = DataFrame({col: [] for col in complete_columns})
        dataframe = dataframe.set_index(index_columns)
        if index_names:
            dataframe.index.names = dataframe.index.names
        else:
            dataframe.index.names = [None] * len(dataframe.index.levels)
    else:
        dataframe = DataFrame()
    return dataframe


def empty_profile(schema: dict, index_names: bool = False, column_names: bool = False):
    """Helper function to build an empty profile dataframe.

    Parameters
    ----------
    schema : dict
        Schema of the profile.
    index_names : bool, optional
        If True, index will have names, by default False.
    column_names : bool, optional
        If True, default time columns will be added, by default False.

    Returns
    -------
    pandas.DataFrame
        Empty DataFrame.
    """
    if column_names:
        dataframe = DataFrame({col["name"]: [] for col in schema["fields"]})
        dataframe = dataframe.set_index(
            [
                fk["fields"][0]
                for fk in schema["foreignKeys"]
                if fk["fields"][0].lower() != "timedata"
            ]
        )
        if index_names:
            dataframe.index.names = dataframe.index.names
        else:
            dataframe.index.names = [None] * len(dataframe.index.levels)
        names = build_timeseries_range()
        profile = DataFrame(columns=names, index=dataframe.index)
    else:
        profile = DataFrame()
    return profile


def infer_sets(template: dict, schemas: dict) -> set:
    """Infer sets from a dictionary of maps, parameters and profiles.

    .. note::

        Both schemas and elements should be lower case.

    Parameters
    ----------
    template : dict
        Template containing a list of the explored elements.
    schemas : dict
        REMix schemas used to get the sets from.

    Returns
    -------
    set : set
        set containing the sets.
    """
    sets = set()
    for parameter in template["parameters"]:
        this_schema = schemas[parameter]
        this_sets = [fk["reference"]["resource"] for fk in this_schema["foreignKeys"]]
        sets.update(this_sets)
    for profile in template["profiles"]:
        this_schema = schemas[profile]
        this_sets = [fk["reference"]["resource"] for fk in this_schema["foreignKeys"]]
        sets.update(this_sets)
    for mp in template["maps"]:
        this_schema = schemas[mp]
        this_sets = [fk["reference"]["resource"] for fk in this_schema["foreignKeys"]]
        sets.update(this_sets)
    for ignored in IGNORED_SETS:
        if ignored in sets:
            sets.remove(ignored)
    sets.update(INCLUDED_SETS)
    return sets


def read_named_remix_csv(filename: str, schema: str) -> DataFrame:
    """Pipe :py:func:`remix.framework.tools.utilities.read_remix_csv` by retrieving the schema from the schema manager.

    Parameters
    ----------
    file : str
        Path to the file.
    schema : str
        Name of the schema for the tabular data.

    Returns
    -------
    df : pandas.DataFrame
        Pandas.DataFrame with the data stored from the .csv file.
    """
    if isinstance(schema, str):
        all_schemas = {k.lower(): v for k, v in sm.get_schemas().items()}
        if schema.lower() in all_schemas:
            return read_remix_csv(filename, schema=all_schemas[schema.lower()])
        else:
            raise KeyError(f"The schema {schema} was not found. Check for typos")

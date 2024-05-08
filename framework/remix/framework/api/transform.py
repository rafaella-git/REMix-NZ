import re
from pathlib import Path

from remix.framework.api.dataset import DataPackage
from remix.framework.api.dataset import DirectoryDataset
from remix.framework.api.instance import Instance
from remix.framework.backend.polars import PolarsDataset

default_transform_remix_kwargs = {
    "datadir": "./",
    "outputdir": "transform",
    "outformat": "dat",
    "mapformat": "dat",
    "profileformat": "wide",
    "indexnames": "0",
    "supportsets": "0",
    "frictionless": "0"
}

default_transform_remix_help = {
    "datadir": "Root directory of the dataset you want to transform.",
    "outputdir": "Output directory of the transformed dataset.",
    "outformat": "Format of the output dataset.",
    "mapformat": "Format of the map, if 'dat', will produce a GAMS-comaptible input, by default 'dat'",
    "profileformat": "Profile format, 'tall' profiles have a single Value column, wide profiles have a column for each timestep, by default 'wide'.",
    "indexnames": "If 0, the indeces will have no names, only works for csv, by default 0.",
    "supportsets": "If 1, data necessary for validation will be written.",
    "frictionless": "If 1, most of the options will be overridden to generate a frictionless-conform dataset."
}


def transform_dataset(datadir, outputdir, outformat="dat", mapformat="dat", profileformat="wide", supportsets=False, frictionless=False, no_sets=False):
    """Transform the input project data into other formats.

    Parameters
    ----------
    datadir : str
        Directory of the input dataset.
    outputdir : str
        Directory of the output dataset.
    outformat : str
        Format of the output dataset, by default "dat".
    mapformat : str
        Format of the map. If "dat", will produce a GAMS comaptible input, by default "dat".
    profile_format : str
        Profile format, "tall" profiles have a single Value column, "wide" profiles have a column for each timestep, by
        default "wide".
    """
    if frictionless:
        outformat = "csv"
        mapformat = "csv"
        profileformat = "tall"
        supportsets = True

    in_format = set([fil.suffix for fil in Path(datadir).iterdir()])
    in_format = [f for f in in_format if f != ""]
    if len(in_format) > 1:
        raise ValueError("The input dataset has mixed data. You need to either have only .csv or .dat files")
    elif len(in_format) == 0:
        raise ValueError("The input folder contains neither .csv nor .dat files")
    in_format = list(in_format)[0]
    # check if we are dealing with non-named csvs
    if in_format == ".csv":
        non_set_map = [fil for fil in Path(datadir).iterdir() if (("set_" not in fil.name) and ("map_" not in fil.name))]
        with open(non_set_map[0], "r") as check_file:
            first_line = check_file.readline()
            no_names = bool(re.search(r"\,\,+", first_line)) or bool(re.search(r"\,\s+\,\s+", first_line))
    else:
        no_names = False

    if (in_format == ".csv") and (not no_names) and (outformat == "csv") and (not no_sets):
        try:
            import polars as pl
            DataClass = PolarsDataset
        except ImportError:
            DataClass = DirectoryDataset
    else:
        DataClass = Instance

    data_package = DataPackage(datadir, dataset_class=DataClass)
    Path(outputdir).mkdir(exist_ok=True, parents=True)
    data_package.export(
        outputdir,
        fileformat=outformat,
        map_format=mapformat,
        profile_format=profileformat,
        support_sets=supportsets,
        no_sets=no_sets
    )

import copy
import warnings
from typing import Optional

from typer import Context
from typer import Option

from remix.framework.api.transform import default_transform_remix_help
from remix.framework.api.transform import default_transform_remix_kwargs
from remix.framework.api.transform import transform_dataset

warnings.filterwarnings("ignore")

typer_args = {
    key: Option(default_transform_remix_kwargs[key], help=default_transform_remix_help[key])
    for key in default_transform_remix_kwargs
}

def program_transform(
    gamskwargs: Context,
    datadir: Optional[str] = typer_args["datadir"],
    outputdir: Optional[str] = typer_args["outputdir"],
    outformat: Optional[str] = typer_args["outformat"],
    mapformat: Optional[str] = typer_args["mapformat"],
    profileformat: Optional[str] = typer_args["profileformat"],
    supportsets: Optional[int] = typer_args["supportsets"],
    frictionless: Optional[int] = typer_args["frictionless"],
):
    """Transform a remix datasetformat."""
    kwargs = copy.deepcopy(locals())

    for arg in gamskwargs.args:
        key, value = arg.strip("--").split("=")
        kwargs[key] = value

    _supportsets = kwargs.pop("supportsets", 0)
    _supportsets = True if _supportsets == 1 else False
    _frictionless = kwargs.pop("frictionless", 0)
    _frictionless = True if _frictionless == 1 else False

    del kwargs["gamskwargs"]

    transform_dataset(supportsets=_supportsets, frictionless=_frictionless, **kwargs)

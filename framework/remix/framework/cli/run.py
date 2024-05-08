import warnings
from typing import Optional

from typer import Context
from typer import Option

from remix.framework import run_remix
from remix.framework.api.run import default_run_remix_help
from remix.framework.api.run import default_run_remix_kwargs

warnings.filterwarnings("ignore")

typer_args = {
    key: Option(default_run_remix_kwargs[key], help=default_run_remix_help[key])
    for key in default_run_remix_kwargs
}

def program_run(
    gamskwargs: Context,
    datadir: Optional[str] = typer_args["datadir"],
    scendir: Optional[str] = typer_args["scendir"],
    resultdir: Optional[str] = typer_args["resultdir"],
    resultfile: Optional[str] = typer_args["resultfile"],
    lo: Optional[int] = typer_args["lo"],
    names: Optional[int] = typer_args["names"],
    postcalc: Optional[int] = typer_args["postcalc"],
    roundts: Optional[int] = typer_args["roundts"],
    keep: Optional[int] = typer_args["keep"],
):
    """Run a REMix model from the command line.

    Example
    -------
    In your favorite terminal, type e.g.:

    .. code::

        remix run --datadir=path/to/datadir

    For more information on the parameters available you can use

    .. code::

        remix run --help

    .. note::

        The arguments `names`, `postcalc`, `roundts`, `keep` are binary selection between 0 (False)
        and 1 (True).
    """
    kwargs = locals()

    for arg in gamskwargs.args:
        key, value = arg.strip("--").split("=")
        kwargs[key] = value

    del kwargs["gamskwargs"]

    run_remix(**kwargs)

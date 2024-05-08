from importlib.metadata import entry_points
from typing import Optional

from typer import Exit
from typer import Option
from typer import Typer
from typer.models import CommandInfo

from remix.framework import __version__

program = Typer()

CLI_PLUGINS = {
    entry_point.name: entry_point.load()
    for entry_point in entry_points().select(group="remix.plugins")
}
CONTEXT_SETTINGS = {
    "allow_extra_args": True,
    "ignore_unknown_options": True
}

for name, entry_point in CLI_PLUGINS.items():
    program.registered_commands.append(
        CommandInfo(
            name=name, callback=entry_point, context_settings=CONTEXT_SETTINGS
        )
    )

def version_callback(value: bool):
    """Function to print the REMix version on :code:`remix --version` command.

    Parameters
    ----------
    value : bool
        True if flag --version is passed with remix command.

    Raises
    ------
    typer.Exit
        Quit the program after version output.
    """
    if value:
        print(f"REMix version: {__version__}")
        raise Exit()


@program.callback()
def program_main(version: Optional[bool] = Option(None, "--version", callback=version_callback)):
    pass

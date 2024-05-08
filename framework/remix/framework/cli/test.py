import os
import warnings
from typing import Optional

from typer import Exit
from typer import Option

from remix.framework import __testingpath__
from remix.framework.api.test import default_test_remix_help
from remix.framework.api.test import default_test_remix_kwargs

warnings.filterwarnings('ignore')

typer_args = {
    key: Option(default_test_remix_kwargs[key], help=default_test_remix_help[key])
    for key in default_test_remix_kwargs
}

def program_test(
    specific: Optional[str] = typer_args["specific"],
    mode: Optional[str] = typer_args["mode"],
    junitxml: Optional[str] = typer_args["junitxml"],
    keep: Optional[int] = typer_args["keep"],
    workers: Optional[int] = typer_args["workers"],
    timelimit: Optional[int] = typer_args["timelimit"]
):
    """Run the REMix test cases. Requires [dev] dependencies."""
    kwargs = locals()
    try:
        import pytest
    except ImportError as e:
        print('## pytest is not available. Please install REMix with [dev] option.')
        raise Exit(code=1) from e
    if not os.path.exists(f"{__testingpath__}/test_instances.py"):
        print('## Test files were not found, please setup your development environment.')
        raise Exit(code=1)
    else:
        arg_list = []
        for k, v in kwargs.items():
            if os.name == "nt" and k == "workers":
                continue
            else:
                arg_list += [f"--{k}", str(v)]

        retcode = pytest.main([f"{__testingpath__}/test_instances.py"] + arg_list)
        raise Exit(code=retcode)

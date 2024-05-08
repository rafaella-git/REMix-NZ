import os
import re
import shutil
import site
import sys
import sysconfig
from subprocess import PIPE
from subprocess import Popen

from remix.framework import __gamscode__


class GAMSPathError(Exception):
    """Exeption for missing GAMS installation."""
    pass


default_run_remix_kwargs = {
    "datadir": "./",
    "resultdir": "./results",
    "scendir": "./",
    "resultfile": "remix",
    "lo": 3,
    "names": 0,
    "postcalc": 1,
    "roundts": 1,
    "keep": 0,
}

default_run_remix_help = {
    "datadir": "Root directory of scenario you want to run.",
    "scendir": "Path (relative to the datadir) to the scenario you want to run.",
    "resultdir": "Directory to write the results to.",
    "resultfile": "Name of the result file.",
    "lo": "'log option' of GAMS; ensures that the output from GAMS will be "
    "visible in the terminal (`lo=3`; `lo=4` will write out a log file).",
    "names": "Give out actual variable names in case of error.",
    "postcalc": "Run a post calculation.",
    "roundts": "Automatically round after-comma digits in large profile files where "
    "necessary to successfully read them.",
    "roundcoefs": "Automatically round profiles and remove values smaller 1e-3.",
    "keep": "Instruct GAMS to keep the scratch directory (225a) after the run is "
    "finished.",
}


def _find_gams_pythonapi():
    """Find the GAMS Python API on the user's machine."""
    # find out if gams is installed
    gams_exe = shutil.which("gams")
    if gams_exe is None:
        msg = (
            "GAMS has not been found in the user's path. Please make sure your "
            "environment variables are configured appropriately."
        )
        raise GAMSPathError(msg)
    version = _get_gams_version(gams_exe)
    sitepkg = site.getsitepackages()
    pydir = sys.base_prefix
    syspath = sys.path

    # for GAMS Python to find the environment with REMix installed
    if os.getenv("GMSPYTHONHOME") is None:
        os.environ["GMSPYTHONHOME"] = pydir

    if os.getenv("GMSPYTHONLIB") is None:
        if os.name == "posix":
            pyver = f"{sys.version_info.major}.{sys.version_info.minor}"
            if sys.platform.startswith("darwin"):
                os.environ["GMSPYTHONLIB"] = os.path.join(
                    sysconfig.get_config_var("LIBDIR"), f"libpython{pyver}.dylib"
                )
            else:
                os.environ["GMSPYTHONLIB"] = os.path.join(
                    sysconfig.get_config_var("LIBDIR"), f"libpython{pyver}.so"
                )
        else:
            pyver = f"{sys.version_info.major}{sys.version_info.minor}"
            os.environ["GMSPYTHONLIB"] = os.path.join(pydir, f"python{pyver}.dll")

    gams_pyver = f"{sys.version_info.major}{sys.version_info.minor}"
    gams_dir = os.path.dirname(gams_exe)
    gamspy = os.path.join(gams_dir, "apifiles", "Python", "gams")
    gamspy_api = os.path.join(gams_dir, "apifiles", "Python", f"api_{gams_pyver}")

    if int(version) < 42:
        if not os.path.isdir(gamspy_api):
            msg = (
                "Could not find the GAMS python API for your Python version. "
                f"The following folder could not be found: {gamspy_api}"
            )
            raise GAMSPathError(msg)

        apifiles = [gamspy_api, gamspy]

        if gamspy_api not in syspath:
            sys.path.append(gamspy_api)
    else:
        # in 42 the GAMS embedded functionality is somewhere else
        gamspy = os.path.join(gams_dir, "GMSPython", "Lib", "site-packages", "gams")
        apifiles = [gamspy]
    # for GAMS Python to find its own Python modules
    pythonpath = os.getenv("PYTHONPATH")
    if pythonpath is None:
        os.environ["PYTHONPATH"] = os.pathsep.join(sitepkg + apifiles)
    elif gamspy not in pythonpath:
        os.environ["PYTHONPATH"] = os.pathsep.join(
            sitepkg + [os.environ["PYTHONPATH"]] + apifiles
        )

    # for your environment Python to find the GAMS Python modules

    if gamspy not in syspath:
        sys.path.append(gamspy)

    return version


def _get_gams_version(executable):
    """Returns the version of the GAMS distribution installed.

    Parameters
    ----------
    executable : str, optional
        Path to the GAMS executable.

    Returns
    -------
    int
        Major version of GAMS executable.
    """
    p = Popen([f"{executable}"], stdout=PIPE)
    response = p.communicate()
    message = response[0].decode("utf-8")
    version = re.search(r"GAMS Release\s*:\s*([0-9]+\.[0-9]+\.[0-9]+)", message)[1]
    return int(version.split(".")[0])


def run_remix(**kwargs):
    """Python API to run a REMix model.

    For a full list of available parameters go to :ref:`this page <cli_label>`.
    Please be aware that only the parameters of the following list will be forwarded as GAMS core parameters to the model execution:
    a, action, docfile, dumpopt, gdx, keep, license, lo, logfile, o, optdir, output, profile

    Returns
    -------
    int
        Result of the command line :code:`run_remix.gms` call
        (:code:`subprocess.Popen.returncode`).

    Example
    -------
    You can import the :code:`run_remix` method on the
    :code:`remix.framework` level:

    >>> from remix.framework import run_remix
    """
    _find_gams_pythonapi()

    remix_parameters = default_run_remix_kwargs.copy()
    remix_parameters.update(**kwargs)
    sourcedir = os.path.join(__gamscode__, "source")
    remix_parameters.update({"sourcedir": sourcedir})

    if kwargs:
        for key, val in kwargs.items():
            remix_parameters[key] = val

    remix_call = ["gams", os.path.join(__gamscode__, "run_remix.gms")]
    gams_options = [
        "a",
        "action",
        "docfile",
        "dumpopt",
        "gdx",
        "keep",
        "license",
        "lo",
        "logfile",
        "o",
        "optdir",
        "output",
        "profile",
    ]

    for key, val in remix_parameters.items():
        if key in gams_options:
            remix_call.append(f"{key}={val}")
        else:
            remix_call.append(f"--{key}={val}")

    cwd = os.getcwd()

    try:
        process = Popen(remix_call, stdout=PIPE)
    except process.returncode != 0:
        return process.returncode

    os.chdir(cwd)

    # Poll process.stdout to show stdout live
    while True:
        output = process.stdout.readline()
        if process.poll() is not None:
            break
        if output:
            print(output.strip().decode("ascii"))

    # success
    return process.returncode

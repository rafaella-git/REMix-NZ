# -*- coding: utf-8
"""REMix energy system optimization framework"""
__version__ = "0.10.0"


import importlib.resources
import os
from pathlib import Path

__gamscode__ = os.path.join(importlib.resources.files("remix.framework"), "model")
__testingpath__ = Path(__gamscode__).parent.parent.parent / "testing"
__remixhome__ = os.path.join(os.path.expanduser("~"), ".remix")
__versionhome__ = os.path.join(os.path.expanduser("~"), ".remix", __version__)

# put imports of api commands to remix.framework top level
from .api.run import run_remix  # noqa: F401
from .api.instance import Instance
from .tools.gdx import GDXEval
from .tools.utilities import read_dat

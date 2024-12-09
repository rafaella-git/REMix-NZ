# SPDX-FileCopyrightText: Copyright (c) 2023 German Aerospace Center (DLR)
# SPDX-License-Identifier: BSD-3-Clause
import os
import shutil
from pathlib import Path

import pytest

from remix.framework import __versionhome__
from remix.framework.api.run import run_remix

LATEST = Path(__versionhome__) / "schema/json"
MINIMAL_LP_PATH = Path(r"testing/instances/minimal_lp/data")


def pytest_addoption(parser):
    parser.addoption("--specific", action="store", default="")
    parser.addoption("--timelimit", action="store", default="7200")
    parser.addoption("--mode", action="store", default=".")
    parser.addoption("--keep", action="store", default="0")


def clean_metadata():
    metadata_exists = LATEST.exists()
    metadata_isempty = not any(LATEST.iterdir()) if metadata_exists else True
    if not metadata_isempty:
        for filename in os.listdir(LATEST):
            file_path = os.path.join(LATEST, filename)
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)


@pytest.fixture
def metadata_cleaner():
    return clean_metadata


@pytest.fixture
def example_gdx(tmp_path_factory):
    results = tmp_path_factory.mktemp("results")
    code = run_remix(
        datadir=MINIMAL_LP_PATH, resultdir=results, resultfile="remix", metadata="0"
    )
    gdx = results / "remix.gdx"
    assert code == 0
    assert gdx.exists
    return str(gdx)


clean_metadata()

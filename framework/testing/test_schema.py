import os
import shutil
from pathlib import Path

import pytest

from remix.framework import __remixhome__
from remix.framework import __versionhome__
from remix.framework.schema.parser import main as parse
from remix.framework.schema.templates import read_named_remix_csv
from remix.framework.schema.templates import sm
from remix.framework.tools.utilities import read_remix_csv

LATEST = Path(__versionhome__) / "schema/json"
LOCAL_PROFILE = Path(__remixhome__).joinpath("specifications").joinpath("table-schema.json")

class TestManager:
    def test_validate_schema_checker(self):
        if LOCAL_PROFILE.exists():
            if os.path.isfile(LOCAL_PROFILE) or os.path.islink(LOCAL_PROFILE):
                os.unlink(LOCAL_PROFILE)
            elif os.path.isdir(LOCAL_PROFILE):
                shutil.rmtree(LOCAL_PROFILE)
        sm.check_latest()

    def test_get_tabular_package_template(self):
        package_template = sm.get_tabular_package_template()
        for key in ["name", "profile", "resources"]:
            assert key in package_template.keys()


class TestTemplates:
    @pytest.mark.skip(reason="Not implemented")
    def test_default_dataframes(self):
        raise NotImplementedError

    @pytest.mark.skip(reason="Not implemented")
    def test_distribute_dictionary(self):
        raise NotImplementedError

    @pytest.mark.skip(reason="Not implemented")
    def test_empty_dataframe(self):
        raise NotImplementedError

    @pytest.mark.skip(reason="Not implemented")
    def test_empty_profile(self):
        raise NotImplementedError

    def test_read_csv_from_schema(self):
        converter_techparam = "testing/input/named_csv_input/data/converter_techparam.csv"
        converter_techparam_schema ={k.lower(): v for k, v in sm.get_schemas().items()}["converter_techparam"]
        df = read_remix_csv(converter_techparam, converter_techparam_schema)
        index_columns = [fk["fields"][0] for fk in converter_techparam_schema["foreignKeys"]]
        assert set(index_columns ) == set(df.index.names), f"Invalid index names {set(df.index.names)}"

    def test_read_named_csv_from_string(self):
        converter_techparam = "testing/input/named_csv_input/data/converter_techparam.csv"
        converter_techparam_schema ={k.lower(): v for k, v in sm.get_schemas().items()}["converter_techparam"]
        df = read_named_remix_csv(converter_techparam, "converter_techparam")
        index_columns = [fk["fields"][0] for fk in converter_techparam_schema["foreignKeys"]]
        assert set(index_columns ) == set(df.index.names), f"Invalid index names {set(df.index.names)}"

    def test_read_named_csv_fail(self):
        converter_techparam = "testing/input/named_csv_input/data/converter_techparam.csv"
        with pytest.raises(KeyError):
            df = read_named_remix_csv(converter_techparam, "converter_techsparams")

    def test_read_csv_fail(self):
        converter_techparam = "testing/input/named_csv_input/data/converter_techparam.csv"
        with pytest.raises(TypeError):
            df = read_remix_csv(converter_techparam, "converter_techsparams")


class TestParser:
    def test_main(self, metadata_cleaner):
        metadata_cleaner()
        LATEST = Path(__versionhome__) / "schema/input/json"
        parse()
        METADATA_EXISTS = LATEST.exists()
        METADATA_ISEMPTY = not any(LATEST.iterdir()) if METADATA_EXISTS else True
        assert not METADATA_ISEMPTY

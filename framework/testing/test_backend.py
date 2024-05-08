from pathlib import Path

import pytest

from remix.framework.api.instance import Instance
from remix.framework.backend.excel import ExcelInstance
from remix.framework.backend.excel import load_workbook
from remix.framework.backend.polars import PolarsDataset
from remix.framework.backend.polars import pl

EXAMPLE_DATA_PATH = Path(r"testing/instances/minimal_lp/data")
NAMED_MODEL_INSTANCE = Path(r"testing/input/named_csv_input/data")

HAS_POLARS = True if pl is not None else False
HAS_OPENPYXL = True if load_workbook is not None else False


@pytest.fixture
def example_dataframes():
    instance = Instance.from_path(path=EXAMPLE_DATA_PATH, inherit=True, index_names=True)
    dataframes = instance.dataframes
    return dataframes


@pytest.fixture
def temporal_excel_file(tmp_path, example_dataframes):
    excel_instance = ExcelInstance.from_dataframes(example_dataframes, workbook=f"{tmp_path}/example.xlsx")
    excel_instance.write(fileformat="xlsx")
    assert Path(f"{tmp_path}/example.xlsx").exists()
    return f"{tmp_path}/example.xlsx"


@pytest.mark.skipif(not HAS_POLARS, reason="Requires polars extra")
class TestPolarsBackend:
    def test_polars_csv_input_output(self, tmp_path):
        scenario = PolarsDataset.from_path(
            NAMED_MODEL_INSTANCE, parent_files=None, version="latest", name="."
        )
        scenario.inherit()
        scenario.write(output_path=tmp_path)
        assert len(list(tmp_path.iterdir())) == 21

    def test_polars_csv_output_tall(self, tmp_path):
        import polars as pl
        scenario = PolarsDataset.from_path(
            NAMED_MODEL_INSTANCE, parent_files=None, version="latest", name="."
        )
        scenario.inherit()
        scenario.write(output_path=tmp_path, profile_format="tall")
        assert len(list(tmp_path.iterdir())) == 21
        converter_activity_profile = pl.read_csv(tmp_path.joinpath("converter_activityprofile.csv"))
        assert converter_activity_profile.columns == ['nodesData', 'years', 'converter_techs', 'profileTypes', 'timeData', 'Value']


@pytest.mark.skipif(not HAS_OPENPYXL, reason="Requires excel extra")
class TestExcelBackend:
    def test_excel_input_output(self, temporal_excel_file, tmp_path):
        instance = ExcelInstance(workbook=temporal_excel_file, datadir=tmp_path)
        instance.write(fileformat="csv", infer=True)
        assert len(list(tmp_path.iterdir())) == 38

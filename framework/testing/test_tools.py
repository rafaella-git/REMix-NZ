import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from remix.framework import __testingpath__
from remix.framework.tools.gdx import GDXEval
from remix.framework.tools.gdx import fill_missing_values_in_timeseries
from remix.framework.tools.gdx import read_gdx_symbol
from remix.framework.tools.gdx import select_from_dataframe
from remix.framework.tools.gdx import write_gdx_symbol_as_csv
from remix.framework.tools.inheritance import main_inheritance
from remix.framework.tools.utilities import to_dat

IS_GITLAB = os.getenv("GITLAB_CI", default="false").lower() in ["true", "1", "t", "y", "yes"]


class TestInheritance:
    def test_main_inheritance(self, tmp_path_factory):

        scr_dir = tmp_path_factory.mktemp("tmp")
        base_path = os.path.join(__testingpath__, "instances", "minimal_lp", "data")
        scen_dir = "."

        main_inheritance(scr_dir, base_path, scen_dir, roundts=1, _logger=print)
        assert len(list(Path(scr_dir).iterdir())) == 39


@pytest.mark.skip(reason="Not implemented")
class TestMGA:
    def test_int_sin_m():
        raise NotImplementedError

    def test_primes():
        raise NotImplementedError

    def test_inverse_increasing():
        raise NotImplementedError

    def test_uniform_hypersphere():
        raise NotImplementedError


@pytest.mark.skip(reason="Not implemented")
class TestTransfer:
    def test_build_cmat():
        raise NotImplementedError


# @pytest.mark.skip(reason="Not implemented")
class TestUtilities:
    # def test_merge_dfs():
    #     raise NotImplementedError

    # def test_concat_ts():
    #     raise NotImplementedError

    # def test_read_dat():
    #     raise NotImplementedError

    def test_to_dat(self, tmp_path_factory):
        temporary = tmp_path_factory.mktemp("tmp")
        df1 = pd.DataFrame({"A": [0.12345678, 0.12345691, 0.12345, 0.1234, 0.123, 0.12, 0.1, 0]})
        to_dat(df1, temporary / "test1.dat")
        df2 = pd.DataFrame({"A": [0.1234567891234, 0.123456, 0.12345, 0.1234, 158888880.123, 0.12, 0.1, np.inf]})
        to_dat(df2, temporary / "test2.dat")
        # TODO: Do a reasonable test of this feature

    # def test_get_limited_df():
    #     raise NotImplementedError

    # def test_read_csv():
    #     raise NotImplementedError

    # def test_read_file():
    #     raise NotImplementedError

    # def test_read_timeseries_range():
    #     raise NotImplementedError


@pytest.mark.skipif(
    IS_GITLAB,
    reason="GitLab has weird interactions with pytest temporary directories, which is why TestGDX won't run there.",
)
class TestGDX:
    def test_read_gdx_symbol(self, example_gdx):
        indicator_accounting = read_gdx_symbol(example_gdx, "indicator_accounting")
        assert "value" in indicator_accounting.columns

    def test_read_gdx_symbol_with_sel_dict(self, example_gdx):
        sel_dict = {
            "accNodesModel": ["CH_IT"],
            "accYears": ["2020"],
            "indicator": ["EmissionTax", "Invest"],
        }
        indicator_accounting = read_gdx_symbol(
            example_gdx, "indicator_accounting", sel_dict=sel_dict
        )
        assert (
            len(indicator_accounting.values) == 2
        ), "Something went wrong with the GDX from minimal-lp"

    def test_select_from_dataframe(self):
        df = pd.DataFrame(
            {
                "Fruits": ["Oranges", "Apples", "Apples", "Tomatoes"],
                "Color": ["orange", "red", "yellow", "red"],
                "Amount": [30, 10, 20, 12],
            }
        ).set_index(["Fruits", "Color"])
        sel_dict = {"Fruits": ["Apples"], "Color": ["yellow"]}
        selected = select_from_dataframe(df, sel_dict=sel_dict)
        assert selected.iloc[0, 0] == 20, "select_from_dataframe is broken"

    def test_write_gdx_symbol_as_csv(self, example_gdx, tmp_path):
        write_gdx_symbol_as_csv(example_gdx, "indicator_accounting", tmp_path)
        assert (
            tmp_path / "indicator_accounting.csv"
        ).exists(), "the normal symbol was not exported correctly"

        write_gdx_symbol_as_csv(example_gdx, "accYears", tmp_path)
        assert (
            tmp_path / "set_accYears.csv"
        ).exists(), "the set symbol was not exported correctly"

        write_gdx_symbol_as_csv(example_gdx, "map_linksModel", tmp_path)
        assert (
            tmp_path / "map_linksModel.csv"
        ).exists(), "the map symbol was not exported correctly"

    def test_fill_missing_values_in_timeseries(self, example_gdx):
        commodity_balance = read_gdx_symbol(example_gdx, "commodity_balance")
        timeModeltoCalc = read_gdx_symbol(example_gdx, "timeModeltoCalc")
        filled = fill_missing_values_in_timeseries(commodity_balance, timeModeltoCalc)
        pv_example = filled.loc[:, "CH_IT", "2020", "PV", "Elec"]
        assert len(pv_example.index) == 8760, "The index was not filled properly"

    def test_fill_missing_values_in_timeseries_error(self, example_gdx):
        commodity_balance = read_gdx_symbol(example_gdx, "commodity_balance_annual")
        timeModeltoCalc = read_gdx_symbol(example_gdx, "timeModeltoCalc")
        with pytest.raises(ValueError):
            filled = fill_missing_values_in_timeseries(
                commodity_balance, timeModeltoCalc
            )

    def test_gdx_eval_init(self, example_gdx, tmp_path):
        gd = GDXEval(example_gdx)
        gd["commodity_balance_annual"]
        assert (
            len(gd.get_metadata(example_gdx)) == 39
        ), "The amount of symbols in the test gdx differs."
        gd.export(tmp_path)
        assert (
            len(list(tmp_path.iterdir())) == 36
        ), "The export csv function did not work properly."

    def test_gdx_eval_init_timeModel(self, example_gdx):
        gd = GDXEval(example_gdx)
        idx_len = len(gd["commodity_balance"].index.levels[0].unique())
        idx_last_digit = list(set(gd["commodity_balance"].index.get_level_values("timeModel")))[-1]
        assert (
            idx_len == idx_last_digit
        ), "The index is not sorted correctly."

    def test_gdx_eval_init_from_dict(self, example_gdx, tmp_path):
        gd = GDXEval({"A": example_gdx, "B": example_gdx})
        assert (
            "scenario" in gd["commodity_balance_annual"].index.names
        ), "The multi scenario reader did not work properly."
        gd.export(tmp_path, stacked=False)
        with pytest.raises(NotImplementedError):
            gd.export(tmp_path, stacked=True)

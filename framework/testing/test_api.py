import os
import re
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
from pandas import DataFrame
from pandas import IndexSlice as idx
from pandas import MultiIndex
from pandas import read_csv
from pandas.errors import IntCastingNaNError
from pandas.testing import assert_frame_equal

from remix.framework.api.dataset import DataPackage
from remix.framework.api.dataset import DirectoryDataset
from remix.framework.api.instance import Instance
from remix.framework.api.instance import container_factory
from remix.framework.api.run import run_remix
from remix.framework.api.transform import transform_dataset
from remix.framework.schema.templates import default_dataframes
from remix.framework.tools.utilities import build_timeseries_range
from remix.framework.tools.utilities import read_dat
from remix.framework.tools.utilities import read_remix_csv

MINIMAL_LP_PATH = Path(r"testing/instances/minimal_lp/data")
MINIMAL_LP_SUBSCENARIOS_PATH = Path(r"testing/instances/minimal_lp_with_scenarios/data")
MINIMAL_LP_CSV_PATH = Path(r"testing/input/minimal_lp_csv/data")
MININAL_LP_NAMED_CSV_PATH = Path(r"testing/input/named_csv_input/data")


@pytest.fixture
def example_dataframes():
    sets = {
        "set_converter_techs": ["tech_a", "tech_b"],
        "set_nodesdata": ["DE_BW"],
        "set_years": ["2050"],
        "set_indicators": ["SystemCost"],
        "set_accnodes": ["Europe"],
    }
    index = MultiIndex.from_product([["DE_BW"], ["2050"], ["tech_a"]])
    accounting_indicatorbounds = DataFrame(
        {"discount": [0.03], "endyear": [1], "obj": [1]},
        index=MultiIndex.from_product([["global"], ["horizon"], ["SystemCost"]]),
    )
    parameters = {
        "converter_capacityparam": DataFrame(
            {"unitsLowerLimit": 0, "unitsUpperLimit": 1}, index=index
        ),
        "accounting_indicatorbounds": accounting_indicatorbounds,
    }
    converter_activityprofile = DataFrame(
        columns=build_timeseries_range(),
        index=MultiIndex.from_product([["DE_BW"], ["2050"], ["tech_a"], ["upper"]]),
        dtype=float,
    ).fillna(0)
    timeseries = {"converter_activityprofile": converter_activityprofile}
    map_aggregatenodesmodel = DataFrame(
        index=MultiIndex.from_tuples([("DE_BW", "DE"), ("DE_BA", "DE")])
    )
    map_accnodes = DataFrame(
        index=MultiIndex.from_tuples([("DE_BW", "Europe"), ("DE_BA", "Europe")])
    )
    maps = {
        "map_aggregatenodesmodel": map_aggregatenodesmodel,
        "map_accnodes": map_accnodes,
    }
    dataframes = {
        "sets": sets,
        "parameters": parameters,
        "timeseries": timeseries,
        "maps": maps,
    }
    return dataframes


@pytest.fixture
def example_dataframes_link_types(example_dataframes):
    transfer_coefperlength = DataFrame(
        columns=["coefPerLength"],
        index=MultiIndex.from_product(
            [["tech_a"], ["2050"], ["commodity_a"], ["sea", "land"]]
        ),
    ).fillna(0)
    example_dataframes["parameters"]["transfer_coefperlength"] = transfer_coefperlength

    return example_dataframes


@pytest.fixture
def example_multiscenario_datfiles(tmp_path_factory, example_dataframes):

    dfs = deepcopy(example_dataframes)
    scen = DirectoryDataset.from_dataframes(dfs)
    temporary_sub_dat_files = tmp_path_factory.mktemp("tmp")
    scen.write(temporary_sub_dat_files, fileformat="dat")

    sets = example_dataframes["sets"]
    set_nodesdata = deepcopy(sets["set_nodesdata"])
    set_nodesdata.extend(["FR", "IT"])
    index = MultiIndex.from_product([set_nodesdata, ["2050"], ["tech_a"]])
    map_aggregatenodesmodel = DataFrame(
        index=MultiIndex.from_tuples([("DE_BW", "DE_BW"), ("DE_BA", "DE_BA")])
    )
    subscenario = {
        "sets": {
            "set_nodesdata": set_nodesdata,
        },
        "parameters": {
            "converter_capacityparam": DataFrame(
                {"unitsLowerLimit": 0, "unitsUpperLimit": 1}, index=index
            ),
        },
        "maps": {"map_aggregatenodesmodel": map_aggregatenodesmodel},
        "timeseries": {},
    }
    scen_sub = DirectoryDataset.from_dataframes(subscenario)
    sub = temporary_sub_dat_files / "sub"
    sub.mkdir()
    scen_sub.write(sub, fileformat="dat")
    return temporary_sub_dat_files


@pytest.fixture
def example_three_level_multiscenario(tmp_path_factory, example_dataframes):

    dfs = deepcopy(example_dataframes)
    scen = DirectoryDataset.from_dataframes(dfs)
    temporary_sub_dat_files = tmp_path_factory.mktemp("tmp")
    scen.write(temporary_sub_dat_files, fileformat="dat")

    sets = example_dataframes["sets"]
    set_nodesdata = deepcopy(sets["set_nodesdata"])
    set_nodesdata.extend(["FR", "IT"])
    subscenario = {
        "sets": {
            "set_nodesdata": set_nodesdata,
        },
        "parameters": {},
        "maps": {},
        "timeseries": {},
    }
    scen_sub = DirectoryDataset.from_dataframes(subscenario)
    sub = temporary_sub_dat_files / "sub"
    sub.mkdir()
    scen_sub.write(sub, fileformat="dat")

    set_nodesdata.extend(["ES"])
    subscenario = {
        "sets": {
            "set_nodesdata": set_nodesdata,
        },
        "parameters": {},
        "maps": {},
        "timeseries": {},
    }
    scen_sub_sub = DirectoryDataset.from_dataframes(subscenario)
    sub_sub = temporary_sub_dat_files / "sub" / "sub"
    sub_sub.mkdir()
    scen_sub_sub.write(sub_sub, fileformat="dat")
    return temporary_sub_dat_files


@pytest.fixture
def example_multiscenario_csvfiles(tmp_path_factory, example_dataframes):

    dfs = deepcopy(example_dataframes)
    scen = DirectoryDataset.from_dataframes(dfs)
    temporary_sub_dat_files = tmp_path_factory.mktemp("tmp")
    scen.write(temporary_sub_dat_files, fileformat="csv")

    sets = example_dataframes["sets"]
    set_nodesdata = deepcopy(sets["set_nodesdata"])
    set_nodesdata.extend(["FR", "IT"])
    map_aggregatenodesmodel = DataFrame(
        index=MultiIndex.from_tuples([("DE_BW", "DE_BW"), ("DE_BA", "DE_BA")])
    )
    subscenario = {
        "sets": {
            "set_nodesdata": set_nodesdata,
        },
        "parameters": {},
        "maps": {"map_aggregatenodesmodel": map_aggregatenodesmodel},
        "timeseries": {},
    }
    scen_sub = DirectoryDataset.from_dataframes(subscenario)
    sub = temporary_sub_dat_files / "sub"
    sub.mkdir()
    scen_sub.write(sub, fileformat="csv")
    return temporary_sub_dat_files


@pytest.fixture
def temporary_dat_files(tmp_path_factory, example_dataframes):
    temporary_dat_files = tmp_path_factory.mktemp("tmp")
    dfs = deepcopy(example_dataframes)
    scen = DirectoryDataset.from_dataframes(dfs)
    scen.write(temporary_dat_files, fileformat="dat")
    return temporary_dat_files


@pytest.fixture
def temporary_csv_files(tmp_path_factory, example_dataframes):
    temporary_csv_files = tmp_path_factory.mktemp("tmp")
    dfs = deepcopy(example_dataframes)
    scen = DirectoryDataset.from_dataframes(dfs)
    scen.write(temporary_csv_files, fileformat="csv")
    return temporary_csv_files


@pytest.fixture
def csv_tall_normal_maps_path(tmp_path_factory):
    this_path = tmp_path_factory.mktemp("data")
    transform_dataset(
        MINIMAL_LP_PATH,
        this_path,
        outformat="csv",
        mapformat="csv",
        profileformat="tall",
    )
    with open(this_path.joinpath("map_aggregatenodesmodel.csv"), "r") as check_file:
        first_line = check_file.readline()
        assert set(first_line.strip().split(",")) == {"nodesData", "nodesModel"}
    with open(this_path.joinpath("converter_activityprofile.csv"), "r") as check_file:
        first_line = check_file.readline()
        assert set(first_line.strip().split(",")) == {
            "nodesData",
            "timeData",
            "profileTypes",
            "converter_techs",
            "Value",
            "years",
        }
    return this_path


def _count_non_empty_files_in_dir(dir):
    return len(
        [
            f
            for f in os.listdir(dir)
            if os.path.isfile(os.path.join(dir, f))
            and os.stat(os.path.join(dir, f)).st_size > 0
        ]
    )


def _get_files_by_prefix_exclusion(path, non_prefixes):
    return [
        fil
        for fil in Path(path).iterdir()
        if all([prefix not in fil.name for prefix in non_prefixes]) and not fil.is_dir()
    ]


class TestDataset:
    def test_dataframe_constructor(self, example_dataframes):
        scenario = DirectoryDataset.from_dataframes(example_dataframes)

        assert (
            scenario.ready
        ), "Something went wrong during the from dataframes construction"
        assert isinstance(
            scenario, DirectoryDataset
        ), "Something went wrong during the from dataframes construction"

        assert scenario.dataframes["sets"]["set_converter_techs"] == [
            "tech_a",
            "tech_b",
        ]
        assert scenario.dataframes["sets"]["set_nodesdata"] == ["DE_BW"]
        assert scenario.dataframes["sets"]["set_years"] == ["2050"]
        assert scenario.dataframes["sets"]["set_indicators"] == ["SystemCost"]

        assert_frame_equal(
            scenario.dataframes["parameters"]["accounting_indicatorbounds"],
            example_dataframes["parameters"]["accounting_indicatorbounds"],
        )
        assert_frame_equal(
            scenario.dataframes["parameters"]["converter_capacityparam"],
            example_dataframes["parameters"]["converter_capacityparam"],
        )

    def test_flat_dataframes(self, example_dataframes):
        only_values = [v for v in example_dataframes.values()]
        flat_input = {}
        for vals in only_values:
            flat_input.update(vals)
        scenario = DirectoryDataset.from_dataframes(flat_input)
        assert (
            scenario.ready
        ), "Something went wrong during the from dataframes construction"

    def test_path_constructor(self):
        scenario = DirectoryDataset.from_path(
            MINIMAL_LP_PATH, parent_files=None, version="latest", name="."
        )
        scenario.inherit()
        assert scenario.ready, "Something went wrong during the from path construction"

        assert "DE" in scenario.dataframes["sets"]["set_nodesdata"]
        assert scenario.dataframes["sets"]["set_years"] == ["2020"]

    def test_write(self, tmp_path):
        scenario = DirectoryDataset.from_path(
            MINIMAL_LP_PATH, parent_files=None, version="latest", name="."
        )
        scenario.inherit()
        tp = tmp_path / "data"
        tp.mkdir()
        scenario.write(tp, inspect=False, fileformat="dat")
        assert _count_non_empty_files_in_dir(tp) == _count_non_empty_files_in_dir(
            MINIMAL_LP_PATH
        )

    def test_write_with_support_sets(self, tmp_path):
        scenario = DirectoryDataset.from_path(
            MINIMAL_LP_PATH, parent_files=None, version="latest", name="."
        )
        scenario.inherit()
        tp = tmp_path / "data"
        tp.mkdir()
        scenario.write(tp, inspect=False, fileformat="dat", support_sets=True)
        files = [Path(f).name for f in tp.iterdir()]
        assert _count_non_empty_files_in_dir(tp) == 42
        assert "set_profiletypes.dat" in files, "set_profiletypes is missing"
        assert "set_timedata.dat" in files, "set_timedata is missing"

    def test_scenario_find_files(self, temporary_dat_files):
        scen = DirectoryDataset.from_path(temporary_dat_files)
        scen._find_files()
        file_names = [Path(f).name for f in scen.files]
        assert set(file_names) == {
            "set_years.dat",
            "map_aggregatenodesmodel.dat",
            "set_nodesdata.dat",
            "set_accnodes.dat",
            "map_accnodes.dat",
            "set_converter_techs.dat",
            "accounting_indicatorbounds.dat",
            "converter_capacityparam.dat",
            "set_indicators.dat",
            "converter_activityprofile.dat",
        }, "The files found in the scenario do not match."

    def test_scenario_arrange_files(self, temporary_dat_files):
        scen = DirectoryDataset.from_path(temporary_dat_files)
        scen._arrange_files()
        scen.arranged_files
        assert all(
            [len(v) == 1 for v in scen.arranged_files.values()]
        ), "There was unexpected duplicates in arranged files."

    def test_split_parameters_sets(self, temporary_dat_files):
        scen = DirectoryDataset.from_path(temporary_dat_files)
        scen._arrange_files()
        scen._split_parameters_sets()
        assert all(
            [len(v) == 1 for v in scen.parameters.values()]
        ), "There was unexpected duplicates in parameter files."
        assert all(
            [len(v) == 1 for v in scen.sets.values()]
        ), "There was unexpected duplicates in set files."
        assert all(
            [len(v) == 1 for v in scen.timeseries.values()]
        ), "There was unexpected duplicates in timeseries files."
        assert all(
            [len(v) == 1 for v in scen.maps.values()]
        ), "There was unexpected duplicates in maps files."

    def test_search_duplicates_raise(self, temporary_dat_files):
        scen = DirectoryDataset.from_path(temporary_dat_files)
        scen._arrange_files()
        scen._split_parameters_sets()
        # Adding some duplicate
        scen.files.append(
            scen.parameters["converter_capacityparam.dat"][0].replace(
                "converter_capacityparam", "converter_capacityParam"
            )
        )
        with pytest.raises(RuntimeError) as exc_info:
            scen._search_duplicates()
        assert (
            exc_info.value.args[0] == "Duplicate files found in the input"
        ), "The duplicate file error was not raised properly"

    def test_search_duplicates_pass(self, temporary_dat_files):
        scen = DirectoryDataset.from_path(temporary_dat_files)
        scen._arrange_files()
        scen._split_parameters_sets()

    def test_inherit_dataframes(self, temporary_dat_files, example_dataframes):
        scen = DirectoryDataset.from_path(temporary_dat_files)
        scen._arrange_files()
        scen._split_parameters_sets()
        scen._inherit_dataframes()
        assert (
            scen.dataframes.keys() == example_dataframes.keys()
        ), "The produced dataframe keys don't match"
        for k, v in scen.dataframes.items():
            assert (
                v.keys() == example_dataframes[k].keys()
            ), f"The contents of {k} don't match"
            for subk, subv in v.items():
                if isinstance(subv, list):
                    assert set(subv) == set(example_dataframes[k][subk])
                elif isinstance(subv, DataFrame):
                    assert set(subv.columns) == set(example_dataframes[k][subk].columns)

    def test_inherit_dataframes_csv(self, temporary_csv_files, example_dataframes):
        scen = DirectoryDataset.from_path(temporary_csv_files)
        scen._arrange_files()
        scen._split_parameters_sets()
        scen._inherit_dataframes()
        assert (
            scen.dataframes.keys() == example_dataframes.keys()
        ), "The produced dataframe keys don't match"
        for k, v in scen.dataframes.items():
            assert (
                v.keys() == example_dataframes[k].keys()
            ), f"The contents of {k} don't match"
            for subk, subv in v.items():
                if isinstance(subv, list):
                    assert set(subv) == set(example_dataframes[k][subk])
                elif isinstance(subv, DataFrame):
                    assert set(subv.columns) == set(example_dataframes[k][subk].columns)

    def test_write_csv_roundtrip(self, tmp_path):
        scenario = DirectoryDataset.from_path(
            MINIMAL_LP_PATH, parent_files=None, version="latest", name="."
        )
        scenario.inherit()
        tp = tmp_path / "data"
        tp.mkdir()
        scenario.write(tp, inspect=False, fileformat="csv")
        assert _count_non_empty_files_in_dir(tp) == _count_non_empty_files_in_dir(
            MINIMAL_LP_PATH
        )

        scenario_csv = DirectoryDataset.from_path(
            tp, parent_files=None, version="latest", name="."
        )
        scenario_csv.inherit()

        for attribute in ["maps", "parameters", "timeseries", "sets"]:
            for key, value in scenario.dataframes[attribute].items():
                if attribute == "sets":
                    assert set(value) == set(scenario_csv.dataframes[attribute][key])
                else:
                    if len(value.index) > 0:
                        assert set(value.index) == set(
                            scenario_csv.dataframes[attribute][key].index
                        )
                    if len(value.columns) > 0:
                        assert set(value.columns) == set(
                            scenario_csv.dataframes[attribute][key].columns
                        )

    def test_named_csvs_input(self):
        scenario = DirectoryDataset.from_path(
            "testing/input/named_csv_input/data",
            parent_files=None,
            version="latest",
            name=".",
        )
        scenario.inherit()

    def test_write_tall_timeseries_with_name(self, tmp_path):
        scenario = DirectoryDataset.from_path(
            "testing/input/named_csv_input/data",
            parent_files=None,
            version="latest",
            name=".",
        )
        scenario.inherit()
        tp = tmp_path / "data"
        tp.mkdir()
        scenario.write(tp, inspect=False, fileformat="csv", profile_format="tall")
        scenario_tall = DirectoryDataset.from_path(
            tp, parent_files=None, version="latest", name=".", inherit=True
        )
        scenario_tall.write(
            output_path=tmp_path, profile_format="tall", fileformat="csv"
        )
        converter_activity_profile = read_csv(
            tmp_path.joinpath("converter_activityprofile.csv")
        )
        assert set(converter_activity_profile.columns) == set(
            [
                "nodesData",
                "years",
                "converter_techs",
                "profileTypes",
                "timeData",
                "Value",
            ]
        )

    def test_write_csv_from_list(self, tmp_path, temporary_dat_files):
        scenario = DirectoryDataset.from_path(
            temporary_dat_files, parent_files=None, version="latest", name="."
        )
        scenario.inherit()
        tp = tmp_path / "data"
        tp.mkdir()
        index = MultiIndex.from_product([["FR"], ["2050"], ["tech_a"]])
        old_param = scenario.dataframes["parameters"]["converter_capacityparam"]
        new_param = [
            old_param,
            DataFrame({"unitsLowerLimit": 0, "unitsUpperLimit": 1}, index=index),
        ]
        scenario.dataframes["parameters"]["converter_capacityparam"] = new_param

        old_map = scenario.dataframes["maps"]["map_aggregatenodesmodel"]
        new_map = [old_map, DataFrame(index=MultiIndex.from_tuples([("FR", "FR")]))]
        scenario.dataframes["maps"]["map_aggregatenodesmodel"] = new_map

        scenario.write(tp, inspect=False, fileformat="csv")
        assert len(list(tp.iterdir())) == 10
        converter_capacityparam = read_remix_csv(f"{tp}/converter_capacityparam.csv")
        map_aggregatenodesmodel = read_remix_csv(f"{tp}/map_aggregatenodesmodel.csv")

        assert tuple(converter_capacityparam.index) == (
            ("DE_BW", "2050", "tech_a"),
            ("FR", "2050", "tech_a"),
        ), "The parameters were not built correctly"
        assert tuple(map_aggregatenodesmodel.index) == (
            ("DE_BA", "DE"),
            ("DE_BW", "DE"),
            ("FR", "FR"),
        ), "The maps were not built correctly"

    def test_multiscenario_dat_package(self, example_multiscenario_datfiles, tmp_path):
        dp = DataPackage(example_multiscenario_datfiles)
        base_path = tmp_path
        sub_path = tmp_path / "sub"
        sub_path.mkdir()
        dp.inherit_scenario(".", base_path, fileformat="dat")
        dp.inherit_scenario("sub", sub_path, fileformat="dat")

        assert len(list(sub_path.iterdir())) == 10

        with open(sub_path / "set_nodesdata.dat", "r") as subset:
            data = subset.read()
            sets = data.strip().split("\n")

        df = read_dat(sub_path / "converter_capacityparam.dat")

        assert set(sets) == {"DE_BW", "FR", "IT"}
        assert tuple(df.index) == (
            ("DE_BW", "2050", "tech_a"),
            ("FR", "2050", "tech_a"),
            ("IT", "2050", "tech_a"),
        )

        aggregate_nodesmodel = read_dat(sub_path / "map_aggregatenodesmodel.dat")

        assert tuple(aggregate_nodesmodel.index) == (
            ("DE_BW", "DE_BW"),
            ("DE_BA", "DE_BA"),
        )

    def test_threescenario_dat_package(
        self, example_three_level_multiscenario, tmp_path
    ):
        dp = DataPackage(example_three_level_multiscenario)
        base_path = tmp_path
        sub_path = tmp_path / "sub"
        sub_path.mkdir()
        sub_sub_path = tmp_path / "sub_sub"
        sub_sub_path.mkdir()

        dp.inherit_scenario(".", base_path, fileformat="dat")
        dp.inherit_scenario("sub", sub_path, fileformat="dat")
        dp.inherit_scenario("sub/sub", sub_sub_path, fileformat="dat")

        assert len(list(sub_path.iterdir())) == 10
        assert len(list(sub_sub_path.iterdir())) == 10

        with open(sub_path / "set_nodesdata.dat", "r") as subset:
            data = subset.read()
            sets = data.strip().split("\n")

        assert set(sets) == {"DE_BW", "FR", "IT"}

        with open(sub_sub_path / "set_nodesdata.dat", "r") as subset:
            data = subset.read()
            sets = data.strip().split("\n")

        assert set(sets) == {"DE_BW", "FR", "IT", "ES"}

    def test_inherit_all_package(self, example_three_level_multiscenario, tmp_path):
        dp = DataPackage(example_three_level_multiscenario)
        dp.inherit_all(tmp_path)

        assert len(list(tmp_path.iterdir())) == 3
        for sub_dir in tmp_path.iterdir():
            assert len(list(Path(sub_dir).iterdir())) == 10

    def test_multiscenario_csv_package(self, example_multiscenario_csvfiles, tmp_path):
        dp = DataPackage(example_multiscenario_csvfiles)
        base_path = tmp_path
        sub_path = tmp_path / "sub"
        sub_path.mkdir()
        dp.inherit_scenario(".", base_path, fileformat="dat")
        dp.inherit_scenario("sub", sub_path, fileformat="dat")

        assert len(list(sub_path.iterdir())) == 10

        with open(sub_path / "set_nodesdata.dat", "r") as subset:
            data = subset.read()
            sets = data.strip().split("\n")

        assert set(sets) == {"DE_BW", "FR", "IT"}

    def test_load_package_from_string(self, temporary_dat_files):
        dp = DataPackage(str(temporary_dat_files))
        dataframes = dp.inherit_scenario(".")

        assert set(dataframes["sets"]["set_nodesdata"]) == {"DE_BW"}

    def test_inherit_scenario_at_init(self, temporary_dat_files):
        scenario = DirectoryDataset.from_path(
            temporary_dat_files,
            parent_files=None,
            version="latest",
            name=".",
            inherit=True,
        )
        assert set(scenario.dataframes["sets"]["set_nodesdata"]) == {"DE_BW"}

    def test_inherit_scenario_from_str(self, temporary_dat_files):
        scenario = DirectoryDataset.from_path(
            str(temporary_dat_files),
            parent_files=None,
            version="latest",
            name=".",
            inherit=True,
        )
        assert set(scenario.dataframes["sets"]["set_nodesdata"]) == {"DE_BW"}


class TestRunner:
    def test_normal_run(self, tmp_path):
        results = tmp_path / "results"
        results.mkdir()
        code = run_remix(datadir=MINIMAL_LP_PATH, resultdir=results, resultfile="remix")

        assert code == 0
        assert (tmp_path / "results" / "remix.gdx").exists

    def test_failing_run(self, tmp_path):
        results = tmp_path / "results"
        results.mkdir()
        code = run_remix(datadir=tmp_path, resultdir=results, resultfile="remix")

        assert code == 3


class TestInstance:
    def test_container_factory_sets_repr(self):
        MyClass = container_factory({"attribute": [1]}, "Sets")
        instance = MyClass()
        assert instance.__repr__() == "Sets with the following items: ['attribute']"
        assert instance.__str__() == "Sets with the following items: ['attribute']"

    def test_container_factory_error(self):
        with pytest.raises(TypeError) as e_info:
            MyClass = container_factory({"attribute": 1}, "MyClass")
            instance = MyClass()

    def test_instance_python_inference(self, example_dataframes):
        # remove sets
        example_dataframes_without_sets = deepcopy(example_dataframes)
        example_dataframes_without_sets["sets"] = {}

        model = Instance.from_dataframes(example_dataframes_without_sets)
        model.infer_set_data()
        model.update_dataframes()
        dfs = model.dataframes
        assert set(dfs["sets"]["set_nodesdata"]) == {"DE_BW", "DE_BA"}
        assert set(dfs["sets"]["set_accnodes"]) == {"Europe"}
        assert set(dfs["sets"]["set_years"]) == {"2050"}
        assert set(dfs["sets"]["set_converter_techs"]) == {"tech_a"}
        assert set(dfs["sets"]["set_indicators"]) == {"SystemCost"}

    def test_dataframes(self):
        instance = Instance.from_path(
            MINIMAL_LP_PATH, parent_files=None, version="latest", name="."
        )
        dataframes = instance.dataframes

        assert "sets" in dataframes
        assert "parameters" in dataframes
        assert "timeseries" in dataframes
        assert "maps" in dataframes

    def test_assign_attributes_descriptor(
        self, temporary_dat_files, example_dataframes
    ):
        inst = Instance.from_path(temporary_dat_files)
        inst._arrange_files()
        inst._split_parameters_sets()
        inst._inherit_dataframes()
        inst.assign_attributes()

        assert (
            inst.set.converter_techs
            == example_dataframes["sets"]["set_converter_techs"]
        )
        assert inst.set.nodesdata == example_dataframes["sets"]["set_nodesdata"]
        assert inst.set.years == example_dataframes["sets"]["set_years"]
        assert inst.set.indicators == example_dataframes["sets"]["set_indicators"]

        assert (
            inst.parameter.converter_capacityparam["unitsLowerLimit"].iloc[0]
            == example_dataframes["parameters"]["converter_capacityparam"][
                "unitsLowerLimit"
            ].iloc[0]
        )
        assert (
            inst.parameter.accounting_indicatorbounds["discount"].iloc[0]
            == example_dataframes["parameters"]["accounting_indicatorbounds"][
                "discount"
            ].iloc[0]
        )

    def test_assign_attributes_roundtrip(self, temporary_dat_files):
        inst = Instance.from_path(temporary_dat_files)
        inst.inherit()

        inst.set.nodesdata.append("FR")

        inst.update_dataframes()
        dfs = inst.dataframes

        assert set(dfs["sets"]["set_nodesdata"]) == set(["DE_BW", "FR"])

        dfs["sets"]["set_nodesdata"].append("ES")
        inst.assign_attributes(dfs)

        assert set(inst.set.nodesdata) == {"DE_BW", "FR", "ES"}

    def test_write_from_path_roundtrip(self, temporary_dat_files, tmp_path_factory):
        inst = Instance.from_path(temporary_dat_files)

        inst.set.nodesdata.append("FR")

        temporary_dat_path = tmp_path_factory.mktemp("tmp")
        inst.datadir = temporary_dat_path
        inst.write(fileformat="dat")

        inst_from_path = Instance.from_path(inst.datadir)

        assert set(inst_from_path.set.serialize()) == set(inst.set.serialize())
        assert set(inst_from_path.map.serialize()) == set(inst.map.serialize())
        assert set(inst_from_path.parameter.serialize()) == set(
            inst.parameter.serialize()
        )
        assert set(inst_from_path.profile.serialize()) == set(inst.profile.serialize())

    def test_writing_columns_with_only_default_values_but_explicitly_set_by_user(
        self, example_dataframes, tmp_path_factory
    ):
        inst = Instance.from_dataframes(example_dataframes)

        temporary_dat_path = tmp_path_factory.mktemp("tmp")
        inst.datadir = temporary_dat_path

        # add a new converter with unitsUpperLimit
        index = MultiIndex.from_product([["DE_SH"], ["2050"], ["tech_a", "tech_b"]])
        df = DataFrame({"unitsUpperLimit": 1}, index=index)
        inst.parameter.add(df, "converter_capacityparam")

        # the written data must have the "unitsLowerLimit" column, as the value 0.0 is explicitly set in the default dataframes
        inst.write(fileformat="dat")
        datcontents = read_dat(
            os.path.join(inst.datadir, "converter_capacityparam.dat")
        )
        assert "unitsLowerLimit" in datcontents.columns

        # now the column "unitsLowerLimit" should vanish: all data nan, all defaults zero
        inst.parameter.converter_capacityparam.loc[
            idx["DE_BW", "2050", "tech_a"], "unitsLowerLimit"
        ] = np.nan
        inst.write(fileformat="dat")
        datcontents = read_dat(
            os.path.join(inst.datadir, "converter_capacityparam.dat")
        )
        assert "unitsLowerLimit" not in datcontents.columns

        # now we want to get back that column, with setting one value that is equal to the default (0.0)
        inst.parameter.converter_capacityparam.loc[
            idx["DE_BW", "2050", "tech_a"], "unitsLowerLimit"
        ] = 0.0
        inst.write(fileformat="dat")
        datcontents = read_dat(
            os.path.join(inst.datadir, "converter_capacityparam.dat")
        )
        assert "unitsLowerLimit" in datcontents.columns

    def test_sorted_inference_of_sets(self, example_dataframes):
        # remove sets
        example_dataframes_without_sets = deepcopy(example_dataframes)
        example_dataframes_without_sets["sets"] = {}

        model = Instance.from_dataframes(example_dataframes_without_sets)

        # make two more additions to the instance
        converter_capacityparam = deepcopy(model.parameter.converter_capacityparam)
        converter_capacityparam.index = converter_capacityparam.index.set_levels(
            ["FR"], level=0
        )
        model.parameter.add(
            deepcopy(converter_capacityparam), "converter_capacityparam"
        )
        converter_capacityparam.index = converter_capacityparam.index.set_levels(
            [["DE"], ["tech_b"]], level=[0, 2]
        )
        model.parameter.add(
            deepcopy(converter_capacityparam), "converter_capacityparam"
        )

        # generate sets and write to dataframes
        model.infer_set_data()
        model.update_dataframes()

        # compare values of model.set and model.dataframe["sets"]
        for _key in model.set._attrs:
            if len(model.set.__dict__[f"_{_key}"]) > 0:
                assert (
                    model.set.__dict__[f"_{_key}"]
                    == model.dataframes["sets"][f"set_{_key}"]
                )

    def test_custom_defaults(self, example_dataframes, tmp_path_factory):
        model = Instance.from_dataframes(
            example_dataframes,
            custom_defaults={"converter_capacityparam": {"unitsBuild": 0.5}},
        )
        assert model.parameter.defaults["converter_capacityparam"]["unitsBuild"] == 0.5

        model.datadir = tmp_path_factory.mktemp("tmp")
        model.write(fileformat="dat")
        datcontents = read_dat(
            os.path.join(model.datadir, "converter_capacityparam.dat")
        )

        assert datcontents.loc[idx["DE_BW", "2050", "tech_a"], "unitsBuild"] == 0.5

    def test_invalid_custom_defaults(self):
        # test errors first level
        with pytest.raises(AttributeError):
            Instance.from_dataframes(
                example_dataframes,
                custom_defaults={"converter_capacitypasam": {"unitsBuild": 0.5}},
            )
        # and in second level, however it is weird that this is an error, any idea why @wetz_mn?
        with pytest.raises(AttributeError):
            Instance.from_dataframes(
                example_dataframes,
                custom_defaults={"converter_capacityparam": {"unitsBuilt": 0.5}},
            )

    def test_wrong_dtypes(self, example_dataframes, tmp_path_factory):
        model = Instance.from_dataframes(example_dataframes)

        model.parameter.converter_capacityparam.loc[:, "noExpansion"] = np.inf
        model.datadir = tmp_path_factory.mktemp("tmp")
        with pytest.raises(IntCastingNaNError):
            model.write(fileformat="dat")

    def test_sorted_sets(self):
        model = Instance()
        sets_to_add = ["C", "B", "10", "A"]
        model.set.add(sets_to_add, "nodesdata")

        assert model.set.nodesdata == sorted(sets_to_add)

    def test_addition_of_identical_sets(self):
        model = Instance()
        model.set.add(["C", "C"], "nodesdata")
        assert model.set.nodesdata == ["C"]

    def test_addition_of_mixed_lists_and_dataframes(self):
        model = Instance()
        model.set.add(["C", "C"], "nodesdata")
        assert model.set.nodesdata == ["C"]

    def test_addition_of_unordered_years_to_sets(self, tmp_path_factory):
        model = Instance()
        model.set.add(["2030", "2050", "2040"], "years")

        temporary_dat_path = tmp_path_factory.mktemp("tmp")

        # save model, reload it and check if order of years is sorted
        model.datadir = temporary_dat_path
        model.write()
        model_from_path = Instance.from_path(temporary_dat_path)

        assert model_from_path.set.years == model.set.years

    def test_adding_nonsupported_type(self):
        inst = Instance()
        with pytest.raises(TypeError):
            inst.parameter.add(
                "adding string datatype should not work", "converter_capacityparam"
            )
        # add a list to the parameters
        # this is currently possible, so it should throw errors, we can revive the test, once list support is dropped
        # with pytest.raises(TypeError):
        #     inst.parameter.add([DataFrame()], "converter_capacityparam")
        # add a dataframe to the sets
        with pytest.raises(TypeError):
            inst.set.add(DataFrame(), "accnodes")

    def test_adding_mismatched_columns(self, example_dataframes):
        inst = Instance()

        modified_df = deepcopy(example_dataframes)
        modified_df["parameters"]["converter_capacityparam"].rename(
            columns={"unitsUpperLimit": "unitsUpperLmit"}, inplace=True
        )
        with pytest.raises(ValueError):
            inst.parameter.add(
                modified_df["parameters"]["converter_capacityparam"],
                "converter_capacityparam",
            )

    @pytest.mark.skip(reason="Not implemented")
    def test_adding_mismatched_indices(self):
        raise NotImplementedError

    def test_adding_nonexistingkeyword(self):
        inst = Instance()
        placeholder = (
            "It should not matter what the data looks like, keyword is checked first!"
        )
        with pytest.raises(AttributeError):
            inst.set.add(placeholder, "name")

        with pytest.raises(AttributeError):
            inst.parameter.add(placeholder, "converter_techparm")

    def test_wrong_dataframes(self):
        defaults_shallow = default_dataframes(flat=True)
        wrong_dataframes = defaults_shallow.copy()
        wrong_dataframes["wrong_name"] = DataFrame()
        with pytest.raises(KeyError) as e_info:
            scenario = Instance.from_dataframes(wrong_dataframes)

    def test_scenario_wrong_keys(self, example_dataframes):
        wrong_dataframes = deepcopy(example_dataframes)
        wrong_dataframes["sets"]["set_wrong"] = ["WRONG"]
        with pytest.raises(KeyError) as e_info:
            scenario = Instance.from_dataframes(wrong_dataframes)

    def test_model_with_years_as_int(self, tmp_path):

        m = Instance.from_path(MINIMAL_LP_PATH)
        some_idx = m.parameter.converter_coefficient.index
        m.parameter.converter_coefficient.index = some_idx.set_levels(
            some_idx.levels[1].astype(int), level=1
        )
        years = m.set.years.copy()
        m.datadir = tmp_path
        m.write()
        assert m.set.years == years


class TestTransform:
    def test_transform_dat_to_no_index_csv(self, tmp_path):
        transform_dataset(MINIMAL_LP_PATH, tmp_path, outformat="csv", no_sets=True)
        assert _count_non_empty_files_in_dir(tmp_path) == _count_non_empty_files_in_dir(
            MINIMAL_LP_PATH
        )
        non_set_map = _get_files_by_prefix_exclusion(tmp_path, ["set_", "map_"])
        for _file in non_set_map:
            csv_in_name = ".csv" in _file.as_posix()
            assert csv_in_name, f"Error at {_file}"
            with open(_file, "r") as check_file:
                first_line = check_file.readline()
                no_names = bool(re.search(r"\,\,+", first_line))
                assert no_names, f"Error at {_file}"

    def test_transform_dat_to_index_csv(self, tmp_path):
        transform_dataset(MINIMAL_LP_PATH, tmp_path, outformat="csv")
        assert _count_non_empty_files_in_dir(tmp_path) == _count_non_empty_files_in_dir(
            MINIMAL_LP_PATH
        )
        non_set_map = _get_files_by_prefix_exclusion(tmp_path, ["set_", "map_"])
        for _file in non_set_map:
            csv_in_name = ".csv" in _file.as_posix()
            assert csv_in_name, f"Error at {_file}"
            with open(_file, "r") as check_file:
                first_line = check_file.readline()
                names = not bool(re.search(r"\,\,+", first_line))
                assert names, f"Error at {_file}"

    def test_transform_index_csv_to_dat(self, tmp_path):
        transform_dataset(MININAL_LP_NAMED_CSV_PATH, tmp_path, outformat="dat")
        assert _count_non_empty_files_in_dir(tmp_path) == 21
        non_set_map = _get_files_by_prefix_exclusion(tmp_path, ["set_", "map_"])
        for _file in non_set_map:
            dat_in_name = ".dat" in _file.as_posix()
            assert dat_in_name, f"Error at {_file}"
            with open(_file, "r") as check_file:
                first_line = check_file.readline()
                names = not bool(re.search(r"\s\s+\d+", first_line))
                assert names, f"Error at {_file}"

    def test_transform_no_index_csv_to_dat(self, tmp_path):
        transform_dataset(MINIMAL_LP_CSV_PATH, tmp_path, outformat="dat")
        assert _count_non_empty_files_in_dir(tmp_path) == _count_non_empty_files_in_dir(
            MINIMAL_LP_CSV_PATH
        )
        non_set_map = _get_files_by_prefix_exclusion(tmp_path, ["set_", "map_"])
        for _file in non_set_map:
            dat_in_name = ".dat" in _file.as_posix()
            assert dat_in_name, f"Error at {_file}"
            with open(_file, "r") as check_file:
                first_line = check_file.readline()
                names = not bool(re.search(r"\s\s+\d+", first_line))
                assert names, f"Error at {_file}"

    def test_add_index_csv(self, tmp_path):
        transform_dataset(MINIMAL_LP_CSV_PATH, tmp_path, outformat="csv")
        assert _count_non_empty_files_in_dir(tmp_path) == _count_non_empty_files_in_dir(
            MINIMAL_LP_CSV_PATH
        )
        non_set_map = _get_files_by_prefix_exclusion(tmp_path, ["set_", "map_"])
        for _file in non_set_map:
            csv_in_name = ".csv" in _file.as_posix()
            assert csv_in_name, f"Error at {_file}"
            with open(_file, "r") as check_file:
                first_line = check_file.readline()
                names = not bool(re.search(r"\,\,+", first_line))
                assert names, f"Error at {_file}"

    def test_remove_index_csv(self, tmp_path):
        transform_dataset(
            MININAL_LP_NAMED_CSV_PATH, tmp_path, outformat="csv", no_sets=True
        )
        assert _count_non_empty_files_in_dir(tmp_path) == 21
        non_set_map = _get_files_by_prefix_exclusion(tmp_path, ["set_", "map_"])
        for _file in non_set_map:
            csv_in_name = ".csv" in _file.as_posix()
            assert csv_in_name, f"Error at {_file}"
            with open(_file, "r") as check_file:
                first_line = check_file.readline()
                no_names = bool(re.search(r"\,\,+", first_line))
                assert no_names, f"Error at {_file}"

    def test_transform_from_tall_normal_maps_to_dat(
        self, csv_tall_normal_maps_path, tmp_path
    ):
        transform_dataset(csv_tall_normal_maps_path, tmp_path, outformat="dat")
        assert _count_non_empty_files_in_dir(tmp_path) == _count_non_empty_files_in_dir(
            csv_tall_normal_maps_path
        )
        non_set_map = _get_files_by_prefix_exclusion(tmp_path, ["set_", "map_"])
        for _file in non_set_map:
            dat_in_name = ".dat" in _file.as_posix()
            assert dat_in_name, f"Error at {_file}"
            with open(_file, "r") as check_file:
                first_line = check_file.readline()
                names = not bool(re.search(r"\s\s+\d+", first_line))
                assert names, f"Error at {_file}"

    def test_transform_dat_to_csv_scenarios(self, tmp_path):
        transform_dataset(MINIMAL_LP_SUBSCENARIOS_PATH, tmp_path, outformat="csv")
        assert _count_non_empty_files_in_dir(tmp_path) == _count_non_empty_files_in_dir(
            MINIMAL_LP_SUBSCENARIOS_PATH
        )
        non_set_map = _get_files_by_prefix_exclusion(tmp_path, ["set_", "map_"])
        scenario_dirs = [fil for fil in Path(tmp_path).iterdir() if fil.is_dir()]
        for _file in non_set_map:
            csv_in_name = ".csv" in _file.as_posix()
            assert csv_in_name, f"Error at {_file}"
            with open(_file, "r") as check_file:
                first_line = check_file.readline()
                no_names = bool(re.search(r"\,\,+", first_line))
                assert not no_names, f"Error at {_file}"

        scen = scenario_dirs[0]
        assert _count_non_empty_files_in_dir(scen) == 6
        non_set_map_scen = _get_files_by_prefix_exclusion(scen, ["set_", "map_"])

        for _file in non_set_map_scen:
            csv_in_name = ".csv" in _file.as_posix()
            assert csv_in_name, f"Error at {_file}"
            with open(_file, "r") as check_file:
                first_line = check_file.readline()
                no_names = bool(re.search(r"\,\,+", first_line))
                assert not no_names, f"Error at {_file}"

    def test_transform_frictionless(self, tmp_path):
        transform_dataset(MINIMAL_LP_PATH, tmp_path, frictionless=True)
        assert _count_non_empty_files_in_dir(tmp_path) == 42
        non_set = _get_files_by_prefix_exclusion(tmp_path, ["set_"])
        for _file in non_set:
            csv_in_name = ".csv" in _file.as_posix()
            assert csv_in_name, f"Error at {_file}"
            with open(_file, "r") as check_file:
                first_line = check_file.readline()
                names = not bool(re.search(r"\,\,+", first_line))
                assert names, f"Error at {_file}"

        with open(
            tmp_path.joinpath("converter_activityprofile.csv"), "r"
        ) as check_file:
            first_line = check_file.readline()
            assert set(first_line.strip().split(",")) == {
                "nodesData",
                "timeData",
                "profileTypes",
                "converter_techs",
                "Value",
                "years",
            }

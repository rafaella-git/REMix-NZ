# SPDX-FileCopyrightText: Copyright (c) 2023 German Aerospace Center (DLR)
# SPDX-License-Identifier: BSD-3-Clause

import json
import os
import time

import numpy as np
import pandas as pd
import pytest
import yaml

from remix.framework import __testingpath__
from remix.framework.api.run import run_remix
from remix.framework.tools.gdx import gt
from remix.framework.tools.utilities import merge_dicts


def pytest_generate_tests(metafunc):
    if 'testcase_path' in metafunc.fixturenames:
        tc_list = [os.path.join(root, cur_dir)
                   for root, dirs, files in os.walk(__testingpath__)
                   for cur_dir in dirs
                   if (
                       os.path.isfile(os.path.join(root, cur_dir, 'valid.json'))
                       and os.path.isfile(os.path.join(root, cur_dir, 'config.yaml')))
                   ]
        tc_filter = metafunc.config.getoption("--specific")
        timelimit = metafunc.config.getoption("--timelimit")
        mode = metafunc.config.getoption("--mode")
        keep = metafunc.config.getoption("--keep")
        if mode == '.' or mode == 'solve':
            mode = 'solve'
            f = open('./metrics.txt', 'w')
            f.write('# TYPE runtime gauge\n')
            f.close()
        tc_list_selected = []
        keep = int(int(keep) == 1) if keep.isnumeric() else 0
        if len(tc_list_selected) > 1:
            keep = 0
        metafunc.parametrize('testcase_path', tc_list)
        metafunc.parametrize('mode', [mode])
        metafunc.parametrize('keep', [keep])
        metafunc.parametrize('tc_filter', [tc_filter])
        metafunc.parametrize('timelimit', [timelimit])


def gdx_to_json(gdxfile, check_symbols, round=1):
    m = gt.Container(gdxfile)

    valid_dict = {}

    for sym, obj in m.data.items():
        if sym not in check_symbols:
            continue

        df = m.data[sym].records

        if df is None:
            continue

        if isinstance(obj, gt.Set):
            df = df.set_index(df.columns.to_list()[:-1])
            if isinstance(df, pd.Series):
                df = df.to_frame()
            df.rename(columns={"element_text": "value"}, inplace=True)
            df["value"] = 1
        if isinstance(obj, gt.Parameter):
            df = df.set_index(df.columns.to_list()[:-1])
        if isinstance(obj, gt.Variable):
            df = df.set_index(df.columns.to_list()[:-5]).stack()
            df.index.names = df.index.names[:-1] + ["lmlus"]
            df = pd.DataFrame(df).rename(columns={0: "value"})

        has_timeseries = False
        idx_time = [c for c in df.index.names if "timeModel" in c]
        if len(idx_time) == 1 and len(df.index.names) > 1:
            has_timeseries = True
        if len(idx_time) > 1:
            raise("Cannot handle timeseries with multiple time dimensions")

        if has_timeseries:
            df = df.unstack(idx_time).fillna(0)

        # Remove default values
        if isinstance(obj, gt.Parameter):
            match_idx = pd.MultiIndex.from_tuples({idx for idx in df.index})
            if len(match_idx) > 0:
                df = df.drop(match_idx[(df.loc[match_idx] == 0).all(axis=1)])

        if isinstance(obj, gt.Variable):
            match_idx = pd.MultiIndex.from_tuples({idx for idx in df.index if idx[-1] in ['upper']})
            if len(match_idx) > 0:
                df = df.drop(match_idx[(df.loc[pd.MultiIndex.from_tuples(match_idx)] == np.inf).all(axis=1)])
            match_idx = pd.MultiIndex.from_tuples({idx for idx in df.index if idx[-1] in ['lower']})
            if len(match_idx) > 0:
                df = df.drop(match_idx[(df.loc[pd.MultiIndex.from_tuples(match_idx)] == -np.inf).all(axis=1)])
            match_idx = pd.MultiIndex.from_tuples({idx for idx in df.index if idx[-1] in ['level', 'marginal', 'lower', 'upper']})
            if len(match_idx) > 0:
                df = df.drop(match_idx[(df.loc[match_idx] == 0).all(axis=1)])
            match_idx = pd.MultiIndex.from_tuples({idx for idx in df.index if idx[-1] in ['scale']})
            if len(match_idx) > 0:
                df = df.drop(match_idx[(df.loc[pd.MultiIndex.from_tuples(match_idx)] == 1).all(axis=1)])

        df_new = pd.DataFrame(index=df.index)
        if round == 1:
            gfunc = np.vectorize(lambda x: float('%.6g' % x))
            df_new["value"] = list(map(gfunc, df.values))
        else:
            df_new["value"] = df["value"]

        dct = df_new.to_dict()["value"]

        if has_timeseries:
            valid_dict[sym] = {".".join(k): list(v) for k, v in dct.items()}
        elif not isinstance(obj, gt.Set):
            valid_dict[sym] = {".".join(k): float(v) for k, v in dct.items()}
        else:
            valid_dict[sym] = {k: float(v) for k, v in dct.items()}

    with open(gdxfile.replace(".gdx", ".json"), "w") as outfile:
        json.dump(valid_dict, outfile, indent=2)


def diff_json(vfile, cfile, eps=1e-8, relEps=1e-5):

    gfunc = lambda x: float('%.8g' % x)

    with open(vfile) as vf:
        valid = json.load(vf)

    with open(cfile) as cf:
        check = json.load(cf)

    diff = {}
    errors = []

    for sym, dct in valid.items():
        if sym not in check:
            errors.append("Symbol {} missing in checkfile".format(sym))
            diff[sym] = "missing"
            continue
        else:
            diff[sym] = {}

        for k, v in dct.items():
            if k not in check[sym]:
                if isinstance(v, list):
                    errors.append("Timeseries {}({}) missing in checkfile".format(sym, k))
                    diff[sym][k] = "missing"
                    continue
                else:
                    errors.append("Element {}({}) missing in checkfile".format(sym, k))
                    diff[sym][k] = "missing"
                    continue

            if isinstance(v, list):
                if len([(i, j) for i, j in zip(v, check[sym][k]) if abs(i - j) > eps]) > 0:
                    errors.append("Timeseries {}({}) differs from valid (absolute)\n".format(sym, k))
                    diff[sym][k] = [gfunc(j - i) for i, j in zip(v, check[sym][k])]
                if len([(i, j) for i, j in zip(v, check[sym][k]) if i != 0 if abs(i - j) / abs(i) > relEps]) > 0:
                    errors.append("Timeseries {}({}) differs from valid (relative)\n".format(sym, k))
                    diff[sym][k] = [gfunc(j - i) for i, j in zip(v, check[sym][k])]
            else:
                if abs(v - check[sym][k]) > eps:
                    errors.append("Element {}({}) differs from valid (absolute)\n".format(sym, k))
                    diff[sym][k] = gfunc(check[sym][k] - v)
                if v != 0:
                    if abs(v - check[sym][k]) / abs(v) > relEps and abs(v - check[sym][k]) > eps:
                        errors.append("Element {}({}) differs from valid (relative)\n".format(sym, k))
                        diff[sym][k] = gfunc(check[sym][k] - v)

    diff = {i: j for i, j in diff.items() if len(j) > 0}
    if len(diff) > 0:
        with open(vfile.replace("valid", "diff"), "w") as outfile:
            json.dump(diff, outfile, indent=2)

    return errors


class ModelRun:
    def __init__(self, tc_path: str, mode: str, keep: int, tc_filter: str, timelimit: int):
        self.tc_path = tc_path
        self.mode = mode
        self.keep = keep
        self.tc_filter = tc_filter
        self.timelimit = timelimit
        self.config_file = None
        self.trusted_result_file = None
        self.gdx_json_parameters = {}
        self.diff_parameters = {}
        self.check_symbols = []
        self.error_messages = {}
        self.symbol_error = ''
        self.data_error = ''
        self.metrics = {}
        self.remix_parameters = {}
        self.search_cfg()
        self.read_cfg()
        self.remix_parameters["datadir"] = os.path.join(tc_path, self.remix_parameters["datadir"])
        self.remix_parameters["resultdir"] = os.path.join(tc_path, self.remix_parameters["resultdir"])
        self.remix_parameters["logfile"] = os.path.join(tc_path, self.remix_parameters["logfile"])
        self.remix_parameters["output"] = os.path.join(tc_path, self.remix_parameters["output"])
        self.remix_parameters["optdir"] = os.path.join(tc_path, self.remix_parameters["optdir"])
        if "instancedir" in self.remix_parameters:
            self.remix_parameters["instancedir"] = os.path.join(tc_path, self.remix_parameters["instancedir"])

    def exec(self):
        self.check_files()
        if self.tc_filter not in self.tc_path:
            pytest.skip("TC filtered by name.")
        if 'runtime' in self.metrics:
            expected_runtime = ''.join([c for c in self.metrics['runtime'] if c in '1234567890.'])
            expected_runtime = float(expected_runtime)
            if expected_runtime > float(self.timelimit):
                pytest.skip("TC filtered due to expected runtime.")
        if self.mode == 'solve':
            self.run_model()
            self.compare_results()
        if self.mode == 'metrics':
            self.compare_metrics()

    def search_cfg(self):
        for file in os.listdir(self.tc_path):
            if file == 'config.yaml':
                self.config_file = 'config.yaml'
            if file == 'valid.json':
                self.trusted_result_file = 'valid.json'

    def check_files(self):
        dat_ctr = 0
        for dat_file in os.listdir(self.remix_parameters["datadir"]):
            if dat_file.endswith('.dat') or dat_file.endswith('.csv'):
                dat_ctr += 1
        if dat_ctr < 5:
            raise RuntimeError("Suspiciously few .dat/.csv files!")

    def read_cfg(self):
        cfg_files = [__testingpath__.as_posix() + '/testing_global.yaml']
        if self.config_file is not None:
            cfg_files.append(self.tc_path + '/' + self.config_file)
        for cfg_file in cfg_files:
            with open(cfg_file) as cur_file:
                config = yaml.safe_load(cur_file)

                for param in [k for k in config.keys() if k not in ["check_symbols", "skip_symbols"]]:
                    data = getattr(self, param)
                    setattr(self, param, merge_dicts(data, config.get(param, {})))
                    if param in ["diff_parameters", "gdx_json_parameters"]:
                        data = getattr(self, param)
                        for k, v in data.items():
                            try:
                                data[k] = float(v)
                            except ValueError:
                                pass

                # read runtime parameters
                if 'check_symbols' in config:
                    self.check_symbols += config['check_symbols']
                    self.check_symbols = list(set(self.check_symbols))
                if 'skip_symbols' in config:
                    self.check_symbols = list(set(self.check_symbols) - set(config['skip_symbols']))

    def run_model(self):
        start_time = time.time()
        run_gams = run_remix(lo=4, keep=self.keep, **self.remix_parameters)
        run_time = time.time() - start_time
        f = open('./metrics.txt', 'a')
        f.write('runtime{{tc="{}"}} {:1.3f}\n'.format(self.tc_path, run_time))
        f.close()
        assert run_gams == 0, "\n## Failed to run model\n{}".format(self.error_messages)

    def compare_results(self):
        if self.trusted_result_file is not None:
            if 'round' in self.gdx_json_parameters:
                gdx_to_json(f"{self.tc_path}/check.gdx", self.check_symbols, round=self.gdx_json_parameters['round'])
            else:
                gdx_to_json(f"{self.tc_path}/check.gdx", self.check_symbols)
            error_list = diff_json(f"{self.tc_path}/valid.json", f"{self.tc_path}/check.json", **self.diff_parameters)

            error_report = ''
            for err in error_list:
                error_report += err

            assert len(error_list) == 0, '\n## Data comparison failed' + error_report

    def compare_metrics(self):
        file = open('./metrics.txt', 'r')
        all_metrics = file.readlines()
        file.close()
        assert len(all_metrics) > 0
        for metric, metric_requirement in self.metrics.items():
            metric_found = False
            for cur_metric in all_metrics:
                if cur_metric.find(metric) == 0 and cur_metric.find('{{tc="{}"}}'.format(self.tc_path)) == len(metric):
                    metric_found = True
                    space_pos = cur_metric.find(' ')
                    recorded_metric_val = cur_metric[space_pos + 1:-1]
                    comp_op = 'eq'
                    try:
                        recorded_metric_val = float(recorded_metric_val)
                    except ValueError:
                        raise ValueError(f"Metric value {recorded_metric_val} in metrics.txt cannot be converted to float.")

                    for i in range(1, 3):
                        try:
                            comp_op = metric_requirement[:i]
                            metric_requirement = float(metric_requirement[i:])
                            break
                        except ValueError:
                            if i == 2:
                                raise Exception(f"Could not translate {metric_requirement} into valid metric criterion")

                    msg = f"Metric requirement for test is {metric_requirement}, recorded value is {recorded_metric_val}."
                    if comp_op in ['eq', '=', '==']:
                        expr = recorded_metric_val == metric_requirement
                        msg += "The values must be equal."
                    elif comp_op in ['<=', 'le']:
                        expr = recorded_metric_val <= metric_requirement
                        msg += "The result must be less than or equal to the requirement."
                    elif comp_op in ['<', 'lt', 'l']:
                        expr = recorded_metric_val < metric_requirement
                        msg += "The result must be less than the requirement."
                    elif comp_op in ['>=', 'ge']:
                        expr = recorded_metric_val >= metric_requirement
                        msg += "The result must be greater than or equal to the requirement."
                    elif comp_op in ['>', 'gt', 'g']:
                        expr = recorded_metric_val > metric_requirement
                        msg += "The result must be greater than the requirement."
                    assert expr, msg

            if not metric_found:
                raise Exception(f'Could not find {metric}{{tc="{self.tc_path}"}} in metrics.txt. Therefore, the conditions imposed by the testcase configuration could not be checked.')


def testcase(testcase_path, mode, keep, tc_filter, timelimit):
    try:
        timelimit = float(timelimit)
    except ValueError:
        timelimit = 7200
    curTC = ModelRun(testcase_path, mode, keep, tc_filter, timelimit)
    curTC.exec()

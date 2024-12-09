from unittest.mock import patch

from typer.testing import CliRunner

from remix.framework.cli.main import program

runner = CliRunner()


class TestMain:
    def test_main_help(self):
        result = runner.invoke(program, ["--help"])
        assert result.exit_code == 0
        for option in ["--version", "run", "test"]:
            assert option in result.stdout

    def test_main_version(self):
        result = runner.invoke(program, ["--version"])
        assert result.exit_code == 0
        assert "REMix version:" in result.stdout


class TestRun:
    def test_run_help(self):
        result = runner.invoke(program, ["run", "--help"])
        assert result.exit_code == 0

        for option in ["--datadir", "--scendir", "--resultdir", "--resultfile",
                       "--lo", "--names", "--postcalc", "--roundts", "--keep"]:
            assert option in result.stdout

    def test_run_min(self, tmp_path_factory):
        tmp = tmp_path_factory.mktemp("tmp")
        tmp_path = tmp.as_posix()
        # run_remix_mock = mocker.patch("remix.framework.cli.run.run_remix")
        with patch("remix.framework.cli.run.run_remix") as run_remix_mock:
            result = runner.invoke(program, ["run", "--datadir=testing/instances/minimal_lp/data",
                                             f"--resultdir={tmp_path}"])
            assert result.exit_code == 0
            assert run_remix_mock.call_args[1]["datadir"] == "testing/instances/minimal_lp/data"
            assert run_remix_mock.call_args[1]["resultdir"] == f"{tmp_path}"

    def test_run_unkown_argument(self):
        tmp_path = "result_directory_mock"
        # run_remix_mock = mocker.patch("remix.framework.cli.run.run_remix")
        with patch("remix.framework.cli.run.run_remix") as run_remix_mock:
            result = runner.invoke(program, ["run", "--datadir=testing/instances/minimal_lp/data",
                                             f"--resultdir={tmp_path}", "--some_unkown_argument=some_value"])
            assert result.exit_code == 0
            assert run_remix_mock.call_args[1]["datadir"] == "testing/instances/minimal_lp/data"
            assert run_remix_mock.call_args[1]["resultdir"] == f"{tmp_path}"
            assert run_remix_mock.call_args[1]["some_unkown_argument"] == "some_value"


class TestTest:
    def test_test_help(self):
        result = runner.invoke(program, ["test", "--help"])
        assert result.exit_code == 0

        for option in ["--specific", "--mode", "--junitxml", "--keep", "--timelimit"]:
            assert option in result.stdout

    def test_test_run(self):
        # Not possible to mock pytest if it is not imported at the top, this gets too meta
        # we run the normal call once and hope it does not take too long
        # remix_test_mock =  mocker.patch("remix.framework.cli.test.pytest")
        result = runner.invoke(program, ["test", "--specific=exact_data_and_inheritance"])
        assert result.exit_code == 0


class TestTransform:
    def test_transform_help(self):
        result = runner.invoke(program, ["transform", "--help"])
        assert result.exit_code == 0

        for option in ["--datadir", "--outputdir", "--outformat", "--mapformat", "--profileformat", "--supportsets", "--frictionless"]:
            assert option in result.stdout

    def test_transform_frictionless(self):
        tmp_path = "result_directory_mock"
        with patch("remix.framework.cli.transform.transform_dataset") as transform_dataset_mock:
            result = runner.invoke(program, ["transform", "--datadir=testing/instances/minimal_lp/data",
                                             f"--outputdir={tmp_path}", "--frictionless=1"])
            assert result.exit_code == 0
            assert transform_dataset_mock.call_args[1]["datadir"] == "testing/instances/minimal_lp/data"
            assert transform_dataset_mock.call_args[1]["outputdir"] == f"{tmp_path}"
            assert transform_dataset_mock.call_args[1]["frictionless"]

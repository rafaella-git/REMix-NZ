# run_scenarios.py
# ------------------------------------------------------------------------------------
# Run REMix optimisation for a base case and (optionally) its scenarios
#
# This script:
#   1. Runs a full optimisation for the base case
#   2. Automatically detects any subfolders inside /data/
#      (each subfolder is treated as a scenario)
#   3. Allows the user to select which scenarios to run, or run all
#   4. Saves all results in the /result/ folder, organized by scenario
#
# Notes:
# - Each scenario folder inside /data/ should contain only the modified files
#   (for example, demand, inflow, or cost variations).
# - The model automatically uses the base files from /data/
#   and overrides them with the scenario files.
# - If the model is infeasible or a run crashes
#   If you build a model on your own, it can always happen that it gets infeasible 
#   or the execution stops out of some other reason. In that case, you can refer 
#   to the `run_remix.lst` file and look for the error marker `****`.
#   You can open that file in an editor and look for the error message.
#   Alternatively, you can also use the commandline tool grep to search for the 
#   pattern "**** Exec Error"  in the file to see what is wrong.
# ------------------------------------------------------------------------------------

import os
import time
from pathlib import Path
import pandas as pd
from remix.framework.api.instance import Instance

# user settings
group_name = "hadi"            # main folder under /project/
base_case = "pypsa-cascade"    # name of the model case folder
run_all_scenarios = True       # if False, the script will ask which to run

# set paths
base_dir = Path(f"../project/{group_name}/{base_case}")
data_dir = base_dir / "data"
result_dir = base_dir / "result"
result_dir.mkdir(parents=True, exist_ok=True)

# shared run options
run_options = dict(
    solver="cplex",      # choose gurobi or cplex for debugging
    datacheck=1,         # `datacheck` : will make CPLEX give out warnings about disproportionate values
    threads=8,           # number of CPU threads to use
    keep=1,              # keep scratch folder “/225a” with exported .csvs for debugging
    lo=4,                # write a .log file for solver output
    names=1,             # include variable names in .lst file
    roundts=1,           # round timeseries to avoid numerical errors
    timeres=1,           # hourly time resolution (1 hour); use 24 for daily aggregation, etc.
    postcalc=1,          # run post-calculation
    pathopt="myopic",    # optimisation mode
)

group_name = "hadi"
base_case  = "pypsa-cascade"       # base case folder under ../project/{group_name}/{base_case}/data
RUN_DISPATCH = True  # set True if you also want the dispatch step

# Folder setup
base_dir = Path(f"../project/{group_name}/{base_case}")
base_data   = base_dir / "data"
base_result = base_dir / "result"

dry_scenario = "dry-year"          # folder under data/
dry_data     = base_data / dry_scenario
dry_result   = base_result / dry_scenario


base_result.mkdir(parents=True, exist_ok=True)
dry_result.mkdir(parents=True, exist_ok=True)


# Notes 
# - Put only modified files under: ../project/{group_name}/{case_name}/data/<scenario>/ 
# - REMix loads base files from data/ and overrides with any file found in data/<scenario>/.

# Explanation for command line arguments to GAMS function call
# `lo=3` : log option of GAMS; ensures that the output from GAMS will be visible in the terminal (`lo=4`
#          will write out a log file)
# `modeldir` : path to model directory (can be set relative or absolute)
# `datadir` : instruction to look for `*.dat` files in a specific directory
# `scendir` : instruction to look for `*.dat` files in a specific directory; REMix default: datadir. Path to the scenario to calculate, relative to the --datadir path
# `solver` : Solver to use for optimising the model; supported options are "copt", "cplex", "gurobi", "highs", "mosek", "xpress", "convert", and "scip"; default: "cplex"

# `solvermethod` : Method of chosen solver to use for optimising the model; supported options are "automatic" (0), "barrier" (1), "simplex" (2), "primal" (3), "dual" (4), "concurrent" (5), "network" (6), "sifting" (7), pdlp (8), and mip (9); there is no solver that would support all these methods, if an unavailable solver method is specified, a hint will be given in the solve log file and the default will be used instead; default: "barrier"/"ipm"
# `resultdir` : instruction to store the resulting `*.gdx` file in a specific directory
# `resultfile` : allows to give a custom name to the output `*.gdx` file; REMix default: "remix"
# `postcalc` : instruction whether to do a post-calculation; REMix default: 1
# `names` : instruction to give out actual variable names in case of error messages; REMix default: 0
# `roundts` : instruction to automatically round after-comma digits in large time series files where necessary to successfully read them in; REMix default: 0
# `keep` : keeps temporary folder with `*.csv` files built during the model run, might be interesting for debugging purposes; REMix default: 0

# solver arguments
# `crossover` : Run the optimisation with crossover (1) or without crossover (0, default)
# `datacheck` : will make CPLEX give out warnings about disproportionate values (which lead to numerical difficulties = non-optimal solutions) at the beginning of the optimisation, default 0

# `barcolnz` : Number of non-zero entries above which CPLEX solver will treat columns as dense, default 0 (which means parameter is determined dynamically)


# function to run a SINGLE optimisation
def run_optimisation(data_folder, result_folder, case_name):
    """
    Runs a full REMix optimisation for a given case or scenario.

    Parameters:
    - data_folder : The folder containing the input data (usually /data or /data/<scenario>/)
    - result_folder : The folder where results will be saved
    - case_name : Name used for result files (e.g. pypsa-cascade_opt)
    """
    print(f"\nRunning optimisation for case: {case_name}")
    print(f"Input data folder: {data_folder}")
    print(f"Output result folder: {result_folder}")

    # Decide whether this is the base case or a scenario
    if data_folder == data_dir:
        # Base case: only use base data
        scenario_name = None
        print("Scenario: none (base case)")
    else:
        # Scenario case: use overrides from /data/<scenario_name>/
        scenario_name = data_folder.name
        print(f"Scenario detected: {scenario_name}")

    # Initialise the model instance
    model = Instance(index_names=False,
                     datadir=data_dir,
                     scendir=scenario_name)

    # Copy solver and model options
    options = run_options.copy()
    options.update(dict(resultdir=result_folder, resultfile=case_name))

    # Run optimisation
    start = time.time()
    result_code = model.run(**options)
    elapsed = time.time() - start

    # Check result
    if result_code != 0:
        print(f"REMix returned error code {result_code}. Stopping execution.")
        raise SystemExit(result_code)

    print(f"Optimization for {case_name} completed successfully.")
    print(f"Time taken: {elapsed/60:.1f} minutes\n")

    return result_folder / f"{case_name}.gdx"


# main run sequence

# record the total start time
t_start = time.time()
print("Starting REMix model runs")

# step 1: run the base case
base_result_file = run_optimisation(data_dir, result_dir, f"{base_case}_opt")

# step 2: check for scenario subfolders inside /data/
scenario_dirs = [f for f in data_dir.iterdir() if f.is_dir()]
if not scenario_dirs:
    print("No scenario folders found in the data directory.")
else:
    print("\nScenarios detected:")
    for i, s in enumerate(scenario_dirs, 1):
        print(f"  {i}. {s.name}")

    # let user decide what to run if not running all
    if not run_all_scenarios:
        selected = input("\nEnter scenario names to run (comma-separated), or type 'all': ").strip()
        if selected.lower() == "all":
            scenarios_to_run = scenario_dirs
        else:
            names = [x.strip() for x in selected.split(",")]
            scenarios_to_run = [f for f in scenario_dirs if f.name in names]
    else:
        scenarios_to_run = scenario_dirs

    # run each selected scenario
    if scenarios_to_run:
        for scen_dir in scenarios_to_run:
            scen_name = scen_dir.name
            scen_result_dir = result_dir / scen_name
            scen_result_dir.mkdir(parents=True, exist_ok=True)
            run_optimisation(scen_dir, scen_result_dir, f"{base_case}_{scen_name}_opt")
    else:
        print("No scenarios selected for execution.")

# print summary
t_end = time.time()
total_time = time.strftime("%Hh %Mm %Ss", time.gmtime(t_end - t_start))
print("All optimisation runs completed successfully.")
print(f"Total runtime: {total_time}")
print(f"Results saved in: {result_dir}")

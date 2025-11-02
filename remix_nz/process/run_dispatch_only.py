# run_dispatch_only.py
# -----------------------------------------------------------------------------
# Run REMix in dispatch-only mode using fixed capacities from a previous run
#
# This script:
#   1. Uses an existing solved model (GDX file) as a capacity reference
#   2. Runs a new dispatch-only simulation with fixed capacities
#   3. Allows the user to test different inflow or demand scenarios
#
# Example use case:
#   - You solved a base case (normal inflow) and want to test dispatch
#     under a "dry-year" inflow scenario using the same installed capacities.
# -----------------------------------------------------------------------------

import os
import time
from pathlib import Path
import pandas as pd
from remix.framework.api.instance import Instance

# -----------------------------------------------------------------------------
# basic user settings
group_name = "hadi"            # main folder under /project/
base_case = "pypsa-cascade"    # base case name (same as optimization run)
solver_choice = "cplex"        # choose solver: cplex, gurobi, etc.

# -----------------------------------------------------------------------------
# folder setup
base_dir = Path(f"../project/{group_name}/{base_case}")
data_dir = base_dir / "data"
result_dir = base_dir / "result"
dispatch_dir = base_dir / "dispatch"

# create dispatch output folder if it does not exist
dispatch_dir.mkdir(parents=True, exist_ok=True)

# -----------------------------------------------------------------------------
# run options
run_options = dict(
    solver=solver_choice,
    datacheck=1,
    threads=8,
    keep=1,
    lo=4,
    names=1,
    roundts=1,
    timeres=1,
    postcalc=1,
    pathopt="myopic", #foresight (perfect foresight), 
		#               miopic (limited foresight - like 2020 to 2030) 
		#               or target (single year -it takes the last year from yearssel-)
)

# -----------------------------------------------------------------------------
# helper functions

def list_gdx_files(folder):
    """Return a list of all GDX files in a folder and its subfolders."""
    return list(folder.rglob("*.gdx"))

def choose_from_list(options, prompt):
    """Ask the user to select an option from a list."""
    print(f"\n{prompt}")
    for i, opt in enumerate(options, 1):
        print(f"  {i}. {opt}")
    while True:
        choice = input("Enter number: ").strip()
        if choice.isdigit() and 1 <= int(choice) <= len(options):
            return options[int(choice) - 1]
        print("Invalid choice. Please enter a valid number.")

# -----------------------------------------------------------------------------
# main function to run dispatch

def run_dispatch(data_folder, fixed_gdx_path, output_folder, case_name):
    """
    Run a dispatch-only REMix simulation using fixed capacities.

    Parameters
    ----------
    data_folder : Path
        Data folder used for this run (/data or /data/<scenario_name>/)
    fixed_gdx_path : Path
        Path to the GDX file containing fixed capacities
    output_folder : Path
        Folder to store dispatch results
    case_name : str
        Name of the result file (e.g. pypsa-cascade_dry-year_dispatch)
    """

    print(f"\nStarting dispatch-only run: {case_name}")
    print(f"Using data from: {data_folder}")
    print(f"Using fixed capacities from: {fixed_gdx_path}")
    print(f"Results will be saved in: {output_folder}")

    # Determine scenario folder name for REMix (None for base, or subfolder name)
    if data_folder == data_dir:
        scenario_name = None
        print("Scenario: none (base data)")
    else:
        scenario_name = data_folder.name
        print(f"Scenario detected: {scenario_name}")

    # Initialize the model instance
    model = Instance(index_names=False, datadir=data_dir, scendir=scenario_name)

    # Freeze expansion (optional redundancy)
    try:
        nodes = list(model.set.nodesdata)
        years = list(model.set.yearssel)
        techs = list(model.set.converters)
        lock = pd.DataFrame(
            index=pd.MultiIndex.from_product([nodes, years, techs],
                                             names=["region", "years", "technology"]),
            columns=["noExpansion"]
        )
        lock["noExpansion"] = 1
        model.parameter.add(lock, "converter_capacityparam")
        print("All capacity expansion disabled (noExpansion = 1).")
    except Exception as e:
        print(f"Warning: could not apply noExpansion programmatically: {e}")

    # Set options for dispatch run
    options = run_options.copy()
    options.update(dict(
        resultdir=output_folder,
        resultfile=case_name,
        fixedcapsfromgdx=str(fixed_gdx_path.as_posix())
    ))

    # Run the dispatch
    start = time.time()
    rc = model.run(**options)
    elapsed = time.time() - start

    if rc != 0:
        print(f"REMix returned error code {rc}. Dispatch run failed.")
        raise SystemExit(rc)

    print(f"Dispatch run '{case_name}' completed successfully.")
    print(f"Time taken: {elapsed/60:.1f} minutes\n")

# -----------------------------------------------------------------------------
# main program 

if __name__ == "__main__":

    print("------------------------------------------------------------------------")
    print("Starting REMix dispatch-only run")
    print("------------------------------------------------------------------------")

    # Step 1: list available GDX files
    all_gdx = list_gdx_files(result_dir)
    if not all_gdx:
        print("No GDX files found in the results directory.")
        raise SystemExit(1)

    # Let the user choose a GDX file (e.g. base or scenario result)
    chosen_gdx = choose_from_list(all_gdx, "Select☺ the GDX file to use for fixed capacities:")

    # Step 2: list available data folders (base and scenarios)
    data_folders = [data_dir] + [f for f in data_dir.iterdir() if f.is_dir()]
    chosen_data = choose_from_list(data_folders, "Select the data configuration for this dispatch run:")

    # Step 3: set output file name
    if chosen_data == data_dir:
        scenario_label = "base"
    else:
        scenario_label = chosen_data.name
    dispatch_name = f"{base_case}_{scenario_label}_dispatch"

    # Step 4: create output folder for this dispatch run
    output_folder = dispatch_dir / scenario_label
    output_folder.mkdir(parents=True, exist_ok=True)

    # Step 5: run the dispatch simulation
    run_dispatch(chosen_data, chosen_gdx, output_folder, dispatch_name)

    # Step 6: final summary
    print("------------------------------------------------------------------------")
    print("Dispatch-only run completed successfully.")
    print(f"Results saved in: {output_folder}")
    print("------------------------------------------------------------------------")

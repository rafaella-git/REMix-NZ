# run_remix.py
import os, time
from pathlib import Path
import pandas as pd
from remix.framework.api.instance import Instance

# Notes:
#  1. RUN_TYPE='opt'        normal full optimization.
#  2. RUN_TYPE='scenario'   optimization using overrides from data/<scenario_name>/.
#  3. RUN_TYPE='dispatch'   dispatch-only simulation using fixed capacities from a GDX file.

# #### If the model is infeasible or a run crashes
# If you build a model on your own, it can always happen that it gets infeasible or the execution stops out of
# some other reason. In that case, you can refer to the `run_remix.lst` file and look for the error marker `****`.
# You can open that file in an editor and look for the error message.
# Alternatively, you can also use the commandline tool grep to search for the pattern "**** Exec Error"  in the file to see what is wrong.

group_name = "hadi"
base_case  = "pypsa"       # base case folder under ../project/{group_name}/{base_case}/data
low_case   = "pypsa-low"   # second case folder for dispatch run

# --- Folder setup -------------------------------------------------------------
base_dir = Path(f"../project/{group_name}/{base_case}")
low_dir  = Path(f"../project/{group_name}/{low_case}")

base_data   = base_dir / "data"
base_result = base_dir / "result"
low_data    = low_dir / "data"
low_result  = low_dir / "result"

base_result.mkdir(parents=True, exist_ok=True)
low_result.mkdir(parents=True, exist_ok=True)

# shared run options
run_opts = dict(
    solver    = "gurobi", #"cplex",
    threads   = 8,
    keep      = 1,                     # keep scratch folder “/225a” with exported .csvs for debugging
    lo        = 4,                     # write a .log file next to the .lst
    names     = 1,                     # use actual variable names in .lst in case of error
    roundts   = 1,                     # round profile files to avoid GAMS line length/precision issues
    timeres   = 1,                     # hourly; use 24 for daily aggregation, etc.
    postcalc  = 1,                     # run post-processing
    pathopt   = "myopic",              # keep if you use myopic/rolling runs
)

os.chdir(Path(__file__).parent.resolve())


# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

def run_opt_case(data_dir, result_dir, result_name):
    """Normal optimization run"""
    print(f"\n=== Running NORMAL optimization for case={data_dir.name} ===")
    m = Instance(index_names=False, datadir=data_dir)

    # FIX: do NOT pass resultdir twice — use run_opts.copy() and inject the directory there
    opts = run_opts.copy()
    opts.update(dict(resultdir=result_dir, resultfile=result_name))

    rc = m.run(**opts)
    if rc != 0:
        print(f"\nREMix returned non-zero code {rc}")
        raise SystemExit(rc)
    print(f"✔ Optimization completed successfully.")
    return result_dir / f"{result_name}.gdx"


def run_dispatch_case(data_dir, result_dir, result_name, fixed_gdx, freeze_expansion=True):
    """Dispatch-only run"""
    print(f"\n=== Running DISPATCH-ONLY for case={data_dir.name} ===")
    print(f"      Using capacities from {fixed_gdx}")

    m = Instance(index_names=False, datadir=data_dir)

    if freeze_expansion:
        try:
            nodes = list(m.set.nodesdata)
            years = list(m.set.yearssel)
            techs = list(m.set.converters)
            lock = pd.DataFrame(
                index=pd.MultiIndex.from_product([nodes, years, techs],
                                                 names=["region","years","technology"]),
                columns=["noExpansion"]
            )
            lock["noExpansion"] = 1
            m.parameter.add(lock, "converter_capacityparam")
            print("     Expansion frozen (noExpansion=1 for all converters)")
        except Exception as e:
            print(f"   ! Couldn’t apply noExpansion programmatically: {e}")

    opts = run_opts.copy()
    opts.update(dict(resultdir=result_dir, resultfile=result_name, fixedcapsfromgdx=str(fixed_gdx.as_posix())))

    rc = m.run(**opts)
    if rc != 0:
        print(f"\nREMix returned non-zero code {rc}")
        raise SystemExit(rc)
    print(f"✔ Dispatch completed successfully.")



t0 = time.perf_counter()

# 1 Run base optimization
opt_gdx = run_opt_case(base_data, base_result, f"{base_case}_opt")

# 2 Run dispatch using capacities from base run
run_dispatch_case(low_data, low_result, f"{low_case}_dispatch", fixed_gdx=opt_gdx)

t1 = time.perf_counter()
print(f"\n✔ All runs completed successfully.")
print(f"Total runtime: {time.strftime('%Hh %Mm %Ss', time.gmtime(t1 - t0))}")
print("Results written to:")
print(f"   Base case : {base_result}")
print(f"   Low case  : {low_result}")



# ---- Notes -------------------------------------------------------------------
# - Put only modified files under: ../project/{group_name}/{case_name}/data/<scenario>/ 
# - REMix loads base files from data/ and overrides with any file found in data/<scenario>/.

# #### Explanation for command line arguments to GAMS function call
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
# `crossover` : Run the optimization with crossover (1) or without crossover (0, default)
# `datacheck` : will make CPLEX give out warnings about disproportionate values (which lead to numerical difficulties = non-optimal solutions) at the beginning of the optimization, default 0

# `barcolnz` : Number of non-zero entries above which CPLEX solver will treat columns as dense, default 0 (which means parameter is determined dynamically)

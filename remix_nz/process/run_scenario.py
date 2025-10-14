# run_remix.py
import os, time
from pathlib import Path
from remix.framework.api.instance import Instance

group_name = "hadi"
case_name  = "pypsa-low" #"h2-domestic_2020-2030-2050"           # this is the *base* scenario directory under project/{group_name}/{case_name}/data
base_case  = "pypsa" 
# optional: list of sub-scenarios *inside* the base data folder to run (each subfolder overrides files from the base)
scenarios = [
    None,              # None = run the base data folder itself (no overrides)
    #"low",#"wind",            # e.g. data/wind/ contains a lower CAPEX file for wind
    # "batt-cheap",
    # "no-h2",
]

# Paths 
data_dir    = Path(f"../project/{group_name}/{case_name}/data")
results_dir = Path(f"../project/{group_name}/{case_name}/result")
results_dir.mkdir(parents=True, exist_ok=True)
fixed_caps_gdx = results_dir/f"{case_name}.gdx"

if not data_dir.exists():
    raise IOError(f"Build step missing: {data_dir} not found")

# one Instance is enough; we pass scendir per run -
m = Instance(index_names=False, datadir=data_dir)

# Small helper for a nice label in the result filename
def tag(scn):
    return "base" if (scn is None or str(scn).strip() == "") else str(scn).replace("/", "_")

# run options (see docs) 
run_args = dict(
    resultdir = results_dir,           # where remix.gdx will be written
    solver    = "gurobi", #"cplex",              # cplex/highs/mosek/xpress/scip also supported
    threads   = 8,
    keep      = 1,                     # keep scratch folder “/225a” with exported .csvs for debugging
    lo        = 4,                     # write a .log file next to the .lst
    names     = 1,                     # use actual variable names in .lst in case of error
    roundts   = 1,                     # round profile files to avoid GAMS line length/precision issues
    timeres   = 1,                     # hourly; use 24 for daily aggregation, etc.
    postcalc  = 1,                     # run post-processing
    iis     = 1,                     # write IIS .gdx on infeasibility (enable when needed)
    # equlist = 1,                     # huge listing of all equations (only for tiny models)
    # gdx     = "default",             # extra symbol dump for deep debugging
    # solvermethod = 1,                # 1=barrier (IPM), 2=simplex, 4=dual, etc.
    # scaletimefrac = 1,               # scale sourcesink annual sums/indicators if using partial year (timestart/timeend)
    # timestart = 1, timeend = 8760,   # limit time window (e.g. for short debug runs)
    fixedcapsfromgdx = "../project/{group_name}/{base_case}/{base_case}.gdx",  # read fixed capacities from a previous run (path relative to run dir)
    pathopt   = "myopic",               # keep if you use myopic/rolling runs
)

# Run base and any scenario overrides 
m = Instance(index_names=False, datadir=data_dir)
t0 = time.perf_counter()

for scn in scenarios:
    label = tag(scn)
    print(f"\n=== Running case={case_name} scenario={label} ===")

    run_kwargs = dict(run_args)
    run_kwargs["resultfile"] = f"{case_name}_{label}"

    # Scenario handling
    if scn is not None:
        scn_path = data_dir / scn
        if not scn_path.exists():
            raise IOError(f"Scenario folder not found: {scn_path}")
        run_kwargs["scendir"] = str(scn)  # override only the files in data/<scenario>

        # if this scenario is "low", import fixed capacities from base results
        if scn == "low":
            if not fixed_caps_gdx.exists():
                raise FileNotFoundError(f"Missing fixed capacity GDX: {fixed_caps_gdx}")
            run_kwargs["fixedcapsfromgdx"] = str(fixed_caps_gdx.as_posix())
            print(f"   Using fixed capacities from {fixed_caps_gdx}")

    # Run REMix
    rc = m.run(**run_kwargs)
    if rc != 0:
        print(f"REMix returned non-zero code {rc} for scenario={label}")

t1 = time.perf_counter()
print(f"\nAll runs finished in {time.strftime('%Hh %Mm %Ss', time.gmtime(t1 - t0))}")
print("Results written to:", results_dir)

# ---- Notes -------------------------------------------------------------------
# - Put only modified files under: ../project/{group_name}/{case_name}/data/<scenario>/
# - REMix loads base files from data/ and overrides with any file found in data/<scenario>/.
# - Keep keep=1 & lo=4 to retain exported CSVs and a detailed .log for debugging.

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


#solver arguments
# `crossover` : Run the optimization with crossover (1) or without crossover (0, default)
# `datacheck` : will make CPLEX give out warnings about disproportionate values (which lead to numerical difficulties = non-optimal solutions) at the beginning of the optimization, default 0
#`iis` : *.gdx file with all infeasible equations (if any), default 0
#`names` : Write out actual variable names into *.lst file in case of error, default 0
#`barcolnz` : .Number of non-zero entries above which CPLEX solver will treat columns as dense, default 0 (which means parameter is determined dynamically)


# #### If the model is infeasible or a run crashes
# If you build a model on your own, it can always happen that it gets infeasible or the execution stops out of
# some other reason. In that case, you can refer to the `run_remix.lst` file and look for the error marker `****`.
# You can open that file in an editor and look for the error message.
# Alternatively, you can also use the commandline tool grep to search for the pattern "**** Exec Error"  in the file to see what is wrong.

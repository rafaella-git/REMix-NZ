# run.py
# ------------------------------------------------------------------------
# Minimal REMix optimisation run (no scenarios)
# Runs the model for one case folder located under ../project/{group_name}/
# ------------------------------------------------------------------------

import os, time
from pathlib import Path
from remix.framework.api.instance import Instance

# -----------------------------------------------------------------------------
# user settings
group_name = "hadi"            # main folder under /project/
case_name  = "pypsa-cascade"   # name of the case folder (contains /data/ and /result/)
solver     = "gurobi"          # choose "gurobi", "cplex", etc.

# -----------------------------------------------------------------------------
# paths
os.chdir(Path(__file__).parent.resolve())
data_dir   = Path(f"../project/{group_name}/{case_name}/data")
result_dir = Path(f"../project/{group_name}/{case_name}/result")
result_dir.mkdir(parents=True, exist_ok=True)

# -----------------------------------------------------------------------------
# run REMix
start = time.perf_counter()
m = Instance(index_names=False, datadir=data_dir)

m.run(
    resultdir=result_dir,
    resultfile=f"{case_name}",
    solver=solver,
    datacheck=1,      # `datacheck` : will make CPLEX give out warnings about disproportionate values
    keep=1,            # keep scratch folder “/225a” with exported .csvs for debugging
    lo=4,              # 3 = only .lst file, 4 = .lst + .log file for solver output
    postcalc=1,        # run post-calculation after solving
    roundts=1,         # round time-series to avoid precision issues / numerical errors
    names=1,           # include variable names in .lst
    timeres=1,         # 1 = hourly; 24 = daily aggregation
    threads=8,         # number of CPU threads
    pathopt="myopic",  # foresight / myopic / target
)

end = time.perf_counter()
print(f"\nRun completed in {time.strftime('%Hh %Mm %Ss', time.gmtime(end-start))}")
print("Results written to:", result_dir)

# -----------------------------------------------------------------------------
# Notes:
# - Ensure the /data/ folder contains all input .csv/.dat files for REMix.
# - Results (.gdx, .lst, .log) are written to the /result/ folder.
# - pathopt options:
#     foresight : perfect foresight over the whole horizon
#     myopic    : limited foresight (e.g. stepwise years)
#     target    : single-year optimisation (last year in yearssel)
# -----------------------------------------------------------------------------

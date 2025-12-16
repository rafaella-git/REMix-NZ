# run_gp_mvp.py
# ------------------------------------------------------------------------
# Minimal REMix optimisation run for the NZ GP MVP case
# ------------------------------------------------------------------------

import os, time
from pathlib import Path
from remix.framework.api.instance import Instance

# -----------------------------------------------------------------------------
# user settings
group_name = "GP-NT-ELEC-BIO-H2"   # main folder under /project/
case_name  = "nz_case_GP_MVP"      # case folder (has /data/ and /result/)
solver     = "cplex"#"gurobi"              # or "cplex", etc.

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
    datacheck=1,
    keep=1,
    lo=4,
    postcalc=1,
    roundts=1,
    names=1,
    timeres=1,          # hourly
    threads=8,
    pathopt="myopic",   # or "myopic"/"foresight" /"target"
)

end = time.perf_counter()
print(f"\nRun completed in {time.strftime('%Hh %Mm %Ss', time.gmtime(end-start))}")
print("Results written to:", result_dir)

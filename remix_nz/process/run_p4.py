# run_all_cases.py
import os, time
from pathlib import Path
from remix.framework.api.instance import Instance

# List all models to run
cases = [

    # -------------------------
    # 2020-2050
    # -------------------------
    
    # ("GP-NT-ELEC-BIO-H2", "nz_case_GP_2020-2050", "cplex"),
    # ("GP-NT-ELEC-BIO-H2", "nz_case_NT_2020-2050", "cplex"),
    # ("GP-NT-ELEC-BIO-H2", "nz_case_ELEC+_2020-2050", "cplex"),
    # ("GP-NT-ELEC-BIO-H2", "nz_case_BIO+_2020-2050", "cplex"),
    # ("GP-NT-ELEC-BIO-H2", "nz_case_H2+_2020-2050", "cplex"),

    # -------------------------
    #  2020-2025-2030-2035-2040-2045-2050
    # -------------------------
    ("GP-NT-ELEC-BIO-H2", "nz_case_GP_2020-2025-2030-2035-2040-2045-2050", "gurobi"),
    ("GP-NT-ELEC-BIO-H2", "nz_case_NT_2020-2025-2030-2035-2040-2045-2050", "gurobi"),
    ("GP-NT-ELEC-BIO-H2", "nz_case_ELEC+_2020-2025-2030-2035-2040-2045-2050", "gurobi"),
    ("GP-NT-ELEC-BIO-H2", "nz_case_BIO+_2020-2025-2030-2035-2040-2045-2050", "gurobi"),
    ("GP-NT-ELEC-BIO-H2", "nz_case_H2+_2020-2025-2030-2035-2040-2045-2050", "gurobi"),
]


# -----------------------------------------------------------------------------
# base path
script_dir = Path(__file__).parent.resolve()
os.chdir(script_dir)

for group_name, case_name, solver in cases:
    print("\n" + "-" * 80)
    print(f"Running case: group='{group_name}', case='{case_name}'")
    print("-" * 80)

    data_dir   = script_dir / ".." / "project" / group_name / case_name / "data"
    result_dir = script_dir / ".." / "project" / group_name / case_name / "result"
    result_dir.mkdir(parents=True, exist_ok=True)

    if not data_dir.exists():
        print(f"  [SKIP] Data directory does not exist: {data_dir}")
        continue

    start = time.perf_counter()
    try:
        m = Instance(index_names=False, datadir=data_dir)

        m.run(
            resultdir=result_dir,
            resultfile=f"{case_name}",
            solver=solver,
            datacheck=1,
            # keep=1,
            lo=4,
            postcalc=1,
            roundts=1,
            names=1,
            timeres=1,      # hourly
            threads=8,
            pathopt="myopic", # "myopic"/"foresight" /"target"
        )

        status = "OK"
    except Exception as e:
        status = f"FAILED ({type(e).__name__}: {e})"

    end = time.perf_counter()
    print(
        f"  Status: {status}\n"
        f"  Runtime: {time.strftime('%Hh %Mm %Ss', time.gmtime(end-start))}\n"
        f"  Results dir: {result_dir}"
    )

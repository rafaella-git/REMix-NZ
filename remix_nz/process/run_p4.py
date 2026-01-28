import os
import time
import shutil
from datetime import datetime
from pathlib import Path

from remix.framework.api.instance import Instance

# -------------------------
# user settings
# -------------------------
save_tempfiles = True

process_dir = Path(r"C:\Local\REMix\remix_nz\process")
lst_src = process_dir / "run_remix.lst"
log_src = process_dir / "run_remix.log"

temp_dir = Path(r"C:\Local\REMix\remix_nz\project\GP-NT-ELEC-BIO-H2\temp_files")
temp_dir.mkdir(parents=True, exist_ok=True)

# -------------------------
# cases
# -------------------------
cases = [
    ("GP-NT-ELEC-BIO-H2", "nz_case_BIO+_2020-2050", "gurobi"),
    #("GP-NT-ELEC-BIO-H2", "nz_case_H2+_2020-2025-2030-2035-2040-2045-2050", "cplex"),
    #("GP-NT-ELEC-BIO-H2", "nz_case_BIO+_2020-2025-2030-2035-2040-2045-2050", "cplex"),
    #("GP-NT-ELEC-BIO-H2", "nz_case_ELEC+_2020-2025-2030-2035-2040-2045-2050", "cplex"),
    #("GP-NT-ELEC-BIO-H2", "nz_case_NT_2020-2025-2030-2035-2040-2045-2050", "cplex"),
    #("GP-NT-ELEC-BIO-H2", "nz_case_GP_2020-2025-2030-2035-2040-2045-2050", "cplex"),
]

# -----------------------------------------------------------------------------
script_dir = Path(__file__).parent.resolve()
os.chdir(script_dir)

for group_name, case_name, solver in cases:
    print("\n" + "-" * 80)
    print(f"Running case: group='{group_name}', case='{case_name}'")
    print("-" * 80)

    data_dir = script_dir / ".." / "project" / group_name / case_name / "data"
    result_dir = script_dir / ".." / "project" / group_name / case_name / "result"
    result_dir.mkdir(parents=True, exist_ok=True)

    if not data_dir.exists():
        print(f"  [skip] Data directory does not exist: {data_dir}")
        continue

    start = time.perf_counter()

    try:
        m = Instance(index_names=False, datadir=data_dir)

        m.run(
            resultdir=result_dir,
            resultfile=case_name,
            solver=solver,
            datacheck=1,
            lo=4,
            postcalc=1,
            roundts=1,
            names=1,
            timeres=1,      # hourly
            threads=8,
            pathopt="myopic",  # "myopic"/"foresight" /"target"
        )

        if save_tempfiles:
            ts = datetime.now().strftime("%Y%m%d-%H%M%S")
            lst_dst = temp_dir / f"{case_name}__{ts}.lst"
            log_dst = temp_dir / f"{case_name}__{ts}.log"

            if lst_src.exists():
                shutil.copy2(lst_src, lst_dst)  # preserves metadata where possible [web:80]
            else:
                print(f"  [warn] Missing .lst: {lst_src}")

            if log_src.exists():
                shutil.copy2(log_src, log_dst)  # preserves metadata where possible [web:80]
            else:
                print(f"  [warn] Missing .log: {log_src}")

        status = "OK"

    except Exception as e:
        status = f"FAILED ({type(e).__name__}: {e})"

    end = time.perf_counter()
    print(
        f"  Status: {status}\n"
        f"  Runtime: {time.strftime('%Hh %Mm %Ss', time.gmtime(end-start))}\n"
        f"  Results dir: {result_dir}"
    )
    if save_tempfiles:
        print(f"  Temp files dir: {temp_dir}")

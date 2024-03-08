# %% [markdown]
# ## Part b: running the model

# %%
# reading model built in Part a
from remix.framework.api.instance import Instance
import pathlib as pt
import os
import time

# %%
# Changing years to be optimised and demand file selected
yrs2run=[2020,2035,2050] # years to be optimised


# Demand files available (different scenarios)
files_lst=["01-battery-distributed","02-battery-overnight", "03-battery-recharging", "04-h2-solar",  "01-h2-distributed", "02-h2-overnight", "03-h2-recharging", "04-battery-solar"]


def run_single(file_list=files_lst):
    yrs_str='-'.join([str(item) for item in yrs2run])
    indx=0
    for file in files_lst:
        indx+=1
        demand_file=files_lst[indx] 
        case_name=f"{demand_file}_{yrs_str}"
        output_dir = pt.Path(f"C:/Local/REMix/remix_nz/output/will/{case_name}")
        data_dir = pt.Path(f"{output_dir}/data")
        result_dir = pt.Path(f"{output_dir}/result")

        if not data_dir.exists():
            raise IOError("You need to run will_build_instance")
        m = Instance.from_path(data_dir)

        # %%
        # running GAMS from Python script
        s1 = time.perf_counter()

        m.run(
            resultdir = result_dir,
            resultfile = case_name,
            lo=4,
            timeres=1,
            names=1,
            roundts=1,
            iis=1,
        )
        print(os. getcwd())

        e1 = time.perf_counter()
        d1=time.strftime("%Hh %Mm %Ss", time.gmtime(e1-s1))
        print(f"------------- Running the model took {d1}.")
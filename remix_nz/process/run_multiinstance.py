# %% [markdown]
# ## Part B: running the model

# Import libraries
from remix.framework.api.instance import Instance
import pathlib as pt
import os
import time

# %% [markdown]
# Define relevant lists for cases to run (years to be considered and demand file)
yrs2run=[2020,2030,2040,2050] # years to be optimised
# Demand files available (different scenarios)
files_lst=["nz-elec","nz-elec","nz_profile_11nodes"] 
yrs2run=[yrs2run]*len(files_lst)
files_493=["medpop_evs","low_pop_out","med_pop_out","high_pop_out"] 


# %% [markdown]
# Define function that reads the model built in `part a`, it takes:
# 1. A list of demand files to access 
# 2. A list of years considered in the creation of the data
# 3. An optional comment fir the resulting file 
def run_multi(files=files_lst, years=yrs2run, timeresolution=1, comment=""):
    s2 = time.perf_counter()
    print(f"================ {len(files)} cases will be run.")
    if not len(files)==len(years):
        years=[years]*len(files) #if only one year list is provided, it uses that one for all the files
    indx=0
    for case in files:
        yrs_str='-'.join([str(item) for item in years[indx]])
        demand_file=files[indx] 
        case_name=f"{demand_file}_{yrs_str}"
        output_dir = pt.Path(f"C:/Local/REMix/remix_nz/output/{case_name}")
        data_dir = pt.Path(f"{output_dir}/data")
        result_dir = pt.Path(f"{output_dir}/result")
        if not data_dir.exists():
            raise IOError("You need to build  the model first! (create the data)")
        m = Instance.from_path(data_dir)
        s1 = time.perf_counter()
        print(f"------------- Currently running {case} for the year(s) {yrs_str}.")
        m.run(
            resultdir=result_dir,
            resultfile=f"{case_name}{comment}",
            lo=4,
            names=1,
            timeres=3,
            roundts=1,
        )
        print(os. getcwd())
        e1 = time.perf_counter()
        d1=time.strftime("%Hh %Mm %Ss", time.gmtime(e1-s1))
        print(f"------------- Running {case} for the year(s) {yrs_str} took {d1}.")
        indx+=1
    e2 = time.perf_counter()
    d2=time.strftime("%Hh %Mm %Ss", time.gmtime(e2-s2))
    print(f"------------- Running all the files in the list took {d2}.")



# %% [markdown]
# ### Run with perfect foresight
# Since perfect foresight is the default run,
# when calling `run_multi` we run perfect foresight
# for all the defined cases
run_multi(["nz_profile_11nodes","nz_profile_11nodes","nz-elec","nz-elec"],[[2020],[2020,2030,2040,2050],[2050],[2020,2030,2040,2050]],3)
# %%

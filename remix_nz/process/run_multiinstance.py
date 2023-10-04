# %% [markdown]
# ## Part b: running the model

# %%
# ### Reading model built in Part a (data created)
from remix.framework.api.instance import Instance
import pathlib as pt
import os
import time


# %%
# Select years to be optimised and data for the run
yrs2run=[2020,2030,2040,2050] # years to be optimised
# Demand files available (different scenarios)
files_lst=["nz_profile_11nodes"] 
files_493=["medpop_evs_base","low_pop_out_base","med_pop_out_base","high_pop_out_base"] #493 as in the course 493, this is for Liv and Sam


def run_multi(files=["nz_profile_11nodes"], years=[2050]):
    yrs_str='-'.join([str(item) for item in years])
    s2 = time.perf_counter()
    indx=0
    for case in files:
        demand_file=files[indx] 
        case_name=f"{demand_file}_{yrs_str}"
        output_dir = pt.Path(f"C:/Local/REMix/remix_nz/output/{case_name}")
        data_dir = pt.Path(f"{output_dir}/data")
        result_dir = pt.Path(f"{output_dir}/result")
        if not data_dir.exists():
            raise IOError("You need to build  the model first! (create the data)")
        m = Instance.from_path(data_dir)
        s1 = time.perf_counter()
        m.run(
            resultdir=result_dir,
            resultfile=case_name,
            lo=4,
            names=1,
            roundts=1,
        )
        print(os. getcwd())
        e1 = time.perf_counter()
        d1=time.strftime("%Hh %Mm %Ss", time.gmtime(e1-s1))
        print(f"------------- Running {case} for the years {yrs_str} took {d1}.")
        indx+=1
    e2 = time.perf_counter()
    d2=time.strftime("%Hh %Mm %Ss", time.gmtime(e2-s2))
    print(f"------------- Running all the files in the list for the years {yrs_str} took {d2}.")

# %%
# ### Running multiple cases
#run_multi(files_493,yrs2run)

yrs_str='-'.join([str(item) for item in yrs2run])
s2 = time.perf_counter()
files_493=["medpop_evs_base","low_pop_out_base","med_pop_out_base","high_pop_out_base"] #493 as in the course 493, this is for Liv and Sam
indx=2
demand_file=files_493[indx] 
case_name=f"{demand_file}_{yrs_str}"
output_dir = pt.Path(f"C:/Local/REMix/remix_nz/output/{case_name}")
data_dir = pt.Path(f"{output_dir}/data")
result_dir = pt.Path(f"{output_dir}/result")
if not data_dir.exists():
    raise IOError("You need to build  the model first! (create the data)")
m = Instance.from_path(data_dir)
s1 = time.perf_counter()
m.run(
    resultdir=result_dir,
    resultfile=case_name,
    lo=4,
    names=1,
    roundts=1,
    )
print(os. getcwd())
e1 = time.perf_counter()
d1=time.strftime("%Hh %Mm %Ss", time.gmtime(e1-s1))
print(f"------------- Running {case} for the years {yrs_str} took {d1}.")
        
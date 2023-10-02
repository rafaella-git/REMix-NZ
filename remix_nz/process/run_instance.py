# %% [markdown]
# ## Part b: running the model

# %%
# reading model built in Part a
from remix.framework.api.instance import Instance
import pathlib as pt
import os
import time

#borrar os.chdir("C:/Local/REMix/remix_nz/process")


# %%
# Changing years to be optimised and demand file selected
yrs2run=[2050] # years to be optimised
indx=0
# Demand files available (different scenarios)
files_lst=["nz_profile_11nodes","medpop_evs","low_pop_out","med_pop_out","high_pop_out"] 
files_493=["medpop_evs_base","low_pop_out_base","med_pop_out_base","high_pop_out_base"] #493 as in the course 493, this is for Liv and Sam
# Do not change these
yrs_str='-'.join([str(item) for item in yrs2run])
demand_file=files_lst[indx] 
case_name=f"{demand_file}_{yrs_str}"
output_dir = pt.Path(f"C:/Local/REMix/remix_nz/output/{case_name}")
data_dir = pt.Path(f"{output_dir}/data")
result_dir = pt.Path(f"{output_dir}/result")

if not data_dir.exists():
    raise IOError("You need to run tutorial 101a first!")


m = Instance.from_path(data_dir)

# %%
# running GAMS from Python script
s1 = time.perf_counter()

m.run(
    resultdir = result_dir,
    resultfile=case_name,
    lo=4,
    names=1,
    roundts=1,
)
print(os. getcwd())

e1 = time.perf_counter()
d1=time.strftime("%Hh %Mm %Ss", time.gmtime(e1-s1))
print(f"------------- Running the model took {d1}.")
# %% [markdown]
# ## Part b: running the model

# %%
# reading in model built in Part a`
from remix.framework.api.instance import Instance
import pathlib as pt
import os

# %%
# Changing years to be optimised and demand file selected
yrs2run=[2050]
indx=0
# Demand files available (different scenarios)
files_lst=["nz_profile_11nodes","medpop_evs","low_pop_out","med_pop_out","high_pop_out"] 
files_493=["medpop_evs_base","low_pop_out_base","med_pop_out_base","high_pop_out_base"] #493 as in the course 493, this is for Liv and Sam
# Do not change these
print(os. getcwd())
yrs_str='-'.join([str(item) for item in yrs2run])
demand_file=files_lst[indx] 
case_name=f"{demand_file}_{yrs_str}"
output_dir = pt.Path(f"C:/Local/REMix/remix_nz/output/{case_name}")

data_dir = pt.Path(f"{output_dir}/data")
result_dir = pt.Path(f"{output_dir}/result").mkdir(parents=False, exist_ok=True)
print(os. getcwd())

if not data_dir.exists():
    raise IOError("You need to run tutorial 101a first!")

#os.chdir("../")
m = Instance.from_path(data_dir)
#print(os. getcwd())

# %%
# running GAMS from Python script
m.run(
    resultdir = result_dir,
    resultfile=case_name,
    lo=4,
    names=1,
    roundts=1,
)
print(os. getcwd())

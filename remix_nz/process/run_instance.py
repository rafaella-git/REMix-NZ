# %% [markdown]
# ## Part b: running the model

# %%
# reading model built in Part a
from remix.framework.api.instance import Instance
import pathlib as pt
import os
import time
from pathlib import Path

#borrar os.chdir("C:/Local/REMix/remix_nz/process")


# %%
# Changing years to be optimised and demand file selected
# yrs2run=[2020,2050] # years to be optimised
# indx=0
# # Demand files available (different scenarios)
# files_lst=["nz-h2-stor","nz-h2-noexp","nz-h2-v2","nz-h2", "nz-hydro","nz-elec","nz_profile_11nodes","medpop_evs","low_pop_out","med_pop_out","high_pop_out"] 
# # files_493=["medpop_evs_base","low_pop_out_base","med_pop_out_base","high_pop_out_base"] #493 as in the course 493, this is for Liv and Sam
# # Do not change these
# yrs_str='-'.join([str(item) for item in yrs2run])
# demand_file=files_lst[indx] 
# case_name=f"{demand_file}_{yrs_str}"
# output_dir = pt.Path(f"C:/Local/REMix/remix_nz/output/{case_name}")
# data_dir = pt.Path(f"{output_dir}/data")
# result_dir = pt.Path(f"{output_dir}/result")
scenario_dict = {       
    "will": [["00-test-elec"], [2020]],
    #"will": [will_lst, [2020, 2035, 2050]],
    #"sdewes-ap": [sdewes_lst, [2020, 2030, 2040, 2050]]
}
group_name="will"
files_lst = scenario_dict[group_name][0]
yrs2run = scenario_dict[group_name][1] 
yrs_str='-'.join([str(item) for item in yrs2run])
yrs_mentioned = [2000, 2020] #[2000, 2020, 2025, 2030, 2035, 2040, 2045, 2050]
indx=0


# Define paths/directories
path_base = "C:/Local/REMix"
path_input = f"{path_base}/remix_nz/input"
path_demand = f"{path_input}/demand/{group_name}"
path_profiles = f"{path_input}/profiles"      # renewables
path_brownfield = f"{path_input}/brownfield"  # info hydro and existing power plants database
path_geo = f"{path_input}/shapefiles"         # geojson
geofile="11regionsNZ.geojson"
demand_file=files_lst[indx] 
case_name=f"{demand_file}_{yrs_str}"
path_output = f"{path_base}/remix_nz/output/{group_name}"
data_dir = Path(f"{path_output}/{case_name}/data").mkdir(parents=True, exist_ok=True)
results_dir = Path(f"{path_output}/{case_name}/result").mkdir(parents=True, exist_ok=True)

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
    timeres=1,
    names=1,
    roundts=1,
    iis=1,
)
print(os. getcwd())

e1 = time.perf_counter()
d1=time.strftime("%Hh %Mm %Ss", time.gmtime(e1-s1))
print(f"------------- Running the model took {d1}.")
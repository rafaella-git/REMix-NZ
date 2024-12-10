# %% [markdown]
# ## Running the model
#
# The written `*.dat` files are used as inputs to the GAMS optimization.
# To successfully run this file, you first need to build the model.


# %%
# Reading in model built previously
from remix.framework.api.instance import Instance
from pathlib import Path


group_name="dlr"
indx=0

# ### Define the years to run the optimisation and the demand file
# The demand file and the years run determine the name of the case and its results

# Defining demand file options per folder
europe=["h2-lut-domestic", "h2-lut-exports", "h2-lut-exports-v2", "h2-pypsa","h2-pypsa-exports-domestic", "h2-pypsa-exports-20","h2-pypsa-exports-40","h2-pypsa-exports-200"]
dlr=["h2-domestic"]

folder_dict = { 
    "europe": [europe, [2020, 2030,2050]],
    "dlr": [dlr, [2020, 2030,2050]]
}

files_lst = scenario_dict[group_name][0]
yrs_sel = scenario_dict[group_name][1] # [2020, 2025, 2030, 2035, 2040, 2045, 2050]
yrs_str='-'.join([str(item) for item in yrs_sel])
yrs_to_calc = [2020, 2025, 2030, 2035, 2040, 2045, 2050]



case_name=f"dlr_2020-2030-2050"
data_dir = Path(f"../output/{scenario_dict[group_name][0]}/{case_name}/data")

# Defining the directory the model data is written in (folder "data/" in the project directory)
if not data_dir.exists():
    raise IOError("You need to build the model first!")
m = {i: Instance(datadir=data_dir,index_names=False) for i in ["wind"]}
m["Base"] = Instance(index_names=False,datadir=data_dir)


# %% [markdown]
# ### Executing optimization model in GAMS
#
# Here the model is executed and optimized using GAMS, using the function `run` from the Instance. 

# %%
# running GAMS from Python script
scen = "wind+"

s1 = time.perf_counter()
m["Base"].run(
    resultdir=results_dir,
    #resultfile=f"{case_name}",
    resultfile=f"{case_name}_{scen.replace('/', '_')}",
    scendir=f"{scen}",
    solver="gurobi",
    threads=8,
    lo=4,
    timeres=1,
    names=1,
    roundts=1,
    iis=1, # Irreducible Infeasible Subsystem
    # profile=1,
    # gdx="default",
    pathopt="myopic",
    barorder=1
)
print(os. getcwd())

e1 = time.perf_counter()
d1=time.strftime("%Hh %Mm %Ss", time.gmtime(e1-s1))
print(f"------------- Running {case_name} took {d1}.")




# %% [markdown]
# #### Explanation for command line arguments to GAMS function call
#
# `lo=3` : log option of GAMS; ensures that the output from GAMS will be visible in the terminal (`lo=4`
#          will write out a log file)
#
# `modeldir` : path to model directory (can be set relative or absolute)
#
# `datadir` : instruction to look for `*.dat` files in a specific directory
#
# `scendir` : instruction to look for `*.dat` files in a specific directory; REMix default: datadir
#
# `resultdir` : instruction to store the resulting `*.gdx` file in a specific directory
#
# `resultfile` : allows to give a custom name to the output `*.gdx` file; REMix default: "remix"
#
# `postcalc` : instruction whether to do a post-calculation; REMix default: 1
#
# `names` : instruction to give out actual variable names in case of error messages; REMix default: 0
#
# `roundts` : instruction to automatically round after-comma digits in large time series files where necessary to
# successfully read them in; REMix default: 0
#
# `keep` : keeps temporary folder with `*.csv` files built during the model run, might be interesting for debugging
# purposes; REMix default: 0
#
#
# #### If the model is infeasible or a run crashes
#
# If you build a model on your own, it can always happen that it gets infeasible or the execution stops out of
# some other reason. In that case, you can refer to the `run_remix.lst` file and look for the error marker `****`.
#
# You can open that file (it will be stored in the model_dir) in an editor and look for the error message.
# Alternatively, you can also use the commandline tool grep to search for the pattern "**** Exec Error"  in the
# file to see what is wrong.
#
#
# ### Evaluation of results
#
# This is everything needed to run a REMix model. If you are interested in an example on how to evaluate and interpret
# results, head over to part c of this tutorial. In part d, the bonus tasks, you will get an introduction to GAMS error
# handling and model extensions.

# %% [markdown]
# ## Running the model
#
# The written `*.dat` files are used as inputs to the GAMS optimization.
# To successfully run this file, you first need to build the model.


# %%
# Import dependencies
import os
import time
from remix.framework.api.instance import Instance
from pathlib import Path





# ### Define the years to run the optimisation and the demand file
# The demand file and the years run determine the name of the case and its results

# Defining demand file options per folder
europe=["h2-lut-domestic", "h2-lut-exports", "h2-lut-exports-v2", "h2-pypsa","h2-pypsa-exports-domestic", "h2-pypsa-exports-20","h2-pypsa-exports-40","h2-pypsa-exports-200"]
dlr=["h2-domestic"]
paper2=["no-h2"]
madison=["base_input"]

folder_dict = { 
    "europe": [europe, [2020, 2030,2050]],
    "dlr": [dlr, [2020, 2030,2050]],
    "paper2": [paper2, [2020, 2030,2050]],
    "madison": [madison, [2020, 2030,2050]]
}



group_name="madison"
case_name=f"base_input"#"separate-demand"
scenario = "base_input"#"not specified"# "wind"


# Defining the directory the model data is written in (folder "data/" in the project directory)
data_dir = Path(f"../project/{group_name}/{case_name}/data")
results_dir = Path(f"../project/{group_name}/{case_name}/result")
if not data_dir.exists():
    raise IOError("You need to build the model first!")

# %% [markdown]
# Executing optimization model in GAMS
s1 = time.perf_counter()
m = {i: Instance(datadir=data_dir,index_names=False) for i in ["wind"]}
m["Base"] = Instance(index_names=False,datadir=data_dir)
m["Base"].run(
    resultdir=results_dir,
    resultfile=f"{case_name}",#_{scenario.replace('/', '_')}",
    #scendir=f"{scenario}",
    solver="gurobi",#"cplex",
    threads=8,
    lo=4,
    timeres=1,
    names=1,
    roundts=1,
    barrier=0,
    iis=1, # Irreducible Infeasible Subsystem
    # profile=1,
    # gdx="default",
    pathopt="myopic",
    barorder=1
)
print(os. getcwd())
e1 = time.perf_counter()
d1=time.strftime("%Hh %Mm %Ss", time.gmtime(e1-s1))
print(f"------------- Running {case_name} (scenario {scenario}) took {d1}.")

# #### Explanation for command line arguments to GAMS function call
# `lo=3` : log option of GAMS; ensures that the output from GAMS will be visible in the terminal (`lo=4`
#          will write out a log file)
# `modeldir` : path to model directory (can be set relative or absolute)
# `datadir` : instruction to look for `*.dat` files in a specific directory
# `scendir` : instruction to look for `*.dat` files in a specific directory; REMix default: datadir
# `resultdir` : instruction to store the resulting `*.gdx` file in a specific directory
# `resultfile` : allows to give a custom name to the output `*.gdx` file; REMix default: "remix"
# `postcalc` : instruction whether to do a post-calculation; REMix default: 1
# `names` : instruction to give out actual variable names in case of error messages; REMix default: 0
# `roundts` : instruction to automatically round after-comma digits in large time series files where necessary to successfully read them in; REMix default: 0
# `keep` : keeps temporary folder with `*.csv` files built during the model run, might be interesting for debugging purposes; REMix default: 0
#

# #### If the model is infeasible or a run crashes
# If you build a model on your own, it can always happen that it gets infeasible or the execution stops out of
# some other reason. In that case, you can refer to the `run_remix.lst` file and look for the error marker `****`.
# You can open that file in an editor and look for the error message.
# Alternatively, you can also use the commandline tool grep to search for the pattern "**** Exec Error"  in the file to see what is wrong.
# %%

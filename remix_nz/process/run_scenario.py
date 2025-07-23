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
hadi=["pypsa"]

dict = { 
    "europe": [europe, [2020, 2030,2050]],
    "dlr": [dlr, [2020, 2030,2050]],
    "paper2": [paper2, [2020, 2030,2050]],
    "madison": [madison, [2020, 2030,2050]],
    "hadi": [hadi, [2020, 2030]]
}



group_name="hadi"
case_name=f"pypsa"#"separate-demand"
#scenario = "hydro-eq"#"not specified"# "wind"


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
    #fixedcapsfromgdx = f"{results_dir}/{case_name}",
    solver="gurobi",#"cplex", # to debug an unfeasible problem, run with CPLEX, you will need a license and select barrier 1
    threads=8,
    lo=4,
    timeres=1,
    names=1,
    roundts=1,
    barrier=1,
    iis=1, # Irreducible Infeasible Subsystem
    # profile=1,
    # gdx="default",
    pathopt="myopic",
    barorder=1,

)
print(os. getcwd())
e1 = time.perf_counter()
d1=time.strftime("%Hh %Mm %Ss", time.gmtime(e1-s1))
print(f"------------- Running {case_name} took {d1}.")# (scenario {scenario}) 

# #### Explanation for command line arguments to GAMS function call
# `lo=3` : log option of GAMS; ensures that the output from GAMS will be visible in the terminal (`lo=4`
#          will write out a log file)
# `modeldir` : path to model directory (can be set relative or absolute)
# `datadir` : instruction to look for `*.dat` files in a specific directory
# `scendir` : instruction to look for `*.dat` files in a specific directory; REMix default: datadir. Path to the scenario to calculate, relative to the --datadir path
# `solver` : Solver to use for optimising the model; supported options are "copt", "cplex", "gurobi", "highs", "mosek", "xpress", "convert", and "scip"; default: "cplex"

# `solvermethod` : Method of chosen solver to use for optimising the model; supported options are "automatic" (0), "barrier" (1), "simplex" (2), "primal" (3), "dual" (4), "concurrent" (5), "network" (6), "sifting" (7), pdlp (8), and mip (9); there is no solver that would support all these methods, if an unavailable solver method is specified, a hint will be given in the solve log file and the default will be used instead; default: "barrier"/"ipm"
# `resultdir` : instruction to store the resulting `*.gdx` file in a specific directory
# `resultfile` : allows to give a custom name to the output `*.gdx` file; REMix default: "remix"
# `postcalc` : instruction whether to do a post-calculation; REMix default: 1
# `names` : instruction to give out actual variable names in case of error messages; REMix default: 0
# `roundts` : instruction to automatically round after-comma digits in large time series files where necessary to successfully read them in; REMix default: 0
# `keep` : keeps temporary folder with `*.csv` files built during the model run, might be interesting for debugging purposes; REMix default: 0


#solver arguments
# `crossover` : Run the optimization with crossover (1) or without crossover (0, default)
# `datacheck` : will make CPLEX give out warnings about disproportionate values (which lead to numerical difficulties = non-optimal solutions) at the beginning of the optimization, default 0
#`iis` : *.gdx file with all infeasible equations (if any), default 0
#`names` : Write out actual variable names into *.lst file in case of error, default 0
#`barcolnz` : .Number of non-zero entries above which CPLEX solver will treat columns as dense, default 0 (which means parameter is determined dynamically)


# #### If the model is infeasible or a run crashes
# If you build a model on your own, it can always happen that it gets infeasible or the execution stops out of
# some other reason. In that case, you can refer to the `run_remix.lst` file and look for the error marker `****`.
# You can open that file in an editor and look for the error message.
# Alternatively, you can also use the commandline tool grep to search for the pattern "**** Exec Error"  in the file to see what is wrong.
# %%

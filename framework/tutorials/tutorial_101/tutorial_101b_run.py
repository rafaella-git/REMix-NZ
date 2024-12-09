# %% [markdown]
# ## Part b: running the model
#
# In the first script of this tutorial (part a), we configured our model and
# wrote its content to files that can be read in by GAMS.
# In the following steps, the written files are used as inputs to the
# GAMS optimization.
#
# N.B.: part a of this tutorial needs to be executed before this part can be
# successfully run.

# %%
# reading in model built in `tutorial_101a_build`
from remix.framework import Instance
import pathlib as pt

_path_tut101_data = pt.Path("../tutorial_101/data")

if not _path_tut101_data.exists():
    raise IOError("You need to run tutorial 101a first!")

m = Instance.from_path(_path_tut101_data)
# %% [markdown]
# ### Executing optimization model in GAMS
#
# In this section, the model we built above is finally executed and optimized
# using GAMS.
# We do that from this Python script using the function `run` from the Instance.
# The model can also be optimized from a terminal of choice. To do this instead,
# open up a terminal and call `run_remix` AFTER writing the files with part a of
# the tutorial.
#
# ```
# cd /path/to/tutorial_101/
# run_remix
# ```
#
# The (optional) arguments used in the GAMS call are explained below.

# %%
# running GAMS from Python script
m.run(
    resultfile="tutorial_101",
    lo=3,
    postcalc=1,
    roundts=1,
)
# %% [markdown]
# #### Explanation for command line arguments to GAMS function call
#
# `lo=3` : log option of GAMS; ensures that the output from GAMS will be visible
#          in the terminal (`lo=4` will write out a log file)
#
# `modeldir` : path to model directory (can be set relative or absolute)
#
# `datadir` : instruction to look for data files in a specific directory
#
# `scendir` : instruction to look for data files in a specific directory;
#             REMix default: datadir
#
# `resultdir` : instruction to store the resulting `*.gdx` file in a specific
#               directory
#
# `resultfile` : allows to give a custom name to the output `*.gdx` file;
#                REMix default: "remix"
#
# `postcalc` : instruction whether to do a post-calculation; REMix default: 1
#
# `names` : instruction to give out actual variable names in case of error
#           messages; REMix default: 0
#
# `roundts` : instruction to automatically round after-comma digits in large
#             time series files where necessary to successfully read them in;
#             REMix default: 0
#
# `keep` : keeps temporary folder with `*.csv` files built during the model run,
#          might be interesting for debugging purposes; REMix default: 0
#
#
# #### If the model is infeasible or a run crashes
#
# If you build a model on your own, it can always happen that it gets infeasible
# or the execution stops out of some other reason.
# In that case, you can refer to the `run_remix.lst` file and look for the error
# marker `****`.
#
# You can open that file (it will be stored in the model_dir) in an editor and
# look for the error message.
# Alternatively, you can also use the commandline tool grep to search for the
# pattern "**** Exec Error" in the file to see what is wrong.
#
#
# ### Evaluation of results
#
# This is everything needed to run a REMix model. If you are interested in an
# example on how to evaluate and interpret results, head over to part c of this
# tutorial.
# In part d, the bonus tasks, you will get an introduction to GAMS error
# handling and model extensions.

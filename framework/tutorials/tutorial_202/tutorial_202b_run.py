# %% [markdown]
# ## Part b: running the model
#
# For further explanations on what happens in this part of the model execution,
# refer to `tutorial_101b_run`.
# %%
# loading model built in `tutorial_202a_build`
from remix.framework.api.instance import Instance
import pathlib as pt
from pandas import IndexSlice as idx
from numpy import inf

_path_tut202_data = pt.Path("../tutorial_202/data")

if not _path_tut202_data.exists():
    raise IOError("You need to run tutorial 202a first!")

m = Instance.from_path(_path_tut202_data)

# %% [markdown]
# To evaluate the effect of the carbon budget, we set up three runs,
# one with a tight budget and one with loose budget and one with infinite
# carbon budget (i.e. no emission limit).
#
# %%
# running GAMS from Python script
m.run(
    resultfile="tutorial_202_loose",
    lo=3,
    names=1,
    postcalc=1,
    roundts=1,
)

m.parameter.accounting_indicatorbounds.loc[
    idx["global", "horizon", "Carbon"], "upperValue"
] = 100000.0
m.write(fileformat="dat")
m.run(
    resultfile="tutorial_202_tight",
    lo=3,
    names=1,
    postcalc=1,
    roundts=1,
)
# %%
m.parameter.accounting_indicatorbounds.loc[
    idx["global", "horizon", "Carbon"], "upperValue"
] = inf
m.write(fileformat="dat")
m.run(
    resultfile="tutorial_202_infinite",
    lo=3,
    names=1,
    postcalc=1,
    roundts=1,
)
# %% [markdown]
# These few links are everything that is needed to run a REMix model once the
# data is set up. For the evaluation of results head over to part c of this
# tutorial.

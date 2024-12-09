# %% [markdown]
# ## Part b: running the model
#
# For further explanations on what happens in this part of the model execution,
# refer to `tutorial_101b_run`.
# %%
# loading model built in `tutorial_103a_build`
from remix.framework import Instance
import pathlib as pt

_path_tut103_data = pt.Path("../tutorial_103/data")

if not _path_tut103_data.exists():
    raise IOError("You need to run tutorial 103a first!")

m = Instance.from_path(_path_tut103_data)

# %%
# running GAMS from Python script
m.run(
    resultfile="tutorial_103",
    lo=3,
    names=1,
    postcalc=1,
    roundts=1,
)

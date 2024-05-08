# %% [markdown]
# ## Part b: running the model
#
# For further explanations on what happens in this part of the model execution, refer to `tutorial_101b_run`.
# %%
# loading model built in `tutorial_201a_build`
from remix.framework.api.instance import Instance
import pathlib as pt

_path_tut201_data = pt.Path("../tutorial_201/data")

if not _path_tut201_data.exists():
    raise IOError("You need to run tutorial 201a first!")

m = Instance.from_path(_path_tut201_data)

# %%
# running GAMS from Python script
m.run(
    resultfile="tutorial_201",
    lo=3,
    names=1,
    postcalc=1,
    roundts=1,
)

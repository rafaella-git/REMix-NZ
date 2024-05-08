# %% [markdown]
# ## Part c: evaluation of results
#
# As in the other tutorials, we want to have a quick look at the optimization results. The first part will be an
# overview of installed capacities again, in the second part we will have a look at the influence that electric
# vehicles have on the generation profiles.

# %%
# importing dependencies
from remix.framework.tools.gdx import GDXEval
import pandas as pd
import matplotlib.pyplot as plt

result_dir = "./results"

# define often-used shortcut
idx = pd.IndexSlice
# %%
# read in the output `*.gdx` file from the optimization in GAMS
results = GDXEval(f"{result_dir}/tutorial_201.gdx")
# %% [markdown]
# ### Evaluating converter capacities
#
# In a first step, we have a look at the installed capacities of the converter technologies in the model.

# %%
# convert converter capacities to a Pandas DataFrame
caps = results["converter_caps"]

ev_uc_caps = caps.loc[idx[:, :, "EVs_UC", :, "total"], :]
ev_cc_caps = caps.loc[idx[:, :, "EVs_CC", :, "total"], :]

print(ev_cc_caps)
print(ev_uc_caps)
# %% [markdown]
# ### Evaluating EV demand profiles

# %%
# convert commodity balances to a Pandas DataFrame
technology_colors = {"EVs": "#FC9A99", "EVs_CC": "#FC9A99", "EVs_UC": "#9dc3e6"}
background_color = "#fffaebff"  # #FC9A99
commodities = results["commodity_balance"]
# random weekly timeslice
timeslice = range(8692 - 168, 8692)

# filter for EV storage
EV_storage = (
    commodities.loc[idx[:, "R1_model", "2030", ["EVs_CC", "EVs_UC"], "EV_Stored"], :]
    .dropna()
    .groupby(["timeModel", "techs"])
    .sum()
    .unstack("techs")
    .fillna(0)
)
EV_storage.columns = EV_storage.columns.get_level_values(1)
EV_storage = EV_storage.loc[:, (EV_storage != 0).any()]  # drop columns with only zeros


fig, ax1 = plt.subplots(figsize=(10, 6))
fig.patch.set_facecolor(background_color)
ax1.set_facecolor(background_color)
ax1.set_title("EVs profiles for CC and UC")
EV_storage.iloc[timeslice].plot.line(ax=ax1, color=technology_colors)
plt.legend(loc=(0.0, 1.05))
plt.ylabel("Demand in GWh_el")

fig.subplots_adjust(bottom=0.1 * EV_storage.index.nlevels)
# %%
# EV demand
EV_demand = (
    commodities.loc[idx[:, "R1_model", "2030", "EVs", "driving_energy"], :]
    .dropna()
    .groupby(["timeModel", "techs"])
    .sum()
    .unstack("techs")
    .fillna(0)
)
EV_demand.columns = EV_demand.columns.get_level_values(1)
EV_demand = EV_demand.loc[:, (EV_demand != 0).any()]  # drop columns with only zeros

fig, ax1 = plt.subplots(figsize=(10, 6))
fig.patch.set_facecolor(background_color)
ax1.set_facecolor(background_color)
ax1.set_title("EV Demand")
EV_demand.iloc[timeslice].plot.line(stacked=True, ax=ax1, color=technology_colors)
plt.legend(loc=(0.0, 1.05))
plt.ylabel("Total driving demand EVs")

fig.subplots_adjust(bottom=0.1 * EV_demand.index.nlevels)
plt.show()
# %% [markdown]
#
# This concludes the tutorial. You should now be familiar with implementing
# electric vehicles in a REMix model, both for a controlled (including V2G) and
# an uncontrolled charging scenario.

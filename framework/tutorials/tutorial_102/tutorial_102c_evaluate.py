# %% [markdown]
# ## Part c: evaluation of results
#
# As in tutorial_101, we want to have a quick look at the optimization results.
# The first part will be an overview of installed capacities again, in the
# second part we will have a look at the influence that the battery has on the
# generation profiles.

# %%
# importing dependencies
from remix.framework import GDXEval
import pandas as pd
import matplotlib.pyplot as plt

result_dir = "./results"

# define often-used shortcut
idx = pd.IndexSlice
# %%
# read in the output `*.gdx` file from the optimization in GAMS
results = GDXEval(f"{result_dir}/tutorial_102.gdx")
# %% [markdown]
# ### Evaluating converter capacities

# %%
# convert converter capacities to a Pandas DataFrame
caps = results["converter_caps"]

print(caps.loc[idx[:, :, :, "Elec", "total"], :])
# %% [markdown]
# ### Evaluating generation profiles with storage possibility
#
# In this tutorial, we are having a look at the generation profiles to evaluate
# the use of batteries in the model we have set up.
# We are using it for the weeks with the highest and lowest usage of batteries.
# Additionally, we are also plotting the respective storage levels during these
# weeks.

# %%
# convert commodity balances to a Pandas DataFrame
commodities = results["commodity_balance"]

# filter for all electricity balances
elec = (
    commodities.loc[idx[:, "R1_model", "2030", :, "Elec"], :]
    .dropna()
    .groupby(["timeModel", "techs"])
    .sum()
    .unstack("techs")
    .fillna(0)
)
elec.columns = elec.columns.get_level_values(1)
elec = elec.loc[:, (elec != 0).any()]  # drop columns with only zeros
# %%
# electricity demand
demand = elec.copy().loc[idx[:], "Demand"]

# generating electricity
positive = elec.copy().drop(columns=["Demand"])
positive = positive[positive > 0].fillna(0)

# consuming electricity
negative = elec.copy().drop(columns=["Demand"])
negative = negative[negative < 0].fillna(0)
# %%
# Finding out the week with the highest / lowest usage of lithium-ion batteries
# can be done by applying a rolling mean with an interval of 168 hours (the
# result of the calculation outputs the last hour of the interval).
hours_per_interval = 168  # 168 hours per week
rolling_mean = positive[["Battery"]].sum(axis=1).rolling(hours_per_interval).mean()
mean_max = rolling_mean.argmax()
mean_min = rolling_mean.argmin()

print(
    f"The week with the highest usage of lithium-ion batteries ends in hour {mean_max}."
)

print(f"The lowest usage is in the week before hour {mean_min}.")

# %%
# reading in data for plotting storage levels in model region
storlvl = results["storage_level_out"]
storlvl = (
    storlvl[storlvl.index.isin(["R1_model"], level=1)].unstack(3).droplevel([1, 2, 3])
)
storlvl.columns = storlvl.columns.get_level_values(1)
# %%
# week with the highest usage of batteries
technology_colors = {
    "CCGT": "#ff5f2d",
    "PV": "#ffc000",
    "WindOnshore": "#9dc3e6",
    "Battery": "#cc99ff",
}
background_color = "#fffaebff"
timeslice = range(mean_max - hours_per_interval, mean_max)

fig, ax1 = plt.subplots(figsize=(10, 6))
fig.patch.set_facecolor(background_color)
ax1.set_title("Week with the highest usage of batteries")
ax1.set_facecolor(background_color)
positive.iloc[timeslice].plot.area(stacked=True, ax=ax1, color=technology_colors)
plt.legend(loc=(0.0, 1.05))
plt.ylabel("Generation in GWh_el")

ax2 = ax1.twinx()
demand.iloc[timeslice].mul(-1).plot(kind="line", stacked=True, ax=ax2, color="black")
plt.legend(loc=(0.8, 1.05))
plt.ylabel("Demand in GWh_el")
ax2.set_ylim(ax1.get_ylim())

fig.subplots_adjust(bottom=0.1 * demand.index.nlevels)
# %%
# plotting the storage level during the week with the highest battery usage
fig, ax1 = plt.subplots(figsize=(10, 6))
fig.patch.set_facecolor(background_color)
ax1.set_title("Storage level during week with the highest usage of batteries")
ax1.set_facecolor(background_color)
storlvl.iloc[timeslice].plot.line(stacked=True, ax=ax1, color=technology_colors)
plt.legend(loc=(0.0, 1.05))
plt.ylabel("Storage level in GWh_el")

fig.subplots_adjust(bottom=0.1 * storlvl.index.nlevels)
# %%
# week with the lowest usage of batteries
timeslice = range(mean_min - hours_per_interval, mean_min)

fig, ax1 = plt.subplots(figsize=(10, 6))
fig.patch.set_facecolor(background_color)
ax1.set_title("Week with the lowest usage of batteries")
ax1.set_facecolor(background_color)
positive.iloc[timeslice].plot.area(stacked=True, ax=ax1, color=technology_colors)
plt.legend(loc=(0.0, 1.05))
plt.ylabel("Generation in GWh_el")

ax2 = ax1.twinx()
demand.iloc[timeslice].mul(-1).plot(kind="line", stacked=True, ax=ax2, color="black")
plt.legend(loc=(0.8, 1.05))
plt.ylabel("Demand in GWh_el")
ax2.set_ylim(ax1.get_ylim())

fig.subplots_adjust(bottom=0.1 * demand.index.nlevels)
# %%
# plotting the storage level during the week with the lowest battery usage
fig, ax1 = plt.subplots(figsize=(10, 6))
fig.patch.set_facecolor(background_color)
ax1.set_title("Storage level during week with the lowest usage of batteries")
ax1.set_facecolor(background_color)
storlvl.iloc[timeslice].plot.line(stacked=True, ax=ax1, color=technology_colors)
plt.legend(loc=(0.0, 1.05))
plt.ylabel("Storage level in GWh_el")

fig.subplots_adjust(bottom=0.1 * storlvl.index.nlevels)
plt.show()
# %% [markdown]
# This concludes this tutorial.
# To get some more insights into REMix, you might want to have a look at the
# bonus tasks.
# In them, some examples are given on how to extend the model.

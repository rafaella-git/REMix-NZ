# %% [markdown]
# ## Part c: evaluation of results
#
# We have now learned how to set up and run a basic energy system model in REMix
# that should have solved correctly.
#
# You can now open up the file `tutorial_101.gdx` in the `results` folder to see
# the results of the optimization.
# Try to find out the total system costs and how the system costs are split up
# between investment and fuel costs as well as how many gas turbines are
# required to fulfil the electricity demand at each hour of the year.
#
# Alternatively, you can use the code in the following sections for accessing
# the output files.
# Apart from examples on how to visualize chosen results with Python, some ideas
# on how to interpret them are given.
#
# First, we need to read in the results from the `*.gdx` file (after importing
# all necessary dependencies).

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

results = GDXEval(f"{result_dir}/tutorial_101.gdx")
# %% [markdown]
# ### Evaluating converter capacities
#
# In a first step, we have a look at the installed capacities of the converter
# technologies in the model.

# %%
# convert converter capacities to a Pandas DataFrame
caps = results["converter_caps"]

print(caps.loc[idx[:, :, :, "Elec", "total"], :])
# %% [markdown]
# ### Evaluating generation and demand profiles
#
# Similar to the capacities we can also check the hourly electricity generation
# and demand profiles.
# Here we load the commodity balance, filter for all electricity generation and
# consumption technologies and plot the generation for one week in spring.

# %%
# convert commodity balances to a Pandas DataFrame
commodities = results["commodity_balance"]

generation = (
    commodities[commodities > 0]
    .loc[idx[:, "R1_model", "2030", :, "Elec"], :]
    .dropna()
    .groupby(["timeModel", "techs"])
    .sum()
    .unstack("techs")
    .fillna(0)
)
generation.columns = generation.columns.get_level_values(1)

demand = (
    commodities[commodities < 0]
    .loc[idx[:, "R1_model", "2030", :, "Elec"], :]
    .dropna()
    .groupby(["timeModel", "techs"])
    .sum()
    .unstack("techs")
    .fillna(0)
)
demand.columns = demand.columns.get_level_values(1)
# %% [markdown]
# Pandas allows us to get a quick statistical look into datasets via the
# .describe() method.
# This is especially helpful for profile data.
#
# Try to find out the average utilization of the renewable energy technologies
# and compare their generation to the theoretically possible maximum feed-in.

# %%
print(generation.describe())
# %% [markdown]
# Another helpful tool for analyzing the system operation is running averages.
# Here we use Pandas to calculate the running average per week (168 hours) in
# order to find out the week with the highest and lowest PV feed-in and plot
# both to get a visual impression of the hourly system operation.

# %%
# Finding out the week with the highest / lowest feed-in can be done by applying
# a rolling mean with an interval of 168 hours (the result of the calculation
# outputs the last hour of the interval)
hours_per_interval = 168  # 168 hours per week
rolling_mean = generation["PV"].rolling(hours_per_interval).mean()
mean_max = rolling_mean.argmax()
mean_min = rolling_mean.argmin()

print(f"The week with the highest electricity feed-in ends in hour {mean_max}.")
print(f"The week with the lowest electricity feed-in ends in hour {mean_min}.")
# %%
# plotting week with highest electricity feed-in
technology_colors = {"CCGT": "#ff5f2d", "PV": "#ffc000", "WindOnshore": "#9dc3e6"}
background_color = "#fffaebff"
timeslice = range(mean_max - hours_per_interval, mean_max)

fig, ax1 = plt.subplots(figsize=(10, 6))
fig.patch.set_facecolor(background_color)
ax1.set_title("Week with the highest renewable feed-in")
ax1.set_facecolor(background_color)
generation.iloc[timeslice].plot.area(stacked=True, ax=ax1, color=technology_colors)
plt.legend(loc=(0.0, 1.05))
plt.ylabel("Generation in GWh_el")

ax2 = ax1.twinx()
demand.iloc[timeslice].mul(-1).plot(kind="line", stacked=True, ax=ax2, color="black")
plt.legend(loc=(0.8, 1.05))
plt.ylabel("Demand in GWh_el")
ax2.set_ylim(ax1.get_ylim())

fig.subplots_adjust(bottom=0.1 * demand.index.nlevels)
# %%
# plotting week with highest electricity feed-in
timeslice = range(mean_min - hours_per_interval, mean_min)

fig, ax1 = plt.subplots(figsize=(10, 6))
fig.patch.set_facecolor(background_color)
ax1.set_title("Week with the lowest renewable feed-in")
ax1.set_facecolor(background_color)
generation.iloc[timeslice].plot.area(stacked=True, ax=ax1, color=technology_colors)
plt.legend(loc=(0.0, 1.05))
plt.ylabel("Generation in GWh_el")

ax2 = ax1.twinx()
demand.iloc[timeslice].mul(-1).plot(kind="line", stacked=True, ax=ax2, color="black")
plt.legend(loc=(0.8, 1.05))
plt.ylabel("Demand in GWh_el")
ax2.set_ylim(ax1.get_ylim())

fig.subplots_adjust(bottom=0.1 * demand.index.nlevels)
plt.show()
# %% [markdown]
# This concludes the first basic REMix tutorial.
# To get some more insights into REMix, you might want to have a look at the
# bonus tasks.
# In there, an introduction on how to understand and handle error messages of
# GAMS is given.
# Also some examples are given on how to extend the model.

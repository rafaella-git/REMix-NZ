# %% [markdown]
# ## Part c: evaluation of results
#
# As in the first basic tutorial, we want to have a quick look at the
# optimization results.
# The first part will be an overview of installed capacities again, in the
# second part we will have a look at the influence that the battery has on the
# generation profiles.

# %%
# importing dependencies
from remix.framework import GDXEval
from remix.framework.tools.plots import plot_network, plot_choropleth
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

result_dir = "./results"

# define often-used shortcut
idx = pd.IndexSlice

# deciding on coordinate reference system to plot transfer on maps
crs = 3857
# %%
# read in the output `*.gdx` file from the optimization in GAMS
results = GDXEval(f"{result_dir}/tutorial_103.gdx")
# %% [markdown]
# ### Evaluating converter capacities

# %%
# convert converter capacities to a Pandas DataFrame
caps = results["converter_caps"]
caps = caps[caps > 0.01].dropna()  # Remove all capacities with less than 10 MW

caps.loc[idx[:, "2030", :, "Elec", "total"], :].round(2)
# %% [markdown]
# We can now check the installed connection capacities between the different
# model nodes.
# This is done in a similar fashion as with the converter capacities.
# %%
transfer_caps = results["transfer_caps"]

transfer_caps.loc[idx[:, :, :, "2030", :, "Elec", "total"], :].round(2)
# %% [markdown]
# To get a feeling on where we benefit from the electrical network, we can check
# the annual transferred energy between model nodes.
# Since the network contains information on the direction of flows, we need to
# also account for the direction.
# The energy flow from model region A (nodesModel) to B (nodesModel_a) is
# defined as positive, whereas the flow from B to A is accounted negative.
# With the "balanceType" entry we can check for the individual flows from A to
# B (positive), flows from B to A (negative), annual sum of directed flows
# (net = positive + negative), annual sum of energy transferred
# (gross = positive - negative), or link utilization
# (flh = gross / link capacity).

# %%
transfer_flows = results["transfer_flows_annual"]
transfer_flows = transfer_flows[
    transfer_flows.abs() > 0.1
].dropna()  # Remove all flows with less than 0.1 GWh

transfer_flows = (
    transfer_flows.loc[idx[:, :, :, :, :, "Elec", "net"], :].div(1e3).round(2)
)  # Convert to TWh

transfer_flows
# %% [markdown]
#  We identified `R3_model` as the main importing model region, so we can
# further check the behavior of the hourly electricity supply.

# %%
# visualization of commodity balance for model region R3_model
commodities = results["commodity_balance"]

elec_R3 = (
    commodities.loc[idx[:, "R3_model", "2030", :, "Elec"], :]
    .groupby(["timeModel", "techs"])
    .sum()
    .unstack("techs")
)

demand_R3 = elec_R3.loc[:, idx[:, "Demand"]]
demand_R3.columns = demand_R3.columns.get_level_values(1)

positive_R3 = elec_R3.drop(columns=("value", "Demand"))
positive_R3 = positive_R3[positive_R3 > 0].fillna(0)
positive_R3.columns = positive_R3.columns.get_level_values(1)

negative_R3 = elec_R3.drop(columns=("value", "Demand"))
negative_R3 = negative_R3[negative_R3 < 0].fillna(0)
negative_R3.columns = negative_R3.columns.get_level_values(1)
# %% [markdown]
# Finding out the week with the highest / lowest usage of HVDC
# (The result of the calculation outputs the last hour of the interval)

# %%
hours_per_interval = 168  # 168 hours per week
rolling_mean = positive_R3[["HVDC"]].sum(axis=1).rolling(hours_per_interval).mean()
last_hour_max = rolling_mean.argmax()
last_hour_min = rolling_mean.argmin()
# %%
# visualization of generation in week with highest transmission to model region R3
technology_colors = {
    "CCGT": "#ff5f2d",
    "PV": "#ffc000",
    "WindOnshore": "#9dc3e6",
    "HVDC": "#70ad47",
}
background_color = "#fffaebff"
timeslice = range(last_hour_max - hours_per_interval, last_hour_max)

fig, ax1 = plt.subplots(figsize=(10, 6))
fig.patch.set_facecolor(background_color)
ax1.set_facecolor(background_color)
positive_R3.iloc[timeslice].plot.area(stacked=True, ax=ax1, color=technology_colors)

plt.legend(loc="upper left")
plt.ylabel("Generation in GWh_el")
plt.title("Generation in week with highest transmission to model region R3")

ax2 = ax1.twinx()
demand_R3.iloc[timeslice].mul(-1).plot(kind="line", stacked=True, ax=ax2, color="black")
plt.legend(loc="upper right")
plt.ylabel("Demand in GWh_el")
ax2.set_ylim(ax1.get_ylim())

fig.subplots_adjust(bottom=0.1 * demand_R3.index.nlevels)
# %%
# map defaults and options
shp_file = "../_input/tutorial.geojson"
shp_attrcol = "ISO"

lat = [35.0, 57.0]
lon = [-6.0, 20.0]

centroids = {
    "FR": (2.34, 47.15),
    "DE": (10.18, 51.35),
    "IT": (12.80, 42.61),
    "CH": (8.17, 46.86),
}

plt.rcParams.update({"figure.autolayout": True})  # use tight_layout
plt.rcParams.update({"figure.titlesize": 20})  # size subtitle
plt.rcParams.update({"figure.dpi": 75})  # size subtitle
plt.rcParams.update({"savefig.dpi": 300})  # dpi
plt.rcParams.update({"font.size": 16})  # default font size
plt.rcParams.update({"axes.titlesize": 20})  # title size per subplot
plt.rcParams.update({"axes.labelsize": 18})  # label size colormap

figsize_dual = (13.0, 6)
figsize_single = (7.5, 6)
# %%
# visualization of renewable generation per model region
df_annual = results["commodity_balance_annual"]

df_annual.index = df_annual.index.set_levels(
    df_annual.index.levels[0].rename_categories(
        {"R3_model": "CH", "R1_model": "DE", "R2_model": "FR", "R4_model": "IT"}
    ),
    level="accNodesModel",
)
df_annual = df_annual.loc[idx[["DE", "CH", "FR", "IT"], "2030", :, "Elec", "net"]]

map_pv = (
    df_annual.loc[idx[:, :, "PV"], idx[:]]
    .groupby("accNodesModel", observed=True)
    .sum()
    .div(1e3)
)
map_wind = (
    df_annual.loc[idx[:, :, "WindOnshore"], idx[:]]
    .groupby("accNodesModel", observed=True)
    .sum()
    .div(1e3)
)

fig = plt.figure(figsize=figsize_dual)
fig.patch.set_facecolor(background_color)
ax1 = fig.add_subplot(121, projection=ccrs.epsg(crs))
ax1.set_facecolor("#F0F8FF")
plot_choropleth(
    map_pv,
    shp_file,
    shp_attrcol,
    fig=fig,
    ax=ax1,
    lat=lat,
    lon=lon,
    title="Power generation from PV",
    clabel="Energy in TWh",
    cmap="Oranges",
)

ax2 = fig.add_subplot(122, projection=ccrs.epsg(crs))
ax2.set_facecolor("#F0F8FF")
plot_choropleth(
    map_wind,
    shp_file,
    shp_attrcol,
    lat=lat,
    lon=lon,
    fig=fig,
    ax=ax2,
    title="Power generation from wind",
    clabel="Energy in TWh",
    cmap="Blues",
)

# %%
# visualization of energy flows between model regions
n2n_flow = results["transfer_flows_annual"]
for level in ["nodesModel_start", "nodesModel_end"]:
    n2n_flow.index = n2n_flow.index.set_levels(
        n2n_flow.index.levels[0].rename_categories(
            {"R3_model": "CH", "R1_model": "DE", "R2_model": "FR", "R4_model": "IT"}
        ),
        level=level,
    )

n2n_flow = n2n_flow[n2n_flow != 0].dropna(how="all")

flow_elec = (
    n2n_flow[n2n_flow != 0]
    .loc[idx[:, :, :, :, :, "Elec", "net"]]
    .dropna()
    .groupby(["nodesModel_start", "nodesModel_end"], observed=True)
    .sum()
    .div(1e3)
    .abs()
)

flow_elec = flow_elec[flow_elec > 0.1]
flow_flh = (
    n2n_flow[n2n_flow != 0]
    .loc[idx[:, :, :, :, :, "Elec", "flh"]]
    .dropna()
    .groupby(["nodesModel_start", "nodesModel_end"], observed=True)
    .sum()
    .abs()
)

flow_flh = flow_flh[flow_flh > 10]

fig = plt.figure(figsize=figsize_dual)
fig.patch.set_facecolor(background_color)
ax1 = fig.add_subplot(121, projection=ccrs.epsg(crs))
ax1.set_facecolor("#F0F8FF")
plot_network(
    data=flow_elec,
    backgr_shp=shp_file,
    backgr_attr=shp_attrcol,
    lat=lat,
    lon=lon,
    centroids=centroids,
    title="Annual net transmission",
    clabel="Energy in TWh",
    cmap="Greens",
    fig=fig,
    ax=ax1,
    cbar=True,
)

ax2 = fig.add_subplot(122, projection=ccrs.epsg(crs))
ax2.set_facecolor("#F0F8FF")
plot_network(
    data=flow_flh,
    backgr_shp=shp_file,
    backgr_attr=shp_attrcol,
    lat=lat,
    lon=lon,
    centroids=centroids,
    title="Annual full load hours",
    clabel="Full load hours",
    cmap="Reds",
    fig=fig,
    ax=ax2,
    cbar=True,
    limits=[0, 8760],
)

plt.show()
# %% [markdown]
#
# This concludes the tutorial. You should now have a fundamental understanding
# on how to use sources, sinks, converters, storage technologies and links in
# REMix. These build the basis for working with energy system models.

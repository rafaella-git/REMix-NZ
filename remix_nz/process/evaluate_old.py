# %% [markdown]
# ## Part c: evaluation of results
# https://dlr-ve.gitlab.io/esy/remix/framework/dev/getting-started/tutorials/tutorial-103.html

technology_colors = {
    "CCGT": "#ff5f2d",
    "PV": "#ffc000",
    "pv_central_fixed": "#ffc000",
    "wind_onshore": "#9dc3e6",
    "HVDC": "#70ad47",
}
background_color = "#fffaebff"


# %% [markdown]
# Overview of installed capacities

# %%
# importing dependencies
from remix.framework.tools.gdx import GDXEval
from remix.framework.tools.plots import plot_network, plot_choropleth
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd



#%%
# ## Define useful functions
# ### 1. Define datafame from gdx results file | (folder of file, name of dataframe)
# ### 2. Plot capacities | (dataframe caps, tech you wnat to plot, year)
# ### 3. Plot transfer (like triangles) | dataframe transfer flows/balance?
# ### 4. Plot ISO plots x=day y=hour heat=intensity | dataframe balance
# ### 5. Define function to export a figure | (figure, name, destiny folder) 




# %%
# Define case to be evaluated (based on demand file name and the yeas optimised)
indx=0
# Demand files available (different scenarios)
files_lst=["nz_profile_11nodes","medpop_evs_base","low_pop_out_base","med_pop_out_base","high_pop_out_base"] #493 as in the course 493, this is for Liv and Sam
yrs=[2020,2030,2040,2050] # years optimised
yrs_str='-'.join([str(item) for item in yrs])
demand_file=files_lst[indx] 
case_name=f"{demand_file}_{yrs_str}"
idx = pd.IndexSlice #often used shortcut


path_base = "C:/Local/REMix"  
path_input = f"{path_base}/remix_nz/input"
path_output = f"{path_base}/remix_nz/output" 
path_result = f"{path_output}/{case_name}/result" 
path_geo = f"{path_input}/shapefiles"      # geojson
geofile="11regionsNZ.geojson"


# %%
# read in the output `*.gdx` file from the optimization in GAMS
results = GDXEval(f"{path_result}/{case_name}.gdx")


# %% [markdown]
# ### Evaluating converter (PV, wind, battery) capacities
caps = results["converter_caps"]   # convert converter capacities to a Pandas DataFrame
caps = caps[caps > 0.01].dropna()  # remove all capacities with less than 10 MW
caps.loc[idx[:, "2030", :, "Elec", "total"], :].round(2) # accNodesModel=all, accYears=2030 only, techs=all, commodit=only Elec, capType=total


# %% [markdown]
# ### Evaluating installed connection capacities between the different model nodes.
transfer_caps = results["transfer_caps"]
# nodesModel_start=all, nodesModel_end=all, linksModel=all, years=2030 only, techs=all (we only have HV), commodity=Elec, capType=total
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
# (netto = positive + negative), annual sum of energy transferred
# (brutto = positive - negative), or link utilization
# (flh = brutto / link capacity).
transfer_flows = results["transfer_flows_annual"]
transfer_flows = transfer_flows[transfer_flows.abs() > 0.1].dropna()  # Remove all flows with less than 0.1 GWh
transfer_flows = (transfer_flows.loc[idx[:, :, :, :, :, "Elec", "netto"], :].div(1e3).round(2))  # Convert to TWh
transfer_flows
# %% [markdown]
#  We identified `CEN` as the main importing model region, so we can
# further check the behavior of the hourly electricity supply (if we wanted to).


# %% [markdown]
#  ## Maps
nodes_lst=["NIS","AKL","WTO","TRN","BOP","HBY","CEN","WEL","NEL","CAN","OTG"]
nodes_lat=[-35.8758611, -36.9547333, -38.4192333, -39.3342222, -37.9867861, -39.5516472, -40.2813000, -41.1502722, -41.6735722, -43.8619833, -45.481222]
nodes_lon=[174.4669472, 174.8625250, 175.8000111, 174.3204333, 176.8294056, 176.8208167, 175.6404750, 174.9811500, 172.8737917, 171.3427694, 169.3195194]    
nodes_coords=[(nodes_lat[i],nodes_lon[i]) for i in range(0,len(nodes_lon))]
# dictionary with 'Node': [lat,long]
coord_dict = dict.fromkeys(nodes_lst)
for key, value in zip(coord_dict.keys(), nodes_coords):
    coord_dict[key] = value
print(coord_dict)



plt.rcParams.update({"figure.autolayout": True})  # use tight_layout
plt.rcParams.update({"figure.titlesize": 20})  # size subtitle
plt.rcParams.update({"figure.dpi": 75})  # size subtitle
plt.rcParams.update({"savefig.dpi": 300})  # dpi
plt.rcParams.update({"font.size": 16})  # default font size
plt.rcParams.update({"axes.titlesize": 20})  # title size per subplot
plt.rcParams.update({"axes.labelsize": 18})  # label size colormap

figsize_dual = (13.0, 6)
figsize_single = (7.5, 6)

#limits of the plot 
lat = [-47.758239, -34.111702]
lon = [165.692103, 179.050919]


centroids = coord_dict



shp_file = f"{path_geo}/11regionsNZ"
shp_attrcol = "id"

# visualization of renewable generation per model region
df_annual = results["commodity_balance_annual"]
df_annual = df_annual.loc[idx[nodes_lst, "2030", :, "Elec", "netto"]]


map_pv = df_annual.loc[idx[:, :, "pv_central_fixed"], idx[:]].groupby("accNodesModel").sum().div(1e3)
map_wind = (
    df_annual.loc[idx[:, :, "wind_onshore"], idx[:]]
    .groupby("accNodesModel")
    .sum()
    .div(1e3)
)



fig = plt.figure(figsize=figsize_dual)
fig.patch.set_facecolor(background_color)
ax1 = fig.add_subplot(121)
ax1.set_facecolor("#F0F8FF")
plot_choropleth(
    map_pv,   #df
    shp_file,
    shp_attrcol,
    lat,
    lon,
    title="Power generation from PV",
    clabel="Energy in TWh",
    cmap="Oranges",
    ax=ax1,
)

ax2 = fig.add_subplot(122)
ax2.set_facecolor("#F0F8FF")
plot_choropleth(
    map_wind,
    shp_file,
    shp_attrcol,
    lat,
    lon,
    title="Power generation from wind",
    clabel="Energy in TWh",
    cmap="Blues",
    ax=ax2,
)



# %% [markdown]
# ## Part c: evaluation of results
# https://dlr-ve.gitlab.io/esy/remix/framework/dev/getting-started/tutorials/tutorial-103.html

# technology_colors = {
#     "CCGT": "#ff5f2d",
#     'PV open area (fixed)': "#ffc000",
#     'Offshore wind (floating)': "#ffc000",
#     'CCGT (H2)': "#9dc3e6",
#     "HVDC": "#70ad47",
#     'Hydro': "#ffc000",
#     'H2 Storage': "#ffc000",
#     'Offshore wind (foundation)': "#ffc000",
#     'Fuel Cell (H2)': "#ffc000",
#     'PV open area (tracking)': "#ffc000",
#     'GT': "#ffc000",
#     'Electrolyser': "#ffc000",
#     'Coal': "#ffc000",
#     'Onshore wind': "#ffc000",
#     'Battery': "#ffc000",
#     'Geothermal': "#ffc000",
#     'OCGT': "#ffc000",
#     'Diesel': "#ffc000",
#     'PV decentralised': "#ffc000",
#     'Biomass': "#ffc000",
# }

# technology_colors = {
#     "Hydro": "#1F78B4",
#     "Geothermal": "#9A6436",
#     "PV open area (tracking)": "#FFD700",
#     "PV open area (fixed)": "#FFD700",
#     "PV decentralised": "#FFD700",
#     "PV": "#FFD700",
#     "Offshore wind": "#018B8C",
#     "Offshore wind (floating)": "#018B8C",
#     "Offshore wind (foundation)": "#018B8C",
#     "Onshore wind": "#20B1AA",
#     "Battery": "#B597CC",
#     "H2 Storage": "#C263CB",
#     "Fuel Cell (H2)": "#693D9A",
#     "Electrolyser": "#AF9968",
#     "CCGT (H2)": "#FCAE99",
#     "Modified Gas Turbine (H2)": "#FCAE99",
#     "CCGT": "#FF7F00",
#     "OCGT": "#FF7F00",
#     "GT": "#FF7F00",
#     "Gas Turbine": "#FF7F00",
#     "Diesel": "#5F5F5F",
#     "Coal": "#393839",
#     "HVDC": "#000000",
#     "HVAC": "#000000",
#     "Biomass": "#779241",
# }


#from pypsa
technology_colors = {
    "Hydro": '#298c81',
    "Geothermal":  "#ba91b1",
    "PV open area (tracking)": "#f9d002",
    "PV open area (fixed)": "#f9d002",
    "PV decentralised": "#f9d002",
    "PV": "#f9d002",
    "Offshore wind": "#6895dd",
    "Offshore wind (floating)": "#6895dd",
    "Offshore wind (foundation)": "#6895dd",
    "Onshore wind": "#4F7EC9", #: "#235ebc"
    "Battery": "#708090",
    "H2 Storage": '#bf13a0',
    "H2_CCGT": '#bf13a0',
    "Fuel Cell (H2)": '#c251ae',
    "Electrolyser": "#AF9968",
    "CCGT (H2)": "#991f83",
    "Modified Gas Turbine (H2)": "#991f83",
    "CCGT": "#ee8340",
    "OCGT": "#FAA460",
    "GT": "#ee8340",
    "Gas Turbine": "#FAA460",
    "Diesel": '#B5A642',
    "Coal":  "#505050",
    "HVDC": "#000000",
    "HVAC": "#000000",
    "Biomass": "#008000",
}


background_color = "#FFFFFFFF"#fffaebff"

# Import dependencies
import os
import time
import json
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from remix.framework.tools.gdx import GDXEval
from remix.framework.tools.plots import plot_network, plot_choropleth
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd



will_elec = ["00-test-elec","01-battery-distributed", "02-battery-overnight", "03-battery-recharging", "04-battery-solar"]
will_h2 = ["01-h2-distributed", "02-h2-overnight", "03-h2-recharging", "04-h2-solar"]
sdewes_ap = ["base", "high", "med", "ev-med", "low"]
mbie=["base","h2pos", "h2"]
europe=["h2-lut-domestic", "h2-lut-exports", "h2-pypsa","h2-pypsa-exports-domestic", "h2-pypsa-exports-20","h2-pypsa-exports-40","h2-pypsa-exports-200"]

  
scenario_dict = {       
    "will": [will_h2, [2020, 2035, 2050]],
    "sdewes-ap": [sdewes_ap, [2020, 2030, 2040, 2050]],
    "mbie": [mbie, [2020, 2030, 2040, 2050]],
    "europe": [europe, [2020, 2030,2050]]
}
group_name="europe"
files_lst = scenario_dict[group_name][0]
yrs_sel = scenario_dict[group_name][1] # [2020, 2025, 2030, 2035, 2040, 2045, 2050]
yrs_str='-'.join([str(item) for item in yrs_sel])
yrs_to_calc = [2020, 2025, 2030, 2035, 2040, 2045, 2050]
indx=0
yrs_string = [str(x) for x in yrs_sel]

# Define paths/directories
path_base = "C:/Local/REMix"
path_input = f"{path_base}/remix_nz/input"
path_demand = f"{path_input}/demand/{group_name}"
path_profiles = f"{path_input}/profiles"      # renewables
path_brownfield = f"{path_input}/brownfield"  # info hydro and existing power plants database
path_geo = f"{path_input}/shapefiles"         # geojson
geofile="11regionsNZ.geojson"
demand_file=files_lst[indx] 
case_name=f"{demand_file}_{yrs_str}"
path_output = f"{path_base}/remix_nz/output/{group_name}"
data_dir = Path(f"{path_output}/{case_name}/data")
data_dir.mkdir(parents=True, exist_ok=True)
results_dir = Path(f"{path_output}/{case_name}/result")
results_dir.mkdir(parents=True, exist_ok=True)



# Define paths/directories
path_base = "C:/Local/REMix"
path_input = f"{path_base}/remix_nz/input"
path_demand = f"{path_input}/demand/{group_name}"
path_profiles = f"{path_input}/profiles"      # renewables
path_brownfield = f"{path_input}/brownfield"  # info hydro and existing power plants database
path_geo = f"{path_input}/shapefiles"         # geojson
geofile="11regionsNZ.geojson"
demand_file=files_lst[indx] 
case_name=f"{demand_file}_{yrs_str}"
path_output = f"{path_base}/remix_nz/output/{group_name}"
data_dir = Path(f"{path_output}/{case_name}/data")
data_dir.mkdir(parents=True, exist_ok=True)
results_dir = Path(f"{path_output}/{case_name}/result")
results_dir.mkdir(parents=True, exist_ok=True)

# ## Define useful functions
# ### 1. Define datafame from gdx results file | (folder of file, name of dataframe)
results = GDXEval(f"{results_dir}/{case_name}.gdx")
# rename_techs = {
#     "H2_CCGT": "CCGT (H2)",
#     "pv_central_fixed": "PV open area (fixed)",
#     "pv_central_track_azimuth": "PV open area (tracking)",
#     "wind_onshore": "Onshore wind",
#     "wind_offshore_floating": "Offshore wind (floating)",
#     "wind_offshore_foundation": "Offshore wind (foundation)",
#     'geoth': "Geothermal",
#     'COAL': "Coal",
#     'DIE': "Diesel",
#     'H2_storage': "H2 Storage",
#     'pv_decentral': "PV decentralised",
#     'BIO': "Biomass",
#     "H2_FC": "Fuel Cell (H2)",
# }
rename_techs = {
    "H2_CCGT": "Modified Gas Turbine (H2)",
    "pv_central_fixed": "PV",
    "pv_central_track_azimuth": "PV",
    "wind_onshore": "Onshore wind",
    "wind_offshore_floating": "Offshore wind",
    "wind_offshore_foundation": "Offshore wind",
    'geoth': "Geothermal",
    'COAL': "Coal",
    'DIE': "Diesel",
    'H2_storage': "H2 Storage",
    'pv_decentral': "PV",
    'BIO': "Biomass",
    "H2_FC": "Fuel Cell (H2)",
    "CCGT":"Gas Turbine",
    "OCGT":"Gas Turbine",
    "GT": "Gas Turbine",
}


converter_caps = results['converter_caps'].rename(index=rename_techs)
converter_caps = converter_caps[converter_caps > 0.01].dropna()  # remove all capacities with less than 10 MW
# converter_caps.loc[idx[:, f"{year}", f"{tech}", "Elec", "total"], :].round(2) # accNodesModel, accYears, techs, commodity, capType
#print(list(set(converter_caps.index.get_level_values('techs'))))

commodity_balance_annual = results['commodity_balance_annual'].rename(index=rename_techs)
commodity_balance_annual = commodity_balance_annual[(commodity_balance_annual.abs() > 0.01)].dropna() # remove all values smaller than abs(10 MW)
# commodity_balance_annual.loc[idx['global', f"{year}", f"{tech}", "Elec", "net"], :].round(2)  # accNodesModel, accYears, techs, commodity, balanceType

storage_flows = results['storage_flows'].rename(index=rename_techs)
# storage_flows.loc[idx[:, "global", f"{year}", f"{tech}", f"{commodity}", ], :].round(2) # timeModel, accNodesModel, accYears, techs, commodity



# ### 2. Plot capacities | (dataframe caps, tech you wnat to plot, year)
# ### 3. Plot transfer (like triangles) | dataframe transfer flows/balance?
# ### 4. Plot ISO plots x=day y=hour heat=intensity | dataframe balance
# ### 5. Define function to export a figure | (figure, name, destiny folder) 




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
shp_file = f"{path_geo}/11regionsNZ.shp"
shp_attrcol = "id"
idx = pd.IndexSlice #often used shortcut



def carnot_isopleths_gen(res):
    map_dict = {
                # "Electrolyzer": {"title":"Electrolyzer output in GW", "techs":"Electrolyzer", "commodity":"H2", "scaling":1, "cmap": "viridis"},
    #             "Methaniser": {"title":"Methanizer output in GW", "techs":"Methaniser", "commodity":"CH4", "scaling":1, "cmap": "viridis"},
    
                "Hydro reservoirs": {"title":"Hydropower reservoir output in GW", "techs":"Hydro_reservoir", "commodity":"Water_in", "scaling":1, "cmap": "viridis", "text": " "},               
                "PV": {"title":"PV generation in GW", "techs":"PV", "commodity":"Elec", "scaling":1, "cmap": "viridis", "text": " "},
                #"H2 Storage": {"title":"H2 Storage operation in GW", "techs":"H2 Storage", "commodity":"elec", "scaling":1, "cmap": "viridis", "text": "b)"}, 
                "Batteries": {"title":"Battery operation in GW", "techs":"Battery", "commodity":"Elec", "scaling":1, "cmap": "viridis", "text": " "}, 
                #"Hydro _reservois": {"title":"Hydropower operation in GW", "techs":"Hydro_reservoir", "commodity":"Elec", "scaling":1, "cmap": "viridis", "text": " "}, 
                # "Wind": {"title":"Wind generation in GW", "techs":["WindOnshore","WindOffshore"], "commodity":"Elec", "scaling":1, "cmap": "viridis", "text": "a)"},
                # "CSP": {"title":"CSP generation in GW", "techs":"CSP_Powerblock", "commodity":"Elec", "scaling":1, "cmap": "viridis"},
    #            "RES": {"title":"Renewable generation in GW", "techs":["Photovoltaic","WindOnshore","WindOffshore","CSP_Powerblock"], "commodity":"Elec", "scaling":1, "cmap": "viridis"},
    #            "CHP": {"title":"CHP generation in GW", "techs":['TH_ExCCGT_NGas_XL','DH_Engine_NGas_M','DH_Engine_NGas_M', 'DH_FuelCell_H2_M'], "commodity":"Elec", "scaling":1, "cmap": "viridis"},
    #            "GT": {"title":"Reserve generation in GW", "techs":["CCGT", "CCGT_H2", "GT", "GT_H2"], "commodity":"Elec", "scaling":1, "cmap": "viridis"},
            }

    cb = res["commodity_balance"]
    cb = cb.rename(index=rename_techs)

    for i, j in map_dict.items():
        data = cb.loc[idx[:,:,"2050",j["techs"],j["commodity"]], idx[:]]

        data = data.groupby("timeModel").sum().mul(j["scaling"]).values.reshape(-1, 24).T

        fig = plt.figure(figsize=(12.5,4), layout="constrained")
        ax1 = fig.add_subplot(111)

        plt.imshow(data, aspect='auto', origin='upper', cmap=j["cmap"], interpolation='nearest', vmin=0,  vmax=np.max(data))
        cbar = plt.colorbar(pad=.02)
        if "text" in j.keys():
            ax1.text(-0.09, 0.96, j["text"], fontsize=16, transform=ax1.transAxes)

        plt.grid(False)

        # setting ticks positions
        t = np.arange(-0.5, 364.6, 30)
        ax1.xaxis.set_ticks(t)
        ax1.set_xticklabels(((t + 0.5)).astype(int))

        t = np.arange(-0.5, 23.6, 6)
        ax1.yaxis.set_ticks(t)
        ax1.set_yticklabels(list(map(lambda x: x + ":00", (t + 0.5).astype(int).astype(str))))

        cbar.set_label(j["title"], labelpad=15) # configure color bar

        filename = "figures/isopleths/iso_{}".format(i)
        Path(f"{results_dir}/{filename}").parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(f"{results_dir}/{filename}.png" , dpi=300)
        plt.savefig(f"{results_dir}/{filename}.svg")

        plt.show()

    

def plot_year_generation(df, ax, nodesData, years, techs, commodity,  upper=600):
    #commodity_balance_annual.loc[idx["global", yrs_string, 'PV open area (fixed)', "Elec", "net"], idx[:]]
    generation_capacities = df.loc[idx[nodesData, years, techs, commodity, "net"], idx[:]]
    generation_capacities = (
        generation_capacities.reset_index()
        .astype({col: object for col in generation_capacities.index.names})
        .groupby(["accNodesModel", "accYears", "techs"])
        .agg({"value": "sum"})
    )
    if len(generation_capacities.index.levels[0]) == 1:
        generation_capacities = generation_capacities.groupby(
            ["accYears", "techs"]
        ).agg({"value": "sum"})
    generation_capacities = generation_capacities.unstack() / 1000
    generation_capacities.columns = generation_capacities.columns.get_level_values(1)

    generation_capacities.plot.bar(stacked=True, ax=ax, color=technology_colors)
    # plt.xticks(rotation=0)
    plt.grid(axis='y', alpha=0.5)
    plt.ylim(bottom=0, top=upper)




# visualization of renewable generation per model region
df_annual = commodity_balance_annual
df_annual = df_annual.loc[idx[nodes_lst, "2030", :, "Elec", "net"]]

##['CCGT', 'PV open area (fixed)', 'Offshore wind (floating)', 'CCGT (H2)', 'Hydro', 'H2 Storage', 'Offshore wind (foundation)', 'H2_FC', 'PV open area (tracking)', 'GT', 'Electrolyser', 'Coal', 'Onshore wind', 'Battery', 'Geothermal', 'OCGT', 'Diesel', 'PV decentralised']

map_pv = df_annual.loc[idx[:, :, 'PV', :, :], idx[:]].groupby("accNodesModel").sum().div(1e3)
map_wind = (
    df_annual.loc[idx[:, :, 'Onshore wind',:, :], idx[:]]
    .groupby("accNodesModel")
    .sum()
    .div(1e3)
)

#for gen plots

generation_techs = ['Gas Turbine', 'PV', 'Offshore wind', 'Modified Gas Turbine (H2)', 'Hydro', 'H2 Storage', 'Fuel Cell (H2)', 'Electrolyser', 'Coal', 'Onshore wind', 'Battery', 'Geothermal', 'Diesel']
generation_techs_without_coal = [tech for tech in generation_techs if tech != 'Coal']

fig, ax1 = plt.subplots(figsize=(10, 6))
fig.patch.set_facecolor(background_color)
#ax1.set_facecolor(background_color)
plot_year_generation(
    commodity_balance_annual,
    ax1,
    "global",
    yrs_string,
    generation_techs_without_coal,
    "Elec",
    upper=370,
)
plt.grid(axis='y', alpha=0.5)
plt.legend(bbox_to_anchor=(1.0, 1.0))
plt.xlabel("Year")
plt.ylabel("Total annual generation in TWh")
plt.title("Electricity generation in New Zealand")

# fig = plt.figure(figsize=figsize_dual)
# fig.patch.set_facecolor(background_color)
# ax1 = fig.add_subplot(121)
# ax = fig.add_axes([0.02, 0.02, 0.96, 0.96], projection=ccrs.epsg(crs))
# ax1.set_facecolor("#F0F8FF")

"""
    data,
    shp,
    shp_attr,
    fig=None,
    ax=None,
    title: Optional[str] = None,
    background: bool = True,
    lat=None,
    lon=None,
    limits: Optional[list[int | float]] = None,
    cmap: str = "viridis",
    cbar: bool = True,
    clabel: Optional[str] = None,
    patch_coll=None,
    crs: int = 3857,
"""
# plot_choropleth(
#     data=map_pv,   #df
#     shp=shp_file,
#     shp_attr=shp_attrcol,
#     lat=lat,
#     lon=lon,
#     title="Power generation from PV",
#     clabel="Energy in TWh",
#     cmap="Oranges",
# )

# # ax2 = fig.add_subplot(122)
# # ax2.set_facecolor("#F0F8FF")
# plot_choropleth(
#     data=map_wind,
#     shp=shp_file,
#     shp_attr=shp_attrcol,
#     lat=lat,
#     lon=lon,
#     title="Power generation from wind",
#     clabel="Energy in TWh",
#     cmap="Blues",
# )

# plot_network(
#     data=map_pv,
#     shp=shp_file,
#     shp_attr=shp_attrcol,
#     centroids=centroids,
#     lat=lat,
#     lon=lon,
#     cmap = "viridis",
#     cbar=True,
#     clabel=None,
#     patch_coll=None,
#     crs = 3857,
# )


#carnot_isopleths_gen(results)


# example

#generation_capacities =commodity_balance_annual.loc[idx["global", yrs_string, 'PV open area (fixed)', "Elec", "net"], idx[:]]
#generation_techs = ['CCGT', 'PV open area (fixed)', 'Offshore wind (floating)', 'CCGT (H2)', 'Hydro', 'H2 Storage', 'Offshore wind (foundation)', 'Fuel Cell (H2)', 'PV open area (tracking)', 'GT', 'Electrolyser', 'Coal', 'Onshore wind', 'Battery', 'Geothermal', 'OCGT', 'Diesel', 'PV decentralised']




plt.show()



commodities = results["commodity_balance"]
commodities = results["commodity_balance"].rename(index=rename_techs)

generation = (
    commodities[commodities > 0]
    .loc[idx[:, "global", "2050", :, "Elec"], :]
    .dropna()
    .groupby(["timeModel", "techs"])
    .sum()
    .unstack("techs")
    .fillna(0)
)
generation.columns = generation.columns.get_level_values(1)

demand = (
    commodities[commodities < 0]
    .loc[idx[:, "global", "2050", :, "Elec"], :]
    .dropna()
    .groupby(["timeModel", "techs"])
    .sum()
    .unstack("techs")
    .fillna(0)
)


# Assuming 'generation' and 'demand' dataframes have an integer index representing hours
hours_per_week = 168
hours_per_two_weeks = hours_per_week * 2

# Define technology colors and background
#technology_colors = {"CCGT": "#ff5f2d", "PV": "#ffc000", "WindOnshore": "#9dc3e6"}
background_color = "#fffaebff"

# Choose hour-based midpoints for winter and summer periods
winter_mid_hour = 500  # Example: hour 500 for mid-winter
summer_mid_hour = 4500  # Example: hour 4500 for mid-summer

# Define two-week hour ranges around these midpoints
winter_timeslice = range(winter_mid_hour - hours_per_week, winter_mid_hour + hours_per_week)
summer_timeslice = range(summer_mid_hour - hours_per_week, summer_mid_hour + hours_per_week)

# Function to plot the time slice
def plot_two_weeks(timeslice, title):
    fig, ax1 = plt.subplots(figsize=(10, 6))
    fig.patch.set_facecolor(background_color)
    ax1.set_title(title)
    ax1.set_facecolor(background_color)


    for tech in generation.columns:
        if tech not in technology_colors:
            technology_colors[tech] = "#cccccc"  # Gray as a default color

    # Plot generation data
    generation.iloc[timeslice].plot.area(stacked=True, ax=ax1, color=technology_colors)
    plt.legend(loc=(0.0, 1.05))
    plt.ylabel("Generation in GWh_el")

    # Plot demand data on secondary y-axis
    ax2 = ax1.twinx()
    demand.iloc[timeslice].mul(-1).plot(kind="line", stacked=True, ax=ax2, color="black")
    plt.legend(loc=(0.8, 1.05))
    plt.ylabel("Demand in GWh_el")
    ax2.set_ylim(ax1.get_ylim())

    fig.subplots_adjust(bottom=0.1 * demand.index.nlevels)
    plt.show()

# Plotting winter and summer two-week periods
plot_two_weeks(winter_timeslice, "Two Weeks of Winter (Centered Around Hour 500)")
plot_two_weeks(summer_timeslice, "Two Weeks of Summer (Centered Around Hour 4500)")


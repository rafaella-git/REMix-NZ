# %% [markdown]
# import dependencies
import os
import time
import yaml
import fiona
import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path
import matplotlib
from remix.framework.tools.gdx import GDXEval
from remix.framework.tools.plots import plot_network, plot_choropleth
import geopandas as gpd
import os 
import matplotlib.pyplot as plt
from IPython.display import display
import warnings
idx = pd.IndexSlice

# ## Define cases
# Global path variables

will_elec = ["01-battery-distributed", "02-battery-overnight", "03-battery-recharging", "04-battery-solar"]
will_h2 = ["01-h2-distributed", "02-h2-overnight", "03-h2-recharging", "04-h2-solar"]
sdewes_ap = ["base", "high", "low", "med", "ev-med"]
        
scenario_dict = {
    "will": [will_h2, [2020, 2035, 2050]],
    "sdewes-ap": [sdewes_ap, [2020, 2030, 2040, 2050]]
    }

group_name="sdewes-ap"
indx=0

files_lst = scenario_dict[group_name][0]
yrs_sel = scenario_dict[group_name][1] 
yrs_sel_str = [f"{year}" for year in yrs_sel]
yrs_str='-'.join([str(item) for item in yrs_sel])
yrs_to_calc = [2020, 2025, 2030, 2035, 2040, 2045, 2050]
    
# Define paths as global variables
path_base = "C:/Local/REMix"
path_input = f"{path_base}/remix_nz/input"
path_demand = f"{path_input}/demand/{group_name}"
path_profiles = f"{path_input}/profiles"      # renewables
path_brownfield = f"{path_input}/brownfield"  # info hydro and existing power plants database
path_geo = f"{path_input}/shapefiles"         # geojson
geofile="11regionsNZ.geojson"
shp_file = f"{path_geo}/11regionsNZ.shp"


demand_file=files_lst[indx] 

case_name=f"{demand_file}_{yrs_str}"
path_output = f"{path_base}/remix_nz/output/{group_name}"
data_dir = Path(f"{path_output}/{case_name}/data")
results_dir = Path(f"{path_output}/{case_name}/result")
csv_dir = Path(f"{path_output}/csv")
csv_dir.mkdir(parents=True, exist_ok=True)

# read in the output `*.gdx` file from the optimization in GAMS
results = GDXEval(f"{results_dir}/{case_name}.gdx")
rename_techs = {
    "pv_central_fixed": "PV", #open area (fixed)",
    "pv_central_track_azimuth": "PV", #open area (tracking)",
    "wind_onshore": "Onshore wind",
    "wind_offshore_floating": "Offshore wind",# (floating)",
    "wind_offshore_foundation": "Offshore wind",# (foundation)",
    "BIO": "Biomass",
    "DIE": "Diesel",
    "GT": "Gas",
    "Battery": "Battery",
    "CCGT": "Gas",
    "COAL": "Coal",
    "OCGT": "Gas",
    "geoth": "Geothermal",
    "Hydro": "Hydro"
    }

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
shp_attrcol = "id"
idx = pd.IndexSlice 
lat = [-47.758239, -34.111702]
lon = [165.692103, 179.050919]

centroids = coord_dict
shp_attrcol = "id"


# %% [markdown]
# Export results csv
def results_to_csv(file=demand_file):
    
    # Converter capacities
    converter_caps = results["converter_caps"]   # convert converter capacities to a Pandas DataFrame
    converter_caps = converter_caps[converter_caps > 0.01].dropna()  # remove all capacities with less than 10 MW
    converter_caps.to_csv(f"{results_dir}/{case_name}_converter_caps.csv")
    # Example of slicing:
    # accNodesModel=all, accYears=2050 only, techs=all, commodit=only Elec, capType=total
    # caps.loc[idx[:, "2050", :, "Elec", "total"], :].round(2) 

    # Storage capacities
    storage_caps = results["storage_caps"]   # convert storage capacities to a Pandas DataFrame
    storage_caps = storage_caps[storage_caps > 0.01].dropna()  # remove all capacities with less than 10 MW
    storage_caps.to_csv(f"{results_dir}/{case_name}_storage_caps.csv")

    # Commodity balance annual
    commodity_balance_annual = results["commodity_balance_annual"]   # convert commodity_balance to a Pandas DataFrame
    commodity_balance_annual.to_csv(f"{results_dir}/{case_name}_commodity_balance_annual.csv")

    # Marginals
    marginals_sourcesink_profile = results["marginals_sourcesink_profile"]   # convert marginals to a Pandas DataFrame
    marginals_sourcesink_profile.to_csv(f"{results_dir}/{case_name}_marginals_sourcesink_profile.csv")

    # Accounting indicators
    indicator_accounting = results["indicator_accounting"]   # convert  to a Pandas DataFrame
    indicator_accounting.to_csv(f"{results_dir}/{case_name}_indicator_accounting.csv")
#results_to_csv(demand_file)

def manage_converter_caps(file=demand_file):
    # Converter capacities
    caps = results["converter_caps"]   # convert converter capacities to a Pandas DataFrame
    nodes = [n for n in caps.index.get_level_values(0) if n != "global"]
    caps = (caps[caps.abs() > 0.01]
            .dropna()  # remove all capacities with less than 10 MW
            .rename(index=rename_techs))
    caps = caps.loc[idx[nodes, yrs_sel_str, :, "Elec", "total"], :].round(2) 
    elec = (caps.groupby(level=["accNodesModel", "accYears", "techs"], observed=True)
                  .sum()
                  .round(2))
                  #.unstack(["techs"]))
    return elec
caps=manage_converter_caps()
caps.to_csv(f"{csv_dir}/{case_name}_converter_caps.csv")


def manage_cba(file=demand_file):
    # Converter capacities
    cba = results["commodity_balance_annual"]   # convert converter capacities to a Pandas DataFrame
    nodes = [n for n in cba.index.get_level_values(0) if n != "global"]
    cba = (cba[cba.abs() > 0.01]
           .loc[idx[nodes, yrs_sel_str, :, "Elec", "netto"], :]
           .dropna()  # remove all capacities with less than 10 MW
           .rename(index=rename_techs))

    # Slicing:
    # accNodesModel=all but global, accYears= only selected, techs=all, commodit=only Elec, capType=total
    unstacked = (cba.round(2) 
            .groupby(["accYears", "accNodesModel", "techs"])
            .sum()
            .unstack(["techs"]))
    return unstacked
cba=manage_cba()
cba.to_csv(f"{csv_dir}/{case_name}_commodity_balance_annual.csv")

def manage_commodity_balance_annual(file=demand_file):
        # Converter capacities
    converter_caps = results["converter_caps"]   # convert converter capacities to a Pandas DataFrame
    converter_caps = converter_caps[converter_caps > 0.01].dropna()  # remove all capacities with less than 10 MW
    converter_caps.to_csv(f"{results_dir}/{case_name}_converter_caps.csv")
    # Example of slicing:
    # accNodesModel=all, accYears=2050 only, techs=all, commodit=only Elec, capType=total
    # caps.loc[idx[:, "2050", :, "Elec", "total"], :].round(2) 

    # Storage capacities
    storage_caps = results["storage_caps"]   # convert storage capacities to a Pandas DataFrame
    storage_caps = storage_caps[storage_caps > 0.01].dropna()  # remove all capacities with less than 10 MW
    storage_caps.to_csv(f"{results_dir}/{case_name}_storage_caps.csv")

    # Commodity balance annual
    commodity_balance_annual = results["commodity_balance_annual"]   # convert commodity_balance to a Pandas DataFrame
    commodity_balance_annual.to_csv(f"{results_dir}/{case_name}_commodity_balance_annual.csv")

    # Marginals
    marginals_sourcesink_profile = results["marginals_sourcesink_profile"]   # convert marginals to a Pandas DataFrame
    marginals_sourcesink_profile.to_csv(f"{results_dir}/{case_name}_marginals_sourcesink_profile.csv")

    # Accounting indicators
    indicator_accounting = results["indicator_accounting"]   # convert  to a Pandas DataFrame
    indicator_accounting.to_csv(f"{results_dir}/{case_name}_indicator_accounting.csv")
#results_to_csv(demand_file)



def select_techs(year,tech,file=demand_file):
    # Converter capacities
    caps = results["converter_caps"]   # convert converter capacities to a Pandas DataFrame
    
    caps = (caps[caps.abs() > 0.01]
            .dropna()  # remove all capacities with less than 10 MW
            .rename(index=rename_techs))

    # techs = list(set(cba.index.get_level_values(2)))
    # techs_pv = [i for i in techs if i.lower().startswith("pv")]
    # techs_wind = [i for i in techs if i.lower().startswith("wind")]

    sliced_df = (caps.loc[(slice(None), f"{str(year)}", f"{tech}"), :]
               .groupby(level="accNodesModel")
               .sum()
               .round(2)
               .loc[lambda x: x.index != 'global'])
    return sliced_df

pv_2020=select_techs(2020,"PV")
pv_2030=select_techs("2030","PV")
pv_2040=select_techs("2040","PV")
pv_2050=select_techs("2050","PV")
#wind_2020=select_techs(2020,"PV")



def plot_choropleth(data, shp_file, shp_attrcol, lat, lon, title="", clabel="", cmap="Oranges", cbar=False, ax=None):
    gdf = gpd.read_file(shp_file)
    if ax is None:
        fig, ax = plt.subplots()
    gdf.plot(ax=ax, color='lightgrey', edgecolor='black')
    data.plot(ax=ax, column=data.rename, legend=cbar, cmap=cmap)
    ax.set_title(title)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    if cbar:
        ax.set_aspect('equal', adjustable='datalim')
        ax.get_figure().colorbar(data.plot(column=data.name, cmap=cmap, ax=ax).collections[0], ax=ax)
    plt.show()




fig = plt.figure(figsize=figsize_dual)
fig.patch.set_facecolor("#F0F8FF")
ax1 = fig.add_subplot(121)
ax1.set_facecolor("#F0F8FF")
plot_choropleth(
    pv_2030,   #df
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
    pv_2040,
    shp_file,
    shp_attrcol,
    lat,
    lon,
    title="Power generation from wind",
    clabel="Energy in TWh",
    cmap="Blues",
    ax=ax2,
)



def plot_cap(year,tech):
    df=select_techs(year,tech)
    fig = plot_choropleth(
        df,
        shp_file,
        shp_attrcol,
        lat,
        lon,
        title=f"{tech} {year}",
        clabel="Installed capacity in GW",
        cmap="Oranges",
        #cbar=bar_bol,
        #ax=ax1
        #,limits=[0,max_pv*1.01]
        )
    return fig
plot_cap(2030, "PV")
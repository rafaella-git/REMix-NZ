# %% [markdown]
# # Part C: evaluation of results
# This script helps you visualise the results obtained from the model run in part B.
# [Acces DLR tutorial](https://dlr-ve.gitlab.io/esy/remix/framework/dev/getting-started/tutorials/tutorial-103.html)

# git test

# %%
# Import dependencies
from remix.framework.tools.gdx import GDXEval
from remix.framework.tools.plots import plot_network, plot_choropleth
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
import warnings

# %%
# Define paths
path_base = "C:/Local/REMix"  
path_input = f"{path_base}/remix_nz/input"
path_output = f"{path_base}/remix_nz/output" 
path_geo = f"{path_input}/shapefiles"      # geojson

# %%
# Define useful files
geofile="11regionsNZ.geojson"
shp_file = f"{path_geo}/11regionsNZ"

# %%
# Define useful shortcuts
shp_attrcol = "id"
idx = pd.IndexSlice #often used shortcut


#%% [markdown]
# ## User interactive
# There are two ways to access the results folders we are interested in
# Option 1 : defining a list with demand files and a list of 
# either a lists of lists of years, or a single list of years
demand_files=["nz_profile_11nodes"]
yrs=[2020,2030,2040,2050] # years optimised
#demand_files=["medpop_evs","low_pop_out","med_pop_out","high_pop_out"] #493 as in the course 493, this is for Liv and Sam
#yrs=[2030,2040,2050] 

#%% [markdown]
# Option 2 (might come in handier later): is creating a list with all the folders
# if folder_lst is empty, the code takes option 1 as default
folder_lst=[]


#%%
# Define figure [arameters]
lat = [-47.758239, -34.111702]
lon = [165.692103, 179.050919]
# params
plt.rcParams.update({"figure.autolayout": True})  # use tight_layout
plt.rcParams.update({"figure.titlesize": 20})  # size subtitle
plt.rcParams.update({"figure.dpi": 75})  # size subtitle
plt.rcParams.update({"savefig.dpi": 300})  # dpi
plt.rcParams.update({"font.size": 16})  # default font size
plt.rcParams.update({"axes.titlesize": 20})  # title size per subplot
plt.rcParams.update({"axes.labelsize": 18})  # label size colormap
#figsize_dual = (13.0, 6)
#figsize_single = (7.5, 6)

# %% [markdown]
# ### Define useful functions
# 1. Define datafame from gdx results file
# 1.a Define a list of folder names aka case names

def get_case_lst(cases=demand_files,years=yrs):
    # if folder_lst is defined as a not empty list it takes it
    if bool(folder_lst):
        case_lst=folder_lst
    # if we need to create the folder list from demand files and years
    else:
        case_lst=[]
        for element in years:
            # Check if each element is also a list
            if not isinstance(element, list):
                yrs_str='-'.join([str(item) for item in yrs])
                for case in cases:
                    case_name=f"{case}_{yrs_str}"
                    case_lst.append(case_name)
            else:
                if len(cases) != len(years):
                    warnings.warn("Lists of years must match number of demand files.")
                indx=0
                for case in cases:
                    yrs_str='-'.join([str(item) for item in yrs[indx]])
                    case_name=f"{case}_{yrs_str}"
                    case_lst.append(case_name)
                    indx+=1
    return case_lst
case_lst=get_case_lst()



#%% [markdown]
# 1.b Define a dictionary that matches the name of the case,
#  with its respective result file `.gdx`

def results_dict(folder_lst):
    # read in the output `*.gdx` file from the optimization in GAMS
    results_dict = {}  # Create a dictionary to store the results
    for name in folder_lst:
        path_result = f"{path_output}/{name}/result" 
        results_dict[name] = GDXEval(f"{path_result}/{name}.gdx")
    return results_dict


results=results_dict(get_case_lst())

# %% [markdown]
# ### 2. Capacities

def caps_dict(results_d=results,create_csv_orig=False,create_csv_grouped=False):
    cap = {}  # Create a new dictionary to store the modified dataframes
    for key, result in results_d.items():
        # Create a new dataframe with only converter cap
        cap_gdx = result["converter_caps"]
        cap[key] = cap_gdx
        path_result = f"{path_output}/{key}/result" 
        if create_csv_orig:
            cap_gdx.to_csv(f"{path_result}/{key}_caps.csv")
        #cap_gdx = cap_gdx[cap_gdx['value'] > 0.01].dropna() # remove all capacities with less than 10 MW
        # Replace 'techs' values based on mapping
        tech_list=list(set(cap_gdx.index.get_level_values('techs')))
        tech_map = {tech: "PV" if tech.startswith("pv") else "Wind" if tech.startswith("wind") else tech for tech in tech_list}
        # Reset the index
        cap_gdx = cap_gdx.reset_index()
        # Apply the mapping to the 'techs' column
        cap_gdx['techs'] = cap_gdx['techs'].map(tech_map)
        # Set the index back
        cap_gdx = cap_gdx.set_index(['accNodesModel', 'accYears', 'techs', 'commodity', 'capType'])
        cap_gdx = cap_gdx.groupby(cap_gdx.index.names).agg({'value': 'sum'})
        # Reset the index to have a flat structure
        cap_gdx.reset_index(inplace=True)
        cap_gdx.set_index(['accNodesModel', 'accYears', 'techs', 'commodity', 'capType'], inplace=True)
        index_labels = cap_gdx.index
        cap_gdx = cap_gdx.loc[idx[:, :, :, "Elec", "total"], :].round(2)
        if create_csv_grouped:
            cap_gdx.to_csv(f"{path_result}/{key}_caps_grouped.csv")
        cap[key] = cap_gdx
        return cap

caps=caps_dict(results,True,True)

#%% [markdown]
# Finds the maximum value for a tech in the caps df

def max_cap(tech,caps_dict):
    max_values = {}  # Create a new dictionary to store maximum values
    nodes=["NIS","AKL","WTO","TRN","BOP","HBY","CEN","WEL","NEL","CAN","OTG"]
    for key, cap in caps_dict.items():
        sliced = cap.loc[idx[nodes, :, f"{tech}", "Elec", "total"], idx[:]].round(2)#.div(1e3)
        max = round(sliced['value'].max(), 2)
        max_values[key] = max
    return max_values

max_pv=max_cap("PV",caps)
max_wind=max_cap("Wind",caps)

# %% [markdown]
# Plot capacities | (dataframe caps, tech you want to plot, year)
# The structure of this fucntion differs from the ones above because
#  it only plots the capacity 
#  idk fix this idk ifit takes the dict or not
def plot_cap(tech,year,caps_dict=caps):
    figs = {}  # Create a new dictionary to store figures
    for key, cap in caps_dict.items():
        # slice capacities data for specified tech
        map_pv=cap.loc[idx[nodes_lst, f"{year}",  f"{tech}", "Elec", "total"],  idx[:]].groupby("accNodesModel").sum().round(2) #.div(1e3)
        fig = plot_choropleth(
        map_pv,   #df
        shp_file,
        shp_attrcol,
        lat,
        lon,
        #title=f"PV {year}",
        clabel="Photovoltaic capacity in GW",
        cmap="Oranges",
        cbar=bar_bol,
        ax=ax1
        #,limits=[0,max_pv*1.01]
        )
        figs[key] = fig
    return figs




def plot_caps_v(df_caps,year,bar_bol=False, case=demand_file):
    # slice capacities dataframe for specified techs
    map_pv=df_caps.loc[idx[nodes_lst, f"{year}",  "PV", "Elec", "total"],  idx[:]].groupby("accNodesModel").sum().round(2)
    map_wind=df_caps.loc[idx[nodes_lst, f"{year}", "Wind", "Elec", "total"],  idx[:]].groupby("accNodesModel").sum().round(2)
    # overall figure
    fig = plt.figure(figsize=figsize_dual_inv_short)
    if bar_bol: #modify to make them match
        fig = plt.figure(figsize=figsize_dual_inv)
    fig.patch.set_facecolor(background_color)
    # subfigure ax1
    ax1 = fig.add_subplot(211)
    #ax1.set_facecolor("#F0F8FF")
    plot_choropleth(
        map_pv,   #df
        shp_file,
        shp_attrcol,
        lat,
        lon,
        #title=f"PV {year}",
        clabel="Photovoltaic capacity in GW",
        cmap="Oranges",
        cbar=bar_bol,
        ax=ax1, 
        limits=[0,4.48*1.01]
        #limits=[0,max_pv*1.01]
    )
    ax2 = fig.add_subplot(212)
    #ax2.set_facecolor("#F0F8FF")
    plot_choropleth(
        map_wind,
        shp_file,
        shp_attrcol,
        lat,
        lon,
        #title=f"Wind {year}",
        clabel="Wind capacity in GW",
        cmap="Blues",
        cbar=bar_bol,
        ax=ax2,
        limits=[0,0.35*1.03]        
        #limits=[0,max_wind*1.03] #1.03
    )
    fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
    ax1.text(0.02, 0.95, f"  {year}", fontsize=22, ha='left', va='top', color='gray', transform=ax1.transAxes)
    ax1.text(0.95, 0.01, f"{case}", fontsize=12, ha='right', va='bottom', color='gray', transform=ax1.transAxes)
    ax2.text(0.02, 0.95, f"  {year}", fontsize=22, ha='left', va='top', color='gray', transform=ax2.transAxes)
    ax2.text(0.95, 0.01, f"{case}", fontsize=12, ha='right', va='bottom', color='gray', transform=ax2.transAxes)
    # define each figure directoy
    mini_fig=f"{minifig_file}_cap_v_{year}.png"
    plt.savefig(mini_fig, bbox_inches='tight', pad_inches=0.1)



# %% [markdown]
# ### 3. Plot transfer (like triangles) | dataframe transfer flows/balance?

# %% [markdown]
# ### 4. Plot ISO plots x=day y=hour heat=intensity | dataframe balance

# %% [markdown]
# ### 5. Define function to export a figure | (figure, name, destiny folder) 


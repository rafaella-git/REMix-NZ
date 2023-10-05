#%% [markdown]
# # Evaluate results for 493


# %%[markdown]
# ### Import dependencies
from remix.framework.tools.gdx import GDXEval
from remix.framework.tools.plots import plot_network, plot_choropleth
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
import os
idx = pd.IndexSlice #often used shortcut

# %%[markdown]
# ### Define relevant lists to locate files
# Demand files available (different scenarios)
files_lst=["medpop_evs_base","low_pop_out_base","med_pop_out_base","high_pop_out_base"] #493 as in the course 493, this is for Liv and Sam
files_lst=["medpop_evs","low_pop_out","med_pop_out","high_pop_out"] 

# Years run
yrs=[2020,2030,2040,2050] 
yrs=[2030,2040,2050]
yrs_str='-'.join([str(item) for item in yrs])

# %%[markdown]
# ### Define directories
path_base = "C:/Local/REMix"  
path_input = f"{path_base}/remix_nz/input"
path_output = f"{path_base}/remix_nz/output" 
path_geo = f"{path_input}/shapefiles"      
geofile="11regionsNZ.geojson"
shp_file = f"{path_geo}/11regionsNZ"
shp_attrcol = "id"

# %%[markdown]
# ### Recover results of the model runs
# The function `results_lsts` allows us to retrieve the capacities
# and the annual balance for the commodities, for each of the files.
# Once those dataframes are defined inside of lists, we can interate
# over the lists to plot.
caps=[]
balance=[]
#FIX part 1: these may be irrelevant for 493
#transfer_caps=[]
#transfer_balance=[]

def results_lsts(lst=files_lst):
    indx=0
    # ### Technology mapping
    # Grouping pv and wind technologies together
    pv_techs= ["pv_central_fixed", "pv_central_track_azimuth", "pv_decentral"]
    wind_techs=["wind_onshore","wind_offshore_floating","wind_offshore_foundation"]
    # Create a mapping dictionary to update techs
    pv_tech_map = {tech: "pv" for tech in pv_techs}
    wind_tech_map = {tech: "wind" for tech in wind_techs}
    for file in files_lst:
        demand_file=file 
        case_name=f"{demand_file}_{yrs_str}"
        path_result = f"{path_output}/{case_name}/result" 
        gdx_file = f"{path_result}/{case_name}.gdx"
        if os.path.exists(gdx_file):
            results = GDXEval(gdx_file)       
            # Create dataframes and append to each respective list 

            # Installed capacity (generation)
            caps_gdx = results["converter_caps"]
            caps_gdx = caps_gdx[caps_gdx > 0.01].dropna()  # Remove all capacities with less than 10 MW
            # Replace techs in the DataFrame based on the mapping
            caps_gdx['techs'] = caps_gdx['techs'].replace(pv_tech_map).replace(wind_tech_map)
            # Group by 'accNodesModel', 'accYears', 'commodity', and 'capType'
            # Sum the 'Value' column for each group
            caps_gdx = caps_gdx.groupby(['accNodesModel', 'accYears', 'techs','commodity', 'capType'])['Value'].sum().reset_index()
            # Get only commodity=electricity and capType=total
            caps_gdx.loc[idx[:, :, :, "Elec", "total"], :].round(2)        
            caps.append(caps_gdx)

            # Commodity annual balance (generation)
            balance_gdx = results["commodity_balance_annual"]
            balance_gdx = balance_gdx[balance_gdx  > 0.01].dropna()  # Remove all capacities with less than 10 MW
            balance_gdx['techs'] = caps_gdx['techs'].replace(pv_tech_map).replace(wind_tech_map)
            balance_gdx = balance_gdx.groupby(['accNodesModel', 'accYears', 'techs','commodity', 'balanceType'])['Value'].sum().reset_index()
            balance_gdx.loc[idx[:, :, :, "Elec", "netto"], :].round(2) 
            balance.append(balance_gdx)

            # FIX part 2: these may be irrelevant for 493
            # Transfer capacity (between regions)
            #transfer_caps_gdx = results["transfer_caps"]
            #transfer_caps.append(transfer_caps_gdx)
            # Annual energy flows between model regions
            #transfer_balance_gdx = results["transfer_flows_annual"]
            #transfer_balance.append(transfer_balance_gdx)   
            print(f"Success: {case_name} is in position {len(caps)} of the list.")
        else:
            print(f"GDX file not found for {case_name}: you need to run that model first.")
        indx+=1   
results_lsts()

if len(caps)==len(files_lst):
    print("Success: the lists are in order and ready to be plotted!")
else:
    print("Oh no! check your list, not all the GAMS files were available")



# %%[markdown]
# ### Run debugger -----


techs_colors = {
    "pv": "#ffc000",
    "wind": "#9dc3e6",
    "HV": "#70ad47",
    "CCGT": "#ff5f2d",
}
background_color = "#fffaebff"



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
#limits of the plot 
lat = [-47.758239, -34.111702]
lon = [165.692103, 179.050919]
centroids = coord_dict


# Set matplotlib parameters
def set_matplotlib_parameters():
    plt.rcParams.update({"figure.autolayout": True})  # use tight_layout
    plt.rcParams.update({"figure.titlesize": 20})  # size subtitle
    plt.rcParams.update({"figure.dpi": 75})  # size subtitle
    plt.rcParams.update({"savefig.dpi": 300})  # dpi
    plt.rcParams.update({"font.size": 16})  # default font size
    plt.rcParams.update({"axes.titlesize": 20})  # title size per subplot
    plt.rcParams.update({"axes.labelsize": 18})  # label size colormap
set_matplotlib_parameters()
figsize_dual = (13.0, 6)
figsize_single = (7.5, 6)




def plot_renewable_gen(balance_df=balance[2],year=2050,case_name=files_lst[2]):

    # Exclude global (we are only plotting nodes)
    balance_df=balance_df.loc[idx[nodes_lst, str(year), :, "Elec", "netto"]]
    # Extract pv and wind
    map_pv = balance_df.loc[idx[:, :, "pv"], idx[:]].groupby("accNodesModel").sum().div(1e3)
    map_wind = balance_df.loc[idx[:, :, "wind"], idx[:]].groupby("accNodesModel").sum().div(1e3)

    fig = plt.figure(figsize=figsize_dual)
    fig.patch.set_facecolor(background_color)
    fig.suptitle(f"Renewable generation for {year}", fontsize=20)

    # Plot solar generation
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
    
    # Plot wind generation
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
    # Add a footnote
    fig.text(0.5, 0.01, f"Case: {case_name}", fontsize=12, ha="center", va="center")
    # Save the figure
    save_path = 'C:/path/to/your/directory/figure.svg'
    plt.savefig(save_path)  
    plt.savefig(f"{file_name}.{format}", format=format)
    

# Plot generation for 2050 with med_pop_base
plot_renewable_gen()
# Plot generation for 2040 with med_pop_base_evs
indx=0
plot_renewable_gen(balance[indx], 2040, files_lst[indx])
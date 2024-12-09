#evaluate
# Import dependencies
import gdxpds
from remix.framework.tools.gdx import GDXEval
from remix.framework.tools.plots import plot_network, plot_choropleth
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
import numpy as np


name="h2pos_2020-2030-2040-2050"

# Define paths
path_base = "C:/Local/REMix"  
path_input = f"{path_base}/remix_nz/input"
path_output = f"{path_base}/remix_nz/output/mbie" 
path_geo = f"{path_input}/shapefiles"      # geojson
# Define useful files
geofile="11regionsNZ.geojson"
shp_file = f"{path_geo}/11regionsNZ"
# Define useful shortcuts
shp_attrcol = "id"
idx = pd.IndexSlice #often used shortcut

# Load data from .gdx file
path_result = f"{path_output}/{name}/result" 
gdx_file = f"{path_result}/{name}.gdx"
#results_dict[name] = GDXEval(f"{path_result}/{name}.gdx")
data = gdxpds.to_dataframes(gdx_file)

# Extract relevant data
converter_caps = data['converter_caps']
commodity_balance_annual = data['commodity_balance_annual']
storage_flows = data['storage_flows']

# Define nodes, technologies, and commodities
nodes_lst = ["NIS","AKL","WTO","TRN","BOP","HBY","CEN","WEL","NEL","CAN","OTG"]
commodities = ["Elec", "H2"]

def plot_installed_capacities(data, years):
    dfs = []
    for year in years:
        df = data.loc[pd.IndexSlice[:, 'global', year, :, :], :]
        df = df.groupby(level=['year', 'technology']).sum().unstack(level=1)
        dfs.append(df)
    combined_df = pd.concat(dfs)
    combined_df.columns = combined_df.columns.droplevel(0)
    combined_df.plot(kind='bar', stacked=True)
    plt.ylabel('Installed Capacities (GW)')
    plt.xlabel('Year')
    plt.legend(title='Technology')
    plt.title('Installed Capacities Over Years')
    plt.show()

# Plot 2: Bar stacked plot for energy output
def plot_energy_output(data, years):
    dfs = []
    for year in years:
        df = data.loc[pd.IndexSlice[:, 'global', year, :, 'net'], :]
        df = df.groupby(level=['year', 'technology']).sum().unstack(level=1)
        dfs.append(df)
    combined_df = pd.concat(dfs)
    combined_df.columns = combined_df.columns.droplevel(0)
    combined_df.plot(kind='bar', stacked=True)
    plt.ylabel('Energy output (TWh)')
    plt.xlabel('Year')
    plt.legend(title='Technology')
    plt.title('Energy Output Over Years')
    plt.show()

# Plot 3: Geospatial plot for installed capacities and commodity balance
def plot_geospatial(data, geofile, year, value_type, title):
    gdf = gpd.read_file(geofile)
    df = data.loc[pd.IndexSlice[nodes_lst, year, :, :, value_type], :]
    df = df.groupby(level='node').sum().reset_index()
    df['coords'] = df['node'].map(dict(zip(nodes_lst, nodes_coords)))
    gdf = gdf.merge(df, left_on='id', right_on='node')
    
    gdf.plot(column='value', cmap='OrRd', legend=True, legend_kwds={'label': title})
    plt.title(f'{title} in {year}')
    plt.show()

# Plot 4: ISO plot for battery operation
def plot_iso(data, year):
    df = data.loc[pd.IndexSlice[:, 'global', year, 'Battery', 'Elec_battery'], :]
    df.index = df.index.get_level_values(0)
    df = df.unstack(level=0)
    
    plt.imshow(df.values, aspect='auto', cmap='viridis')
    plt.colorbar(label='Battery Operation (GW)')
    plt.xlabel('Day of Year')
    plt.ylabel('Hour of Day')
    plt.title(f'Battery Operation in {year}')
    plt.show()

# Define coordinates
nodes_lat = [-35.8758611, -36.9547333, -38.4192333, -39.3342222, -37.9867861, -39.5516472, -40.2813000, -41.1502722, -41.6735722, -43.8619833, -45.481222]
nodes_lon = [174.4669472, 174.8625250, 175.8000111, 174.3204333, 176.8294056, 176.8208167, 175.6404750, 174.9811500, 172.8737917, 171.3427694, 169.3195194]
nodes_coords = [(lat, lon) for lat, lon in zip(nodes_lat, nodes_lon)]


nodes_lst=["NIS","AKL","WTO","TRN","BOP","HBY","CEN","WEL","NEL","CAN","OTG"]
nodes_lat=[-35.8758611, -36.9547333, -38.4192333, -39.3342222, -37.9867861, -39.5516472, -40.2813000, -41.1502722, -41.6735722, -43.8619833, -45.481222]
nodes_lon=[174.4669472, 174.8625250, 175.8000111, 174.3204333, 176.8294056, 176.8208167, 175.6404750, 174.9811500, 172.8737917, 171.3427694, 169.3195194]    
nodes_coords=[(nodes_lat[i],nodes_lon[i]) for i in range(0,len(nodes_lon))]
# dictionary with 'Node': [lat,long]
coord_dict = dict.fromkeys(nodes_lst)
for key, value in zip(coord_dict.keys(), nodes_coords):
    coord_dict[key] = value
print(coord_dict)







# Example usage
years = ['2030', '2040', '2050']
plot_installed_capacities(converter_caps, years)
plot_energy_output(commodity_balance_annual, years)
plot_geospatial(converter_caps, '11regionsNZ.geojson', '2030', 'value', 'Installed Capacities')
plot_iso(storage_flows, '2030')
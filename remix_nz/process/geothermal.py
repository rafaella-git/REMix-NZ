
# 

import pandas as pd
from pathlib import Path


# path definition
path_base = "C:/Local/REMix"
path_input = f"{path_base}/remix_nz/input"
path_brownfield = f"{path_input}/brownfield" 
#path_output = f"{path_base}/remix_nz/output/will" 
#path_demand = f"{path_input}/demand/will"    
#path_profiles = f"{path_input}/profiles"  
#path_geo = f"{path_input}/shapefiles"   





def load_brownfield():
    # Specify columns to be read from the CSV
    columns_used = ['Type','Primary_fuel','Capacity_MW', 'Status', 'Year', 'Connection_type', 'Avg_Ann_Gen_GWh', 'GIP substation', 'lat', 'long']
    
    # Read only the specified columns and turn it into a df
    df = pd.read_csv(Path(f"{path_brownfield}/power-plant-db-PM.csv"), usecols=columns_used)

    # If Type is 'Thermal', change the value in column 'Thermal' for the value in 'Primary_fuel'
    df.loc[df['Type'] == 'Thermal', 'Type'] = df.loc[df['Type'] == 'Thermal', 'Primary_fuel']

    # For every value in the column 'Capacity_MW', divide it by 1000 and then change the name of the column to 'Capacity_GW'
    df['Capacity_GW'] = df['Capacity_MW'] / 1000
    df.drop(columns=['Capacity_MW'], inplace=True)  # Drop the original column
    
    # Isolate rows where 'long' is NaN and 'lat' contains a comma
    mask = df['long'].isna() & df['lat'].str.contains(',')

    # Split lat and long values, keeping original lat values for non-affected rows
    df.loc[mask, ['lat', 'long']] = df.loc[mask, 'lat'].str.split(',\s+', n=1, expand=True)

    # Filter only the ones where "Status" is "Operational"
    


    return df

    # Type: primary fuel
    #  Geothermal: binary, double flash, dry steam, single flash, triple flash, binary
    #  Hydroelectric: Water
    #  Thermal: Biogas, Biomass, Coal, Diesel, Natural gas, Unknown, Waste heat, Wood, Wood waste
    #  Wind: Wind


    # Capacity: from MW to GW

    # Connection type
    #   embedded = connected to a local distribution network
    #   grid = connected directly to the main grid infrastructure
    #   partially embedded = connected to both




load_brownfield()






















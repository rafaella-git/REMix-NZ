
def add_hydro(m):

    # Configure how water exits the system (flowing to the ocean).
    sourcesink_config = pd.DataFrame(index=pd.MultiIndex.from_product([m["Base"].set.nodesdata, m["Base"].set.yearssel, ["Ocean"], ["Water_out"]]))
    # Upper profile constraints: we need negative values to get water out of the system
    sourcesink_config.loc[idx[m["Base"].set.nodesdata, :, :, :], "usesUpperProfile"] = 1 
    sourcesink_config.dropna(inplace=True)
    m["Base"].parameter.add(sourcesink_config, "sourcesink_config")
    sourcesink_config

    hydro_vintage = [1950]
    hydro_years= [2000]+yrs_to_calc
    # hydro_techs = ["Hydroplant"] 
    # hydro_nodes = ["BOP", "CAN", "CEN", "HBY", "NEL", "OTG", "WTO"] 
    
    hydro_techs = ["Hydro_Clyde", "Hydro_Roxburgh"] 
    hydro_nodes = ["OTG"]
    hydro_commodities = ["Elec", "Dunstan_in", "Roxburgh_in", "Water_out"]
    hydro_activities = ["Power_gen","Spill"] 

	# Converter (turbine) lifetime and activity limits
    conv_tech = pd.DataFrame(index=pd.MultiIndex.from_product([hydro_techs, hydro_vintage]))
    conv_tech.loc[idx[:, :], ["lifeTime"]] = 100
    conv_tech.loc[idx[:, :], ["activityUpperLimit"]] = 1
    m["Base"].parameter.add(conv_tech, "converter_techparam")

    conv_cap = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [hydro_nodes, hydro_years, hydro_techs] #FIXME: maybe years_mentioned
        )
    )
    

    # Installed turbine capacities for each node and year
    # conv_cap.loc[idx[["BOP"], [2000], "Hydro"], "unitsBuild"] = 0.17095  # GW_el
    # conv_cap.loc[idx[["CAN"], [2000], "Hydro"], "unitsBuild"] = 1.82683  # GW_el
    # conv_cap.loc[idx[["CEN"], [2000], "Hydro"], "unitsBuild"] = 0.399  # GW_el
    # conv_cap.loc[idx[["HBY"], [2000], "Hydro"], "unitsBuild"] = 0.1422  # GW_el
    # conv_cap.loc[idx[["NEL"], [2000], "Hydro"], "unitsBuild"] = 0.0453  # GW_el
    # conv_cap.loc[idx[["OTG"], [2000], "Hydro"], "unitsBuild"] = 1.664  # GW_el
    # conv_cap.loc[idx[["WTO"], [2000], "Hydro"], "unitsBuild"] = 1.0873  # GW_el
    conv_cap.loc[idx[["OTG"], [2000], "Hydro_Clyde"], "unitsBuild"] = 464 * 0.001 # GW_el   
    conv_cap.loc[idx[["OTG"], [2000], "Hydro_Roxburgh"], "unitsBuild"] = 320 * 0.001   # GW_el   

    conv_cap.loc[idx[:, :, :], "noExpansion"] = 1  # Prevent expansion. Boolean.
    m["Base"].parameter.add(conv_cap, "converter_capacityparam")

    # Turbine conversion coefficients (input-output relationships).
    conv_coef = pd.DataFrame(index=pd.MultiIndex.from_product([hydro_techs, hydro_vintage, hydro_activities,hydro_commodities]))

    # Define coefficients for electricity generation.  , "Roxburgh_in"
    conv_coef.loc[idx["Hydro_Clyde", :, "Power_gen", "Elec"], "coefficient"] = 0.5181 * 0.001 # GW_el
    conv_coef.loc[idx["Hydro_Clyde", :, "Power_gen", "Dunstan_in"], "coefficient"] = -1  # GW_el
    conv_coef.loc[idx["Hydro_Clyde", :, "Power_gen", "Roxburgh_in"], "coefficient"] = 1 # GW_el
    conv_coef.loc[idx["Hydro_Roxburgh", :, "Power_gen", "Elec"], "coefficient"] = 0.4016 * 0.001 # GW_el
    conv_coef.loc[idx["Hydro_Roxburgh", :, "Power_gen", "Roxburgh_in"], "coefficient"] = -1  # GW_el
    conv_coef.loc[idx["Hydro_Roxburgh", :, "Power_gen", "Water_out"], "coefficient"] = 1 # GW_el   
    # #FIXME: DIRTY hack TO MATCH 2024 TO REAL DATA MBIE
    # conv_coef.loc[idx[:, :, "Power_gen", "Elec"], "coefficient"] = 1 # GW_el
    # conv_coef.loc[idx[:, :, "Power_gen", "Water_in"], "coefficient"] = -0.95  # GW_el
    # conv_coef.loc[idx[:, :, "Power_gen", "Water_out"], "coefficient"] = 0.95 # GW_el


    # Define coefficients for bypassing water without generating electricity. Spill is not limited by the capacity of the turbine
    conv_coef.loc[idx["Hydro_Clyde", :, "Spill", "Dunstan_in"], "coefficient"] = -100  # GW_el
    conv_coef.loc[idx["Hydro_Clyde", :, "Spill", "Roxburgh_in"], "coefficient"] = 100 # GW_el
    conv_coef.loc[idx["Hydro_Roxburgh", :, "Spill", "Roxburgh_in"], "coefficient"] = -100  # GW_el
    conv_coef.loc[idx["Hydro_Roxburgh", :, "Spill", "Water_out"], "coefficient"] = 100 # GW_el   
    # conv_coef.loc[idx[:, :, "Spill", "Water_in"], "coefficient"] = -100 # GW_el # Arbitrary large bypass.
    # conv_coef.loc[idx[:, :, "Spill", "Water_out"], "coefficient"] = 100 # GW_el
    m["Base"].parameter.add(conv_coef, "converter_coefficient")

    # Economic parameters for turbines.
    conv_acc = pd.DataFrame(index=pd.MultiIndex.from_product([["Invest", "OMFix"], ["global"], hydro_techs, hydro_vintage]))
    conv_acc.loc[idx["Invest", :, :, :], "perUnitBuild"] = 2560 # million EUR / unit
    conv_acc.loc[idx["Invest", :, :, :], "useAnnuity"] = 1  # Use annuity method. binary yes/no
    conv_acc.loc[idx["Invest", :, :, :], "amorTime"] = 20  #  Amortization time (years).
    conv_acc.loc[idx["Invest", :, :, :], "interest"] = 0.06  # Interest rate (6%)
    conv_acc.loc[idx["OMFix", :, :, :], "perUnitTotal"] = (
        conv_acc.loc[idx["Invest", :, :, :], "perUnitBuild"] *0.0300
    )  # Mio EUR per unit
    m["Base"].parameter.add(conv_acc, "accounting_converterunits")

    
    # Configure Reservoirs (Storage)
    stor_techs = [ "Lake_Hawea"]#, "Lake_Ohau", "Lake_Pukaki", "Lake_Taupo", "Lake_Tekapo", "Lake_Waikaremoana", "Lakes_Manapouri_Te_Anau"] # ["Hydro_reservoir"]

    # Reservoir technical properties.
    stor_tech = pd.DataFrame(index=pd.MultiIndex.from_product([stor_techs, hydro_vintage]))
    stor_tech.loc[idx[:, :], "lifeTime"] = 100
    stor_tech.loc[idx[:, :], "levelUpperLimit"] = 1 # Normalized fill level limit.
    m["Base"].parameter.add(stor_tech, "storage_techparam")

	# Storage reservoir sizes 
    stor_size = pd.DataFrame(index=pd.MultiIndex.from_product([stor_techs,hydro_vintage, ["Water_in"]]))
    stor_size.loc[idx["Lake_Hawea", :, "Dunstan_in"], "size"] = 1  # GWh_ch / unit  
    stor_size.loc[idx["Lake_Hawea", :, "Dunstan_in"], "selfdischarge"] = 0 # TODO: ask if i can do that to all commodities or only water comming in
    # stor_size = pd.DataFrame(index=pd.MultiIndex.from_product([stor_techs,hydro_vintage, ["Water_in"]]))
    # stor_size.loc[idx["Hydro_reservoir", :, "Water_in"], "size"] = 1  # GWh_ch / unit  
    # stor_size.loc[idx[:, :, "Water_in"], "selfdischarge"] = 0 
    m["Base"].parameter.add(stor_size, "storage_sizeparam")


    # Installed reservoir capacities (in Mm3) from reservoir_limits.csv
    stor_res = pd.DataFrame(index=pd.MultiIndex.from_product([hydro_nodes, hydro_years, stor_techs]))
    stor_res.loc[idx[:, :, :], "unitsUpperLimit"] = 3000  # units 
    #TODO: mdify capacity of the storage
    stor_res.loc[idx[["CAN"], [2000], "Lake_Hawea"], "unitsBuild"] = 1141.95 / 1000000000 # cubic meter (m³) = 1,000,000,000 cubic millimeters (mm³)
    # stor_res.loc[idx[["CAN"], [2000], :], "unitsBuild"] = 2517.2429 # GWh_el
    # stor_res.loc[idx[["HBY"], [2000], :], "unitsBuild"] = 154.2635  # GWh_el
    # stor_res.loc[idx[["OTG"], [2000], :], "unitsBuild"] = 729.5595 # GWh_el
    # stor_res.loc[idx[["WTO"], [2000], :], "unitsBuild"] = 587.1371 # GWh_el

    stor_res.loc[idx[:, :, :], "noExpansion"] = 1
    m["Base"].parameter.add(stor_res, "storage_reservoirparam")

    stor_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], ["global"], stor_techs, hydro_vintage]
        )
    )
    stor_acc.loc[idx["Invest", :, :, :], "perUnitBuild"] = 1650  # million EUR / unit
    stor_acc.loc[idx["Invest", :, :, :], "useAnnuity"] = 1  # binary yes/no
    stor_acc.loc[idx["Invest", :, :, :], "amorTime"] = 20  # years
    stor_acc.loc[idx["Invest", :, :, :], "interest"] = 0.06  # percent/100
    stor_acc.loc[idx["OMFix", :, :, :], "perUnitTotal"] = (
        stor_acc.loc[idx["Invest", :, :, :], "perUnitBuild"] * 0.03
    )  # Mio EUR per unit
    m["Base"].parameter.add(stor_acc, "accounting_storageunits")


# %% [markdown]
# # Hydropower inflow profile creation
import pandas as pd
import numpy as np
from IPython.display import display

# Define path
file_path = "C:/Local/REMix/remix_nz/input/brownfield" 


# %% [markdown]
# ## Access 3 relevant databases
# 1. Inflow (in cumecs) data from JADE for every catchment
inflow_df = pd.read_csv(f"{file_path}/hydro/inflows.csv", skiprows=[0,  3, 4, 5], header=0)
inflow_df.rename(columns={"CATCHMENT": "Year", "Unnamed: 1": "Week"}, inplace=True)
inflow_df.head()




# %% [markdown]
# ## Management of inflow 

# ### Manage inflow in water units for a single year
#  Select a year to work with
year_sel = 2020



# %% [markdown]
# 2. Group catchments in the same region
def group_catchments(df=inflow_df):
    df2=df.copy()
    renamed_df = df2[df2['Year'] == year_sel]
    # define columns names corresponding to catchments in each region
    column_mapping = {
        'BOP': ['Lake_Matahina'],
        'CAN': ['Lake_Aviemore', 'Lake_Benmore', 'Lake_Coleridge', 'Lake_Ohau', 'Lake_Pukaki', 'Lake_Roxburgh', 'Lake_Tekapo', 'Lake_Waitaki'],
        'CEN': ['Lake_Moawhango', 'Lake_Rotoaira', 'Mangahao_head'],
        'HBY': ['Lake_Waikaremoana'],
        'NEL': ['Lake_Cobb'],
        'OTG': ['Lake_Dunstan', 'Lake_Hawea', 'Lake_Wanaka'], #, 'Lakes_Manapouri_Te_Anau'],
        'WTO': ['Lake_Arapuni', 'Lake_Aratiatia', 'Lake_Atiamuri', 'Lake_Karapiro', 'Lake_Maraetai', 'Lake_Ohakuri', 'Lake_Taupo', 'Lake_Waipapa', 'Lake_Whakamaru']}
        
    # Rename columns based on the dictionary
    for new_column, original_columns in column_mapping.items():
        for original_column in original_columns:
            if original_column in renamed_df.columns:
                renamed_df.rename(columns={original_column: new_column}, inplace=True)
    # Group by columns with the same name and sum the values
    grouped_df = renamed_df.groupby(renamed_df.columns, axis=1).sum()
    
    # rearange cols
    cols = grouped_df.columns.tolist()
    new_order = cols[-2:] + cols[:-2]
    new_df = grouped_df[new_order]
    new_df2 =new_df.drop(new_df.columns[1], axis=1)
    new_df2.set_index('Week', inplace=True)
    return new_df2.T
new_df=group_catchments()
display(new_df)


# %% [markdown]
# Inflow: energy units in hourly resolution

def repeat_columns(df, repeat_factor=168, year_used=year_sel):

    # Assuming df is your DataFrame

    # Initialize an empty DataFrame to store the repeated columns
    repeated_df = pd.DataFrame()

    # Iterate over each column in the original DataFrame
    for col in df.columns:
        # Repeat each column 168 times and concatenate to the result
        repeated_df = pd.concat([repeated_df, pd.concat([df[col]] * repeat_factor, axis=1)], axis=1)

    return repeated_df
hourly_inflow_tiny = repeat_columns(new_df)
display(hourly_inflow_tiny)


#%%
def repeat_columns_and_add_average(df, target_columns=8760):
    # Assuming df is your DataFrame

    # Create a new column that is the average of the first and last column
    average_column = (df.iloc[:, 0] + df.iloc[:, -1]) / 2

    # Repeat the entire DataFrame 168 times along the columns axis
    repeated_df = pd.concat([df, pd.concat([average_column] * (target_columns - len(df.columns)), axis=1)], axis=1)

    # Rename the columns with sequential numbers from 1 to 8760
    repeated_df.columns = np.arange(1, target_columns + 1)
    repeated_df = repeated_df.multiply(-1)
    # Specify the path and name of the CSV file
    file_destiny = f"{file_path}/hydro_inflow.csv"
    # Export the DataFrame to a CSV file, with index
    repeated_df.to_csv(file_destiny, index=True)

    return repeated_df

hourly_inflow = repeat_columns_and_add_average(hourly_inflow_tiny)
display(hourly_inflow)
# %%


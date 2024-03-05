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
# 2. Existing powerplant database (to obtain installed capacity, average yearly generation)
existing_gen_df = pd.read_csv(f"{file_path}/power-plant-nz-database.csv",  usecols=[1,7,21,22,23,34] , header=0)

# %% [markdown]
# 3. Weight that every catchment will have in each region
weight_df= pd.read_csv(f"{file_path}/hydro/inflow2energy.csv", header=0)
weight_dict = dict(zip(weight_df['Catchment'], weight_df['Cap_percent']))
# print(weight_dict)



# %% [markdown]
# ## Management of inflow 

# ### Manage inflow in water units for a single year
#  Select a year to work with
year_sel = 2022

# %% [markdown]
# 1. Multiply each catchment times its weight 
def weight_catchment(year_used=year_sel, inflows=inflow_df, dict=weight_dict):
    # Create a copy of the inflow_df dataframe
    weighted_catchment = inflows.copy()
    for column, scalar in dict.items():
        if column in weighted_catchment.columns:
            weighted_catchment[column] = weighted_catchment[column] * scalar
    weighted_catchment_year = weighted_catchment[weighted_catchment['Year'] == year_used]
    return weighted_catchment_year
weighted_catchment = weight_catchment()

# %% [markdown]
# 2. Group catchments in the same region
def group_catchments(df=weighted_catchment):
    renamed_df=df.copy()
    # define columns names corresponding to catchments in each region
    column_mapping = {
        'BOP': ['Lake_Matahina'],
        'CAN': ['Lake_Aviemore', 'Lake_Benmore', 'Lake_Coleridge', 'Lake_Ohau', 'Lake_Pukaki', 'Lake_Roxburgh', 'Lake_Tekapo', 'Lake_Waitaki'],
        'CEN': ['Lake_Moawhango', 'Lake_Rotoaira', 'Mangahao_head'],
        'HBY': ['Lake_Waikaremoana'],
        'NEL': ['Lake_Cobb'],
        'OTG': ['Lake_Dunstan', 'Lake_Hawea', 'Lake_Wanaka', 'Lakes_Manapouri_Te_Anau'],
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

    df_reordered = grouped_df[new_order]

    return df_reordered
grouped_catchments=group_catchments()



# %% [markdown]
# ### Calculate coefficient 'k=E/(Qt)'
# 1. We get E for each region from the dataframe 'existing_hydro'. 
def existing_gen_data():
    # From the existing power plants database we collect information per region on existing capacity, average anual generation and maximum storage

    # Filter rows where "Type" is "Hydroelectric"
    hydro_df = existing_gen_df[(existing_gen_df["Type"] == "Hydroelectric")] 

    # Group by "Node" and sum the specified columns
    grouped_df = hydro_df.groupby("Node").agg({
        "Capacity_MW": "sum",
        "Avg_Ann_Gen_GWh": "sum",
        "Hydro_max_storage_m3": "sum"
    }).reset_index()  # Reset index to make "Node" a regular column
    # Change units of capacity to fit units in REMix
    grouped_df.rename(columns={"Capacity_MW": "Capacity_GW"}, inplace=True)
    grouped_df["Capacity_GW"] /= 1000

    # Add a row at the end with the total sum for each column
    total_row = pd.Series({
        "Node": "Total",
        "Capacity_GW": grouped_df["Capacity_GW"].sum(),
        "Avg_Ann_Gen_GWh": grouped_df["Avg_Ann_Gen_GWh"].sum(),
        "Hydro_max_storage_m3": grouped_df["Hydro_max_storage_m3"].sum()
    })

    # Concatenate the DataFrame with the total row
    grouped_df_total = pd.concat([grouped_df, total_row.to_frame().T], ignore_index=True)
    return grouped_df
existing_hydro=existing_gen_data()

# %% [markdown]
# 2. Value of Q used to get k will be the average Q per region, obtained from the dataframe 'inflow_df'/
#  (first getting the average per catchment and then applying the catchment weight to get the average)/
# and then grouping per region getting a weighted average according to the catchments weighted sum
def get_Q_region():
    # Calculate the column averages and store them in a dictionary to use ltr
    averages_dict = (inflow_df.mean()).to_dict()
    # print(averages_dict)
    # Duplicate catchment info
    avg_Q_df=weight_df.copy()
    # Create a new column by mapping values from the dictionary to the "Catchment" column
    avg_Q_df['Overall_avg_Q'] =avg_Q_df['Catchment'].map(averages_dict)
    # Create a new column "Avg_weighted" by multiplying columns 3 and 4
    avg_Q_df['Q_weighted_cumecs'] =avg_Q_df.iloc[:, 3] * avg_Q_df.iloc[:, 4]    
    summarized_df = avg_Q_df.groupby("Node").sum().iloc[:,[4]]
    summarized_df = summarized_df.reset_index()
    summarized_df.rename(columns={'index': 'Node'}, inplace=True)
    return summarized_df  
avg_Q_df=get_Q_region()
# display(avg_Q_df)

# %% [markdown]
# 3. Value of k used is 'k=E/(Qt)'
def get_k():
    filtered_existing_hydro = pd.merge(existing_hydro, avg_Q_df['Node'], how='inner', on='Node')
    filtered_existing_hydro['Q_weighted_cumecs'] =avg_Q_df['Q_weighted_cumecs'].copy()
    with_k_df=filtered_existing_hydro.copy()
    with_k_df['k=E/(Qt)'] = with_k_df.iloc[:, 2] / ( with_k_df.iloc[:, 4] * 8760 )
    with_k_df= with_k_df.iloc[:, [0,2,3,5]] 
    return with_k_df
coef_k_df = get_k()
# display(coef_k_df)






# %% [markdown]
# ### Convert inflow from water units to potential energy
# The energy will be E=Q*t*k Where Q is the number we have hourly


def multiply_columns_by_scalar(coef_k_df, grouped_catchments):
    # Create a copy of the grouped_catchments DataFrame to avoid modifying the original
    result_df = grouped_catchments.copy()

    # Iterate over the "Nodes" column in coef_k_df
    for node in coef_k_df['Node']:
        # Check if the node is a column in grouped_catchments
        if node in grouped_catchments.columns:
            # Multiply the corresponding column by the scalar in "k=E/(Qt)"
            result_df[node] *= coef_k_df.loc[coef_k_df['Node'] == node, 'k=E/(Qt)'].values[0]
    new_df =result_df.drop(result_df.columns[1], axis=1)
    new_df.set_index('Week', inplace=True)

    # Add a total row
    #new_df.loc['Total'] = new_df.sum()
    return new_df.T

# %% [markdown]
# Inflow: energy units in weekly resolution
weekly_inflow = multiply_columns_by_scalar(coef_k_df, grouped_catchments)
# display(weekly_inflow )


# %% [markdown]
# Inflow: energy units in hourly resolution

def repeat_columns(df, repeat_factor=168):
    # Assuming df is your DataFrame

    # Initialize an empty DataFrame to store the repeated columns
    repeated_df = pd.DataFrame()

    # Iterate over each column in the original DataFrame
    for col in df.columns:
        # Repeat each column 168 times and concatenate to the result
        repeated_df = pd.concat([repeated_df, pd.concat([df[col]] * repeat_factor, axis=1)], axis=1)

    return repeated_df
hourly_inflow_tiny = repeat_columns(weekly_inflow)
# display(hourly_inflow_tiny)

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

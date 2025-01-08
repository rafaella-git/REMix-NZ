# %% [markdown]
# # Hydropower inflow profile creation
import pandas as pd
import numpy as np
from IPython.display import display

# Define path
file_path = "C:/Local/REMix/remix_nz/input/brownfield" 
case_name = "_dividedby1000test"

# %% [markdown]
# ## Access 3 relevant databases
# 1. Inflow (in cumecs) data from JADE for every catchment
inflow_df = pd.read_csv(f"{file_path}/hydro/inflows.csv", skiprows=[0,  3, 4, 5], header=0)
inflow_df.rename(columns={"CATCHMENT": "Year", "Unnamed: 1": "Week"}, inplace=True)
inflow_df.head()
display(inflow_df)



# %% [markdown]
# ## Management of inflow 

# ### Manage inflow in water units for a single year
#  Select a year to work with
year_sel = 2020



# %% [markdown]
# 2. Group catchments
def group_catchments(df=inflow_df):
    df2=df.copy()
    renamed_df = df2[df2['Year'] == year_sel]
    # Divide all numeric values by 1000, to make the unit more proportionate (now its m3/1000)
    renamed_df.iloc[:, 2:] = renamed_df.iloc[:, 2:] / 1000

    # define columns names corresponding to catchments in each region
    column_mapping = {
        # 'BOP': ['Lake_Matahina'],
        # 'CAN': ['Lake_Aviemore', 'Lake_Benmore', 'Lake_Coleridge', 'Lake_Ohau', 'Lake_Pukaki', 'Lake_Roxburgh', 'Lake_Tekapo', 'Lake_Waitaki'],
        # 'CEN': ['Lake_Moawhango', 'Lake_Rotoaira', 'Mangahao_head'],
        # 'HBY': ['Lake_Waikaremoana'],
        # 'NEL': ['Lake_Cobb'],
        # 'OTG': ['Lake_Dunstan', 'Lake_Hawea', 'Lake_Wanaka'], #, 'Lakes_Manapouri_Te_Anau'],
        # 'WTO': ['Lake_Arapuni', 'Lake_Aratiatia', 'Lake_Atiamuri', 'Lake_Karapiro', 'Lake_Maraetai', 'Lake_Ohakuri', 'Lake_Taupo', 'Lake_Waipapa', 'Lake_Whakamaru']}
        'BOP': ['Lake_Matahina'],
        'CAN': ['Lake_Aviemore', 'Lake_Benmore', 'Lake_Coleridge', 'Lake_Ohau', 'Lake_Pukaki', 'Lake_Roxburgh', 'Lake_Tekapo', 'Lake_Waitaki'],
        'CEN': ['Lake_Moawhango', 'Lake_Rotoaira', 'Mangahao_head'],
        'HBY': ['Lake_Waikaremoana'],
        'NEL': ['Lake_Cobb'],
        'OTG': ['Lake_Dunstan', 'Lake_Hawea', 'Lake_Wanaka'],
        'Manapouri_in': ['Lakes_Manapouri_Te_Anau'],
        'WTO': ['Lake_Arapuni', 'Lake_Aratiatia', 'Lake_Atiamuri', 'Lake_Karapiro', 'Lake_Maraetai', 'Lake_Ohakuri', 'Lake_Taupo', 'Lake_Waipapa', 'Lake_Whakamaru']}
        #        
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
    # Create a new column that is the average of the first and last column
    average_column = (df.iloc[:, 0] + df.iloc[:, -1]) / 2

    # Repeat the entire DataFrame 168 times along the columns axis
    repeated_df = pd.concat([df, pd.concat([average_column] * (target_columns - len(df.columns)), axis=1)], axis=1)

    # Rename the columns with sequential numbers from 1 to 8760
    repeated_df.columns = np.arange(1, target_columns + 1)
    repeated_df = repeated_df.multiply(-1)
    # Specify the path and name of the CSV file
    file_destiny = f"{file_path}/hydro_inflow{case_name}.csv"
    # Export the DataFrame to a CSV file, with index
    # repeated_df.to_csv(file_destiny, index=True)
    return repeated_df

hourly_inflow = repeat_columns_and_add_average(hourly_inflow_tiny)
display(hourly_inflow)
# %%


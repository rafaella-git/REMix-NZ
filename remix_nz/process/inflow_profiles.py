# %% [markdown]
# # Hydropower inflow profile creation
import pandas as pd
import numpy as np
import os
from IPython.display import display

# Define path
file_path = "C:/Local/REMix/remix_nz/input/brownfield" 



# %% [markdown]
# ## Access 3 relevant databases
# 1. Inflow (in cumecs) data from JADE for every catchment
inflow_df = pd.read_csv(f"{file_path}/hydro/inflows.csv", skiprows=[0,  3, 4, 5], header=0)
inflow_df.rename(columns={"CATCHMENT": "Year", "Unnamed: 1": "Week"}, inplace=True)


# %% [markdown]
# 2. Existing powerplant database (to obtain installed capacity, average yearly generation)
existing_gen_df = pd.read_csv(f"{file_path}/power-plant-nz-database.csv",  usecols=[1,7,21,22,23,34] , header=0)

# %% [markdown]
# 3. Weight that every catchment will have in each region
weight_df= pd.read_csv(f"{file_path}/hydro/inflow2energy.csv", header=0)
weight_dict = dict(zip(weight_df['Catchment'], weight_df['Cap_percent']))
# print(weight_dict)




# %% [markdown]
# ## Calculate coefficient 'k=E/(Qt)'
# 1. We get E (avg energy generation) for each region from the dataframe 'existing_hydro'. 
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
# ## Management of inflow 

# ### Manage inflow in water units for all the years
# 1. Multiply each catchment times its weight
# 2. Group weighted catchments in the same region
# 3. Return df with cumecs per region, for every year

def group_weighted_catchments(inflows=inflow_df, dict=weight_dict):
    # 1. Multiply each catchment times its weight 
    weighted_catchment = inflows.copy()
    for column, scalar in dict.items():
        if column in weighted_catchment.columns:
            weighted_catchment[column] = weighted_catchment[column] * scalar
    weighted_catchment_year = weighted_catchment

    # 2. Group weighted catchments in the same region
    renamed_df=weighted_catchment_year.copy()
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
    result_df = grouped_df[new_order]

    # 3. Return df with cumecs per region, for every year
    return result_df 
display(group_weighted_catchments())


# %% [markdown]
# ## Select relevant years
# Driest
# Most average

def analyze_years(df, num_years):
    # Calculate yearly sum and set "Year" as index
    sum_df = df[df.columns[2:]].groupby(df['Year']).sum().reset_index(drop=False)

    # Add a column named "Total" with the sum of all the remaining columns
    sum_df['Sum'] = sum_df.iloc[:, 1:].sum(axis=1)  # Exclude "Year" column

    # Display the result 
    #display(sum_df)

    new_df=sum_df[["Year","Sum"]]

    #display(new_df)

    # Remove the row where the year is 1932
    new_df = new_df[new_df['Year'] != 1932]
    
    # Convert Year values to integers
    new_df['Year'] = new_df['Year'].astype(int)
    # Find the specified number of years with the lowest values in the "Sum" column
    lowest_sum_years = new_df.nsmallest(num_years, 'Sum')['Year'].tolist()

    # Calculate the average of the "Sum" column
    average_sum = new_df['Sum'].mean()

    # Find the specified number of years closest to the average "Sum" across all years
    closest_to_average_years = new_df.iloc[(new_df['Sum'] - average_sum).abs().argsort()[:num_years]]['Year'].tolist()

    return lowest_sum_years, closest_to_average_years

num_years = 7
driest_years, average_years = analyze_years(group_weighted_catchments(), num_years)

print(f"{num_years} years with the lowest values:", driest_years)
print(f"{num_years} years closest to the average:", average_years)








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

def repeat_columns(df, repeat_factor=168):
    # Assuming df is your DataFrame

    # Initialize an empty DataFrame to store the repeated columns
    repeated_df = pd.DataFrame()

    # Iterate over each column in the original DataFrame
    for col in df.columns:
        # Repeat each column 168 times and concatenate to the result
        repeated_df = pd.concat([repeated_df, pd.concat([df[col]] * repeat_factor, axis=1)], axis=1)

    return repeated_df

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





# %% [markdown]
# Ultimate dataframe depending on year selected
def yearly_inflow_df(year_sel):
    doc_name=f"inflow{year_sel}"
    # Get the weekly inflow where "Year" is year_sel}
    historic_df=group_weighted_catchments()
    df= historic_df[(historic_df["Year"] == year_sel)] 
    # Inflow: energy units in weekly resolution
    weekly_inflow = multiply_columns_by_scalar(coef_k_df, df)
    # Inflow: energy units in hourly resolution
    hourly_inflow_tiny = repeat_columns(weekly_inflow)
    hourly_inflow = repeat_columns_and_add_average(hourly_inflow_tiny)
    hourly_inflow.insert(0, 'year', year_sel)
    hourly_inflow.insert(1, 'sector', 'Inflow')
    hourly_inflow.insert(2, 'carrier', 'Water_in')
    display(hourly_inflow)
    hourly_inflow.to_csv(f"{file_path}/hydro/inflow_profile/{doc_name}.csv", index=True)
    return hourly_inflow

yearly_inflow_df(2020)


# %% [markdown]
def inflow_csv(years):
    # years is either driest_years or average_years
    doc_name="average_years"
    if years == driest_years:
        doc_name="driest_years"
    # keep a list of the resulting dataframes
    dfs = []
    for year in years:
        df = yearly_inflow_df(year)
        # Append the generated DataFrame to the list
        dfs.append(df)
    # concatenate all DataFrames in the list into one DataFrame
    concatenated_df = pd.concat(dfs, ignore_index=False)

    # List of values to repeat
    values_to_repeat = [2020, 2025, 2030, 2035, 2040, 2045, 2050]

    # Repeat the list 7 times
    repeated_values = np.repeat(values_to_repeat, 7)

    # Make sure the length of repeated values matches the length of the DataFrame
    repeated_values = repeated_values[:len(concatenated_df)]

    # Replace the values in the "Year" column with the repeated values
    concatenated_df['year'] = repeated_values

    concatenated_df.to_csv(f"{file_path}/hydro/inflow_profile/{doc_name}.csv", index=True)
    return concatenated_df

inflow_csv(driest_years)
inflow_csv(average_years)


# %%

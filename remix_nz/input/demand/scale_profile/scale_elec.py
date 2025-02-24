import pandas as pd
import numpy as np


file_path = 'original_elec.csv'  
df = pd.read_csv(file_path)

# Example dictionary for total yearly demand (in GWh)
total_yearly_demand = {
    2020: 39725.56483, 
    2030: 57778.20452,  
    2040: 400,   
    2050: 500   
}



# Function to scale the demand values to match the total yearly demand
def scale_demand(df, total_yearly_demand):
    # Create a copy of the DataFrame to avoid modifying the original
    df_scaled = df.copy()
    
    # Iterate over each unique year in the DataFrame
    for year in df_scaled['year'].unique():
        # Filter the DataFrame for the current year
        year_mask = df_scaled['year'] == year
        year_data = df_scaled[year_mask]
        
        # Extract the demand columns (time-series data)
        demand_columns = [col for col in year_data.columns if col.startswith('t')]
        demand_values = year_data[demand_columns].values
        
        # Calculate the sum of the original demand values for the year
        original_sum = np.sum(demand_values)
        
        # Get the total yearly demand for the current year
        total_demand = total_yearly_demand.get(year, 0)  # Default to 0 if year not in dictionary
        
        # Compute the scaling factor
        scaling_factor = total_demand / original_sum
        
        # Scale the demand values
        scaled_values = demand_values * scaling_factor
        
        # Update the DataFrame with the scaled values
        df_scaled.loc[year_mask, demand_columns] = scaled_values
    
    return df_scaled

# Apply the scaling function to the DataFrame
df_scaled = scale_demand(df, total_yearly_demand)

# Save the scaled DataFrame to a new CSV file (optional)
df_scaled.to_csv('scaled_nz-elec.csv', index=False)

# # Example demonstration with simplified data
# # Create a small example DataFrame
# example_data = {
#     'node': ['A', 'A', 'B', 'B'],
#     'year': [2020, 2020, 2030, 2030],
#     'sector': ['Wholesale', 'Wholesale', 'Retail', 'Retail'],
#     'carrier': ['Electricity', 'Electricity', 'Electricity', 'Electricity'],
#     't0001': [10, 20, 30, 40],
#     't0002': [15, 25, 35, 45],
#     't0003': [20, 30, 40, 50]
# }
# example_df = pd.DataFrame(example_data)

# # Example total yearly demand for the simplified data
# example_total_yearly_demand = {
#     2020: 100,  # 100 GW for 2020
#     2030: 200   # 200 GW for 2030
# }

# # Apply the scaling function to the example DataFrame
# example_scaled = scale_demand(example_df, example_total_yearly_demand)

# # Display the example DataFrame before and after scaling
# print("Original Example DataFrame:")
# print(example_df)
# print("\nScaled Example DataFrame:")
# print(example_scaled)
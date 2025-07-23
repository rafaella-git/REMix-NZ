import pandas as pd


# File path
csv_path = "C:/Local/REMix/remix_nz/process/data-for-paper/Toktarova2019-ElectricityProjection-NZ.csv"

# Read the CSV without headers
df_raw = pd.read_csv(csv_path, header=None)

# Extract the years from the first row (index 0)
years = df_raw.iloc[0, 1:].tolist()  # Skip the first column which says 'year'

# Extract hourly data from row 6 (index 5) onwards
hourly_data = df_raw.iloc[5:, 1:].copy()  # Skip first column (e.g., "Hour_1 in MW")

# Clean values: remove commas and convert to float
hourly_data = hourly_data.applymap(lambda x: float(str(x).replace(",", "").replace('"', '')))

# Clean and convert values
hourly_data = hourly_data.applymap(lambda x: float(str(x).replace(",", "").replace('"', '')))

# Replace NaNs with 0
hourly_data = hourly_data.fillna(0)

# Reset index to hours in REMix format
hourly_data.index = [f't{i+1:04d}' for i in range(hourly_data.shape[0])]


# Assign column names (years)
hourly_data.columns = years

# Final DataFrame: each column = year, each row = hourly value
df_hourly = hourly_data


# print(df_hourly.head())
# print(df_hourly.tail())

# Normalize each column by its yearly sum
df_normalized = df_hourly.div(df_hourly.sum(axis=0), axis=1)
df_normalized.columns = df_normalized.columns.astype(int)
df_normalized.to_csv('C:/Local/REMix/remix_nz/process/data-for-paper/Normalized-Load-Profile.csv', index=True)


#pending, fix name of resulting csv
# fixlenght, its giving less that 8760 hours
#fix name in first row
import pandas as pd
import numpy as np

location = location = "C:\\Local\\REMix\\remix_nz\\input\\brownfield\\hydro\\" # Adjust as needed

# Read the CSV file, keeping all columns
df = pd.read_csv(f"{location}inflow_week.csv")

# Separate the first 2 columns and the remaining columns
first_2_cols = df.iloc[:, :2]
remaining_cols = df.iloc[:, 2:]

# Repeat the remaining columns 168 times
repeated_cols = pd.DataFrame(np.tile(remaining_cols.values, (1, 168)), columns=remaining_cols.columns.tolist() * 168)

# Repeat the last column 4 more times
last_col_repeated = pd.DataFrame(np.tile(remaining_cols.iloc[:, -1:], (1, 4)), columns=[f"{remaining_cols.columns[-1]}_extra{i+1}" for i in range(4)])

# Combine all columns
result_df = pd.concat([first_2_cols, repeated_cols, last_col_repeated], axis=1)

# Save the result to a CSV file
result_df.to_csv(f"{location}repeated_columns.csv", index=False)

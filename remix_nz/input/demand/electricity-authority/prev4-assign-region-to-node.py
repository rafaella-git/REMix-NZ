# This code takes the info from previously created annual hourly demand (8760 hours for each of the nodes) for 2019 
# It verifies that the original data has the correct amount of trading periods, if it doesnt, it fills them to then unify all months in a yearly file

# There is 1 main outcome:
# 1. an annual hourly demand (8760 hours for each of the nodes) for 2019
# 

import pandas as pd

# --- File paths --- doe
input_path = "C:/Local/REMix/remix_nz/input/demand/electricity-authority/"
demand_path = f"{input_path}hourly_demand_all_2019_with_coords.csv"
mapping_path = "C:/Local/REMix/remix_nz/input/shapefiles/match-sites-16regions.csv"
output_path = f"{input_path}hourly_demand_all_2019_with_regions.csv"

# --- Load datasets ---
demand_df = pd.read_csv(demand_path, parse_dates=["Period start", "Period end"])
mapping_df = pd.read_csv(mapping_path)

# --- Clean and align keys ---
mapping_df.columns = mapping_df.columns.str.strip()
mapping_df['Site'] = mapping_df['Site'].str.strip()

# --- Merge on Region Code ↔ Site ---
enriched_df = demand_df.merge(
    mapping_df,
    how="left",
    left_on="Region Code",
    right_on="Site"
)

# --- Drop 'Site' helper column if not needed ---
enriched_df.drop(columns=["Site"], inplace=True)

# --- Check for unmatched region codes ---
unmatched = enriched_df[
    mapping_df.columns.difference(['Site']).tolist()
].isna().all(axis=1)

if unmatched.any():
    unmatched_codes = enriched_df.loc[unmatched, "Region Code"].drop_duplicates()
    print("⚠️ Region Codes with no sectoral mapping:")
    print(unmatched_codes.to_string(index=False))
else:
    print("✅ All Region Codes matched with sectoral classifications.")

# --- Save final enriched file ---
enriched_df.to_csv(output_path, index=False)
print(f"\n📁 Final file with sectoral groups saved to:\n{output_path}")

# --- Preview ---
print("\n📋 Preview of new sectoral columns added:")
print(enriched_df[
    ['Region Code', 'Region Name'] + mapping_df.columns.difference(['Site']).tolist()
].drop_duplicates().head(10))

print(enriched_df)
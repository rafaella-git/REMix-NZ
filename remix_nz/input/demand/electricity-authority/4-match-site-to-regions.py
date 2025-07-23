# This code takes the info from previously created yearly hourly demand for each node, a verified and patched version with 8760 hours per node
# It matches the node data to regions, there is a 16 region option (admin regions, originally from pypsa's data), and the 11 regions from remix 

# There is 1 main outcome:
# 1. an annual hourly demand (8760 hours for each of the nodes) for 2019 WITH ADDED region assignment

import pandas as pd

# --- 1. File paths ---
input_path = "C:/Local/REMix/remix_nz/input/demand/electricity-authority/"
demand_path = f"{input_path}hourly_demand_all_2019_cleaned.csv"
assignment_path = f"{input_path}site_region_assignment_named.csv"
output_path = f"{input_path}hourly_demand_all_2019_with_regions.csv"

# --- 2. Load data ---
print("📥 Loading data...")
demand_df = pd.read_csv(demand_path, parse_dates=["Period start", "Period end"])
assignment_df = pd.read_csv(assignment_path)

print("\n🔎 Demand data preview:")
print(demand_df.head())

print("\n🔎 Site-region assignment data preview:")
print(assignment_df.head())

# --- 3. Prepare region assignment data ---
assignment_lookup = assignment_df[["MXLOCATION", "Administrative Region", "id"]].copy()
assignment_lookup.rename(columns={
    "Administrative Region": "16Regions",
    "id": "11Regions"
}, inplace=True)

# --- 4. Apply manual patch to Region Code where no match is expected ---
manual_region_site_map = {
    "MKE": "MKT",  # McKee → McKee Tee
    "NPL": "JRD"   # New Plymouth → Junction Road
}

demand_df["MXLOCATION_final"] = demand_df["Region Code"].map(manual_region_site_map)
demand_df["MXLOCATION_final"] = demand_df["MXLOCATION_final"].fillna(demand_df["Region Code"])

# --- 5. Merge demand with site-region assignment using patched values ---
merged_df = demand_df.merge(
    assignment_lookup,
    how="left",
    left_on="MXLOCATION_final",
    right_on="MXLOCATION"
)

# --- 5a. Add Island column (NI or SI) based on 16Regions ---
north_island_regions = [
    "Auckland", "Bay of Plenty", "Gisborne", "Hawkes Bay",
    "Manawatu-Whanganui", "Northland", "Taranaki", "Waikato", "Wellington"
]
south_island_regions = [
    "Canterbury", "Marlborough", "Nelson", "Otago",
    "Southland", "Tasman", "West Coast"
]

merged_df["Island"] = merged_df["16Regions"].apply(
    lambda x: "NI" if x in north_island_regions else ("SI" if x in south_island_regions else None)
)

# --- 6. Show merged result ---
print("\n🔗 Merged demand data (with region names):")
print(merged_df[["Region Code", "MXLOCATION_final", "16Regions", "11Regions"]].drop_duplicates().head(10))

# --- 7. Report unmatched Region Codes ---
unmatched = merged_df[merged_df["16Regions"].isna() | merged_df["11Regions"].isna()]
if not unmatched.empty:
    print("\n⚠️ Warning: Some region codes could not be matched to sites:")
    print(unmatched[["Region Code", "Region Name"]].drop_duplicates().to_string(index=False))
else:
    print("\n✅ All Region Codes successfully mapped to both 16 and 11 region names.")

# --- 8. Save result ---
merged_df.to_csv(output_path, index=False)
print(f"\n📁 Final demand file with region names exported to:\n{output_path}")

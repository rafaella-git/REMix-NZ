# This code takes the info from previously created hourly demand for each month and each node (originally processed from Electricity Authority - EMI (market statistics and tools) https://www.emi.ea.govt.nz/Wholesale/Reports/R_NSPL_DR?_si=v|3)
# It verifies that the original data has the correct amount of trading periods, if it doesnt, it fills them to then unify all months in a yearly file

# There is 1 main outcome:
# 1. an annual hourly demand (8760 hours for each of the nodes) for 2019

import pandas as pd



# Load site data
input_path = "C:/Local/REMix/remix_nz/input/demand/electricity-authority/"
sites_path = f"{input_path}sites-transpower.csv"
demand_path = f"{input_path}hourly_demand_all_2019_cleaned.csv"
output_path = f"{input_path}hourly_demand_all_2019_with_coords.csv"

# --- Load Data ---
sites_df = pd.read_csv(sites_path)
demand_df = pd.read_csv(demand_path, parse_dates=["Period start", "Period end"])

# Rename for clarity
sites_df.rename(columns={"X": "Easting", "Y": "Northing"}, inplace=True)

# Extract unique region metadata
region_meta = demand_df[["Region ID", "Region", "Region Code", "Region Name"]].drop_duplicates()

# --- Initial Merge on Region Code ↔ MXLOCATION ---
merged = region_meta.merge(sites_df, how="left", left_on="Region Code", right_on="MXLOCATION")

# --- Check for unmatched ---
unmatched = merged[merged["Easting"].isna()]
manual_region_site_map = {
    "MKE": "MKT",  # McKee → McKee Tee
    "NPL": "JRD"   # New Plymouth → Junction Road
}

if not unmatched.empty:
    print(f"\n🔁 Attempting manual fix for {len(unmatched)} unmatched regions...\n")

    # Apply patch column
    region_meta["MXLOCATION_patch"] = region_meta["Region Code"].map(manual_region_site_map)
    region_meta["MXLOCATION_final"] = region_meta["MXLOCATION_patch"].fillna(region_meta["Region Code"])

    # Merge using patched MXLOCATION
    merged = region_meta.merge(sites_df, how="left", left_on="MXLOCATION_final", right_on="MXLOCATION")

    # Report applied patches
    applied_patches = region_meta[region_meta["Region Code"].isin(manual_region_site_map.keys())]
    print("🛠️ Manual patches applied:")
    for _, row in applied_patches.iterrows():
        print(f"  - Region Code {row['Region Code']} mapped to site '{manual_region_site_map[row['Region Code']]}'")

# --- Final Check: Regions still unmatched ---
still_unmatched = merged[merged["Easting"].isna()]
if not still_unmatched.empty:
    print("\n⚠️ Regions with no matching site even after patch:")
    print(still_unmatched[["Region ID", "Region Code", "Region Name"]].to_string(index=False))
else:
    print("✅ All regions now have a site match.")

# --- Sites never used ---
used_sites = merged["MXLOCATION"].dropna().unique()
unused_sites = sites_df[~sites_df["MXLOCATION"].isin(used_sites)]
if not unused_sites.empty:
    print("\n⚠️ Sites not matched to any region:")
    print(unused_sites[["MXLOCATION", "description"]].to_string(index=False))
else:
    print("✅ All sites are used in region matching.")

# --- Merge into full demand dataframe ---
demand_df = demand_df.merge(
    merged[["Region ID", "Easting", "Northing"]],
    on="Region ID",
    how="left"
)

# --- Save Output ---
demand_df.to_csv(output_path, index=False)
print(f"\n📁 Final demand file with coordinates exported to:\n{output_path}")

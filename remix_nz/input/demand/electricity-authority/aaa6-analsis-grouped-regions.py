# --- Profile Matching and Sectoral Aggregation Script ---

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import normalize

# --- 1. File Paths ---
input_path = "C:/Local/REMix/remix_nz/input/demand/electricity-authority/"
demand_path = f"{input_path}hourly_demand_all_2019_with_regions.csv"
match_results_path = "node_profile_matches.csv"

# --- 2. Load Data ---
print("📥 Loading data...")
demand_df = pd.read_csv(demand_path, parse_dates=["Period start", "Period end"])
match_df = pd.read_csv(match_results_path)

# --- 3. Verify Unique Region Codes ---
print("\n🔍 Unique Region Codes:")
print(f"Demand data: {demand_df['Region Code'].nunique()} unique nodes")
print(f"Match data: {match_df['Region Code'].nunique()} unique nodes")

# --- 4. Assign Sector Labels ---
def classify_sector(profile):
    if "ResRegion Codeential" in profile:
        return "ResRegion Codeential"
    elif "Heavy Ind." in profile:
        return "Industrial"
    elif "Ind & Comm." in profile:
        return "Commercial"
    return "Unclassified"

match_df["Sector"] = match_df["Best Match Profile"].fillna("Unclassified").apply(classify_sector)

# Merge sector info into demand dataframe
demand_df = demand_df.merge(match_df[["Region Code", "Sector"]], on="Region Code", how="left")
demand_df["Sector"] = demand_df["Sector"].fillna("Unclassified")

# --- 5. Aggregate by Region and Sector ---
print("\n🔄 Aggregating demand by region and sector...")

# For 16 Regions
agg_16 = (
    demand_df.groupby(["16Regions", "Sector", "Period start"])
    .agg({"Demand (GWh)": "sum"})
    .reset_index()
)

# For 11 Regions
agg_11 = (
    demand_df.groupby(["Region Code", "Sector", "Period start"])
    .agg({"Demand (GWh)": "sum"})
    .reset_index()
)

print("\n📊 Sample of Aggregated Data (16 Regions):")
print(agg_16.head())

print("\n📊 Sample of Aggregated Data (11 Regions):")
print(agg_11.head())

# --- 6. Plot Total Demand by Sector per Region ---
print("\n📈 Plotting total demand by sector per region...")

# Plotting total yearly demand per sector (16 Regions)
total_16 = (
    agg_16.groupby(["16Regions", "Sector"])
    .agg({"Demand (GWh)": "sum"})
    .reset_index()
    .pivot(index="16Regions", columns="Sector", values="Demand (GWh)")
    .fillna(0)
)

total_16.plot(kind="bar", stacked=True, figsize=(14, 6))
plt.title("Total Yearly Demand by Sector per Administrative Region (16 Regions)")
plt.ylabel("Total GWh")
plt.xlabel("Region")
plt.tight_layout()
plt.show()

# Plotting total yearly demand per sector (11 Regions)
total_11 = (
    agg_11.groupby(["Region Code", "Sector"])
    .agg({"Demand (GWh)": "sum"})
    .reset_index()
    .pivot(index="Region Code", columns="Sector", values="Demand (GWh)")
    .fillna(0)
)

total_11.plot(kind="bar", stacked=True, figsize=(14, 6))
plt.title("Total Yearly Demand by Sector per Region (11 Regions)")
plt.ylabel("Total GWh")
plt.xlabel("Region Code")
plt.tight_layout()
plt.show()


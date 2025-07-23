import pandas as pd

# --- Parameters ---
csv_path = "C:/Local/REMix/remix_nz/input/demand/electricity-authority/hourly_demand_all_2019_with_regions.csv"
region_level = "16Regions"  # ← Change to "11Regions", "Region Code", etc.
output_path = f"aggregated_demand_by_{region_level}.csv"

# --- Load original dataset ---
print("📥 Loading hourly demand data...")
df = pd.read_csv(csv_path, parse_dates=["Period start"])

# --- Check and reassign TUI ---
if "TUI" in df["Region Code"].unique():
    tui_rows = df[df["Region Code"] == "TUI"]
    tui_total_demand = tui_rows["Demand (GWh)"].sum()
    previous_region = tui_rows[region_level].unique()

    print(f"\n🔄 TUI node detected.")
    print(f"   • Current region assignment(s): {previous_region}")
    print(f"   • Total yearly demand for TUI: {tui_total_demand:.6f} GWh")
    print("   🛠️ Reassigning TUI to 'Gisborne' for 16Regions.")

    df.loc[df["Region Code"] == "TUI", region_level] = "Gisborne"
else:
    print("ℹ️ No TUI node found in data.")

# --- Validate expected columns ---
expected_cols = ["Hour of Year", "Demand (GWh)", region_level]
for col in expected_cols:
    if col not in df.columns:
        raise ValueError(f"❌ Missing expected column: '{col}'")

# --- Group and aggregate demand ---
print(f"\n🔄 Aggregating demand by '{region_level}' and Hour of Year...")
agg_df = (
    df.groupby(["Hour of Year", region_level])["Demand (GWh)"]
    .sum()
    .unstack(fill_value=0)
    .reset_index()
)

# --- Clean up and round ---
agg_df = agg_df.rename(columns={"Hour of Year": "hour"})
agg_df.set_index("hour", inplace=True)
agg_df = agg_df.round(6)

# --- Save result ---
agg_df.to_csv(output_path)
print(f"\n✅ Aggregated demand saved to: {output_path}")
print(f"📊 Shape of result: {agg_df.shape} (hours × regions)")

# --- Total demand from aggregated ---
agg_total = agg_df.sum().sum()

# --- Total demand from original ---
original_total = round(df["Demand (GWh)"].sum(), 6)

# --- Print summary ---
print(f"\n⚡ Total energy demand from aggregated data: {agg_total:,.6f} GWh")
print(f"⚡ Total energy demand from original data:   {original_total:,.6f} GWh")

if abs(agg_total - original_total) < 1e-6:
    print("✅ Totals match exactly.")
else:
    print("⚠️ Totals do NOT match! Check for missing or duplicate data.")

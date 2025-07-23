import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

import pandas as pd
import matplotlib.pyplot as plt

# --- Load data ---
demand_path = "hourly_demand_all_2019_with_regions.csv"
match_path = "node_profile_matches.csv"
profiles_path = "normalized_emi_profiles.csv"

# --- Load ---
demand_df = pd.read_csv(demand_path, parse_dates=["Period start"])
match_df = pd.read_csv(match_path)
emi_profiles = pd.read_csv(profiles_path)

# --- Define filters ---
island = "NI"
profile_label = "NI Residential"
opacity = 0.01

# --- Get nodes with this best-match profile ---
selected_nodes = match_df[
    (match_df["Island"] == island) & 
    (match_df["Best Match Profile"] == profile_label)
]["Region Code"].values

print(f"🎯 {len(selected_nodes)} nodes matched as '{profile_label}'")

# --- Prepare hourly matrix ---
pivot = demand_df[
    demand_df["Region Code"].isin(selected_nodes)
].pivot_table(index="Region Code", columns="Period start", values="Demand (GWh)")

# --- Normalize each node ---
pivot = pivot.dropna()
pivot_norm = ((pivot.T - pivot.T.min()) / (pivot.T.max() - pivot.T.min())).T

# --- Define summer and winter windows ---
# These can be tuned — below is generic NZ approximation
summer_start = pd.to_datetime("2019-01-15 00:00")
summer_end   = pd.to_datetime("2019-01-15 23:30")

winter_start = pd.to_datetime("2019-07-15 00:00")
winter_end   = pd.to_datetime("2019-07-15 23:30")

summer_cols = [c for c in pivot_norm.columns if summer_start <= c <= summer_end]
winter_cols = [c for c in pivot_norm.columns if winter_start <= c <= winter_end]

# --- Plotting ---
def plot_day_window(day_cols, title, emi_series):
    plt.figure(figsize=(12, 4))
    for node in pivot_norm.index:
        plt.plot(
            range(len(day_cols)),
            pivot_norm.loc[node, day_cols].values,
            color="gray",
            alpha=opacity
        )
    # Overlay EMI reference profile
    plt.plot(range(len(day_cols)), emi_series, color="black", linewidth=2, label="EMI Reference")
    plt.title(title)
    plt.xlabel("Hour")
    plt.ylabel("Normalized Demand")
    plt.legend()
    plt.tight_layout()
    plt.show()

# Extract EMI reference shape for summer & winter
emi_profile = emi_profiles[profile_label].values
emi_summer = emi_profile[14 * 24 : 15 * 24]  # Jan 15
emi_winter = emi_profile[195 * 24 : 196 * 24]  # Jul 15

# --- Plot Summer ---
plot_day_window(summer_cols, f"Summer Day - {profile_label}", emi_summer)

# --- Plot Winter ---
plot_day_window(winter_cols, f"Winter Day - {profile_label}", emi_winter)

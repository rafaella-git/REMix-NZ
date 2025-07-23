# import pandas as pd
# import matplotlib.pyplot as plt
# from sklearn.cluster import KMeans
# from datetime import datetime

# # === SETTINGS ===
# input_path = "C:/Local/REMix/remix_nz/input/demand/electricity-authority/"
# file_path = f"{input_path}hourly_demand_all_2019_with_coords.csv"
# output_path = f"{input_path}hourly_demand_all_2019_with_inferred_sectors.csv"
# n_clusters = 3

# # === Load data ===
# print("📂 Loading hourly demand data...")
# df = pd.read_csv(file_path, parse_dates=["Period start"])
# print(f"✅ Loaded {len(df):,} rows from {df['Region ID'].nunique()} regions.\n")

# # === Extract time + season ===
# df["Hour"] = df["Period start"].dt.hour
# df["Weekday"] = df["Period start"].dt.dayofweek
# df["IsWeekend"] = df["Weekday"] >= 5
# df["Month"] = df["Period start"].dt.month
# df["Season"] = df["Month"].map({
#     12: "Summer", 1: "Summer", 2: "Summer",
#     6: "Winter", 7: "Winter", 8: "Winter"
# }).fillna("Other")

# # === Aggregate average hourly demand ===
# def average_profile_by(df, condition):
#     return df[condition].groupby(["Region ID", "Hour"])["Demand (GWh)"].mean().unstack()

# weekday_df = df[df["IsWeekend"] == False]
# weekend_df = df[df["IsWeekend"] == True]

# profiles = {
#     "wd_all": average_profile_by(weekday_df, weekday_df["Season"].isin(["Summer", "Winter", "Other"])),
#     "we_all": average_profile_by(weekend_df, weekend_df["Season"].isin(["Summer", "Winter", "Other"])),
#     "wd_winter": average_profile_by(weekday_df, weekday_df["Season"] == "Winter"),
#     "we_winter": average_profile_by(weekend_df, weekend_df["Season"] == "Winter"),
#     "wd_summer": average_profile_by(weekday_df, weekday_df["Season"] == "Summer"),
#     "we_summer": average_profile_by(weekend_df, weekend_df["Season"] == "Summer")
# }

# # === Filter low-demand profiles ===
# valid_regions = (profiles["wd_all"].sum(axis=1) > 0.01) & (profiles["we_all"].sum(axis=1) > 0.01)
# for key in profiles:
#     profiles[key] = profiles[key][valid_regions]

# print(f"✅ Filtered to {len(valid_regions)} valid regions with meaningful demand.\n")

# # === Normalize (shape-only) ===
# weekday_norm = profiles["wd_all"].div(profiles["wd_all"].sum(axis=1), axis=0)
# weekend_norm = profiles["we_all"].div(profiles["we_all"].sum(axis=1), axis=0)

# combined_profile = pd.concat([
#     weekday_norm.add_prefix("wd_"),
#     weekend_norm.add_prefix("we_")
# ], axis=1)

# # Drop profiles with too many NaNs
# min_valid_hours = 20
# valid_profiles = (weekday_norm.notna().sum(axis=1) >= min_valid_hours) & \
#                  (weekend_norm.notna().sum(axis=1) >= min_valid_hours)
# combined_profile = combined_profile.loc[valid_profiles]

# # === KMeans Clustering ===
# print(f"🔍 Running KMeans with {n_clusters} clusters...")
# kmeans = KMeans(n_clusters=n_clusters, random_state=42)
# combined_profile_filled = combined_profile.fillna(0)
# combined_profile["Cluster"] = kmeans.fit_predict(combined_profile_filled)

# # === Manually assign sector labels ===
# cluster_to_sector = {
#     0: "Industrial",
#     1: "Residential",
#     2: "Commercial"
# }
# combined_profile["Inferred Sector"] = combined_profile["Cluster"].map(cluster_to_sector)

# # === Plot with winter/summer curves ===
# fig, axes = plt.subplots(1, n_clusters, figsize=(7 * n_clusters, 5), sharey=True)

# for cluster_id, ax in zip(sorted(combined_profile["Cluster"].unique()), axes):
#     sector = cluster_to_sector[cluster_id]
#     region_ids = combined_profile[combined_profile["Cluster"] == cluster_id].index

#     def plot_avg(season_key, label, style):
#         avg = profiles[season_key].loc[region_ids].mean()
#         ax.plot(range(24), avg.values, label=label, **style)

#     plot_avg("wd_winter", "Weekday Winter", {"color": "blue", "linestyle": "-", "marker": "o"})
#     plot_avg("wd_summer", "Weekday Summer", {"color": "blue", "linestyle": "--", "marker": "o"})
#     plot_avg("we_winter", "Weekend Winter", {"color": "orange", "linestyle": "-", "marker": "x"})
#     plot_avg("we_summer", "Weekend Summer", {"color": "orange", "linestyle": "--", "marker": "x"})

#     ax.set_title(f"{sector} ({len(region_ids)} regions)")
#     ax.set_xlabel("Hour")
#     ax.set_ylabel("Normalized Demand")
#     ax.grid(True)
#     ax.legend()

# plt.suptitle("📊 Clustered Load Shapes by Sector with Seasonal Curves", fontsize=16)
# plt.tight_layout()
# plt.show()

# # === Merge sector labels into full dataset ===
# df["Inferred Sector"] = df["Region ID"].map(combined_profile["Inferred Sector"])
# df["Cluster"] = df["Region ID"].map(combined_profile["Cluster"])
# df["Excluded"] = df["Inferred Sector"].isna()

# # === Save output ===
# df.to_csv(output_path, index=False)

# # === Summary ===
# print("\n📋 Sector distribution (valid regions):")
# print(df[~df["Excluded"]][["Region ID", "Inferred Sector"]].drop_duplicates()["Inferred Sector"].value_counts())

# print("\n⚠️ Excluded (unclustered) regions:", df["Excluded"].sum())
# print(f"\n💾 Final dataset saved to:\n{output_path}")

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans

# === Setup ===
file_path = "C:/Local/REMix/remix_nz/input/demand/electricity-authority/hourly_demand_all_2019_with_coords.csv"
output_dir = "C:/Local/REMix/remix_nz/input/demand/electricity-authority/output_figures"
os.makedirs(output_dir, exist_ok=True)

# === Load data ===
df = pd.read_csv(file_path, parse_dates=["Period start"])
df["Hour"] = df["Period start"].dt.hour
df["Weekday"] = df["Period start"].dt.dayofweek
df["IsWeekend"] = df["Weekday"] >= 5
df["Month"] = df["Period start"].dt.month
df["Season"] = df["Month"].map({6: "Winter", 7: "Winter", 8: "Winter", 12: "Summer", 1: "Summer", 2: "Summer"}).fillna("Other")

# === Profiles by region/hour ===
def avg_profile(data, condition):
    return data.loc[condition].groupby(["Region ID", "Hour"])["Demand (GWh)"].mean().unstack()

profiles = {
    "wd_all": avg_profile(df[df["IsWeekend"] == False], df["Season"].isin(["Winter", "Summer", "Other"])),
    "we_all": avg_profile(df[df["IsWeekend"] == True], df["Season"].isin(["Winter", "Summer", "Other"])),
    "wd_winter": avg_profile(df[df["IsWeekend"] == False], df["Season"] == "Winter"),
    "we_winter": avg_profile(df[df["IsWeekend"] == True], df["Season"] == "Winter"),
    "wd_summer": avg_profile(df[df["IsWeekend"] == False], df["Season"] == "Summer"),
    "we_summer": avg_profile(df[df["IsWeekend"] == True], df["Season"] == "Summer")
}

# === Filter valid nodes ===
valid_regions = (profiles["wd_all"].sum(axis=1) > 0.01) & (profiles["we_all"].sum(axis=1) > 0.01)
excluded_nodes = profiles["wd_all"].index.difference(valid_regions[valid_regions].index)

for key in profiles:
    profiles[key] = profiles[key][valid_regions]

# === Normalize & Cluster ===
weekday_norm = profiles["wd_all"].div(profiles["wd_all"].sum(axis=1), axis=0)
weekend_norm = profiles["we_all"].div(profiles["we_all"].sum(axis=1), axis=0)
combined_profile = pd.concat([weekday_norm.add_prefix("wd_"), weekend_norm.add_prefix("we_")], axis=1)
combined_profile = combined_profile.dropna()

# KMeans
kmeans = KMeans(n_clusters=3, random_state=42)
combined_profile["Cluster"] = kmeans.fit_predict(combined_profile)

# Label Clusters
cluster_labels = {}
for c in combined_profile["Cluster"].unique():
    avg_wd = combined_profile[combined_profile["Cluster"] == c].iloc[:, :24].mean()
    peak_hour = avg_wd.idxmax()
    peak_val = avg_wd.max()
    hour_num = int(peak_hour.split("_")[1])

    if peak_val > 0.05 and 9 <= hour_num <= 17:
        cluster_labels[c] = "Commercial"
    elif avg_wd.std() > 0.005:
        cluster_labels[c] = "Residential"
    else:
        cluster_labels[c] = "Industrial"

combined_profile["Inferred Sector"] = combined_profile["Cluster"].map(cluster_labels)

# Merge sectors into main df
df["Inferred Sector"] = df["Region ID"].map(combined_profile["Inferred Sector"])
df["Excluded"] = df["Inferred Sector"].isna()

# === Demand summaries ===
excluded_df = df[df["Excluded"]]
included_df = df[~df["Excluded"]]

total_gwh = df["Demand (GWh)"].sum()
peak_total = df.groupby("Region ID")["Demand (GWh)"].max().sum()
excluded_gwh = excluded_df["Demand (GWh)"].sum()
excluded_peak = excluded_df.groupby("Region ID")["Demand (GWh)"].max().sum()

print("\n🕵️ EXCLUDED REGION SUMMARY")
print(f"📌 Regions excluded: {excluded_df['Region ID'].nunique()}")
print(f"⚡ Energy excluded: {excluded_gwh:.2f} GWh ({excluded_gwh / total_gwh:.2%})")
print(f"🔺 Peak demand excluded: {excluded_peak:.2f} GWh ({excluded_peak / peak_total:.2%})")

print("\n❌ Nodes not included due to insufficient data:")
print(", ".join(sorted(excluded_nodes)))

# Sector summaries
sector_energy = included_df.groupby("Inferred Sector")["Demand (GWh)"].sum()
sector_peak = included_df.groupby("Inferred Sector").apply(lambda x: x.groupby("Region ID")["Demand (GWh)"].max().sum(),include_groups=False)

print("\n🔍 SECTOR DEMAND SUMMARY")
for sector in sector_energy.index:
    print(f"🏷️ {sector}")
    print(f"   • Total demand: {sector_energy[sector]:,.2f} GWh ({sector_energy[sector] / total_gwh:.2%})")
    print(f"   • Peak demand: {sector_peak[sector]:,.2f} GWh ({sector_peak[sector] / peak_total:.2%})")

# === Aggregated Hourly Plot by Sector ===
hourly_sector = included_df.groupby(["Inferred Sector", "Hour"])["Demand (GWh)"].mean().unstack()

plt.figure(figsize=(10, 5))
for sector in hourly_sector.index:
    plt.plot(hourly_sector.columns, hourly_sector.loc[sector], marker='o', label=sector)
plt.title("Average Hourly Demand by Sector")
plt.xlabel("Hour of Day")
plt.ylabel("Avg Demand (GWh)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(f"{output_dir}/hourly_profile_by_sector.png")
plt.show()

# === Seasonal Curves (multi-panel) ===
sectors = combined_profile["Inferred Sector"].unique()
fig, axs = plt.subplots(1, len(sectors), figsize=(6 * len(sectors), 4), sharey=True)

for i, sector in enumerate(sectors):
    cluster_ids = combined_profile[combined_profile["Inferred Sector"] == sector].index
    wd_win = profiles["wd_winter"].loc[profiles["wd_winter"].index.intersection(cluster_ids)].mean()
    wd_sum = profiles["wd_summer"].loc[profiles["wd_summer"].index.intersection(cluster_ids)].mean()
    we_win = profiles["we_winter"].loc[profiles["we_winter"].index.intersection(cluster_ids)].mean()
    we_sum = profiles["we_summer"].loc[profiles["we_summer"].index.intersection(cluster_ids)].mean()

    ax = axs[i]
    ax.plot(wd_win.index, wd_win / wd_win.sum(), label="Weekday Winter", color="blue", marker="o")
    ax.plot(wd_sum.index, wd_sum / wd_sum.sum(), label="Weekday Summer", color="blue", linestyle="--", marker="o")
    ax.plot(we_win.index, we_win / we_win.sum(), label="Weekend Winter", color="orange", marker="x")
    ax.plot(we_sum.index, we_sum / we_sum.sum(), label="Weekend Summer", color="orange", linestyle="--", marker="x")

    ax.set_title(f"{sector} ({len(cluster_ids)} nodes)")
    ax.set_xlabel("Hour")
    if i == 0:
        ax.set_ylabel("Normalized Demand")
    ax.grid(True)
    ax.legend()

fig.suptitle("Clustered Load Shapes by Sector – Seasonal View", fontsize=14)
plt.tight_layout()
fig.savefig(f"{output_dir}/seasonal_profiles_by_sector.png")
plt.show()

# === Heatmap ===
heatmap_df = included_df.groupby(["Inferred Sector", "Hour"])["Demand (GWh)"].sum().reset_index()
pivot_df = heatmap_df.pivot(index="Inferred Sector", columns="Hour", values="Demand (GWh)")
pivot_df_normalized = pivot_df.div(pivot_df.sum(axis=1), axis=0)

plt.figure(figsize=(12, 5))
sns.heatmap(pivot_df_normalized, cmap="YlGnBu", cbar_kws={"label": "Share of Daily Demand"})
plt.title("Hourly Demand Distribution by Sector (Normalized)")
plt.ylabel("Inferred Sector")
plt.xlabel("Hour of Day")
plt.tight_layout()
plt.savefig(f"{output_dir}/heatmap_hourly_demand_by_sector.png")
plt.show()


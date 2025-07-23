# This code classifies actual demand curves by similarity to representative profiles (residential, industrial, commercial) for NI and SI.
# It also includes deeper diagnostics and enhanced visualizations.

import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler, normalize
from sklearn.metrics.pairwise import cosine_similarity
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime

# --- 1. File paths ---
input_path = "C:/Local/REMix/remix_nz/input/demand/electricity-authority/"
profile_path = f"{input_path}existing-profiles/NI_SI_Res_Ind_Com_Profile- from EMI Retail category Datasets Additional information  Supporting information and analysis  2014.csv"
demand_path = f"{input_path}hourly_demand_all_2019_with_regions.csv"

# --- 2. Load and convert EMI profile to hourly ---
print("\U0001F4E5 Loading EMI profile data...")
raw_profiles = pd.read_csv(profile_path, skiprows=3)

cols = [
    "NI Residential", "NI Heavy Ind.", "NI Ind & Comm.",
    "SI Residential", "SI Heavy Ind.", "SI Ind & Comm."
]
profiles = raw_profiles[cols].copy()
profiles = profiles.apply(pd.to_numeric, errors="coerce")

# Convert to hourly
hourly_profiles = profiles.groupby(profiles.index // 2).sum().reset_index(drop=True)
hourly_profiles = hourly_profiles.iloc[:8760]

# Normalize EMI profiles
scaler = MinMaxScaler()
profile_curves = pd.DataFrame(scaler.fit_transform(hourly_profiles), columns=hourly_profiles.columns)

# --- 3. Load and prepare demand data ---
print("\n\U0001F4E5 Loading hourly demand data...")
demand_df = pd.read_csv(demand_path, parse_dates=["Period start", "Period end"])
gwh_col = next((col for col in demand_df.columns if "gwh" in col.lower()), None)
print(f"✅ Using demand column: '{gwh_col}'")

# --- 4. Pivot and normalize demand ---
pivot_df = demand_df.pivot_table(index="Region Code", columns="Period start", values=gwh_col).dropna()
demand_normalized = normalize(pivot_df, axis=1)

node_island_map = demand_df.drop_duplicates("Region Code").set_index("Region Code")["Island"].to_dict()

# --- 5. Match demand shapes to profiles ---
match_results = []

for node_id, demand_vec in zip(pivot_df.index, demand_normalized):
    island = node_island_map.get(node_id, "Unknown")
    
    if island == "NI":
        reference_profiles = profile_curves[["NI Residential", "NI Heavy Ind.", "NI Ind & Comm."]]
        labels = ["NI Residential", "NI Heavy Ind.", "NI Ind & Comm."]
    elif island == "SI":
        reference_profiles = profile_curves[["SI Residential", "SI Heavy Ind.", "SI Ind & Comm."]]
        labels = ["SI Residential", "SI Heavy Ind.", "SI Ind & Comm."]
    else:
        reference_profiles = pd.DataFrame()
        labels = []

    # Primary match
    sim_scores = cosine_similarity([demand_vec], reference_profiles.T)[0]
    best_match = labels[np.argmax(sim_scores)]
    best_score = np.max(sim_scores)

    # If low confidence, try all profiles
    if best_score < 0.9:
        all_profiles = profile_curves[cols]
        all_labels = cols
        all_scores = cosine_similarity([demand_vec], all_profiles.T)[0]
        alt_match = all_labels[np.argmax(all_scores)]
        alt_score = np.max(all_scores)

        best_match = alt_match
        best_score = alt_score

    match_results.append({
        "Region Code": node_id,
        "Island": island,
        "Best Match Profile": best_match,
        "Similarity Score": best_score
    })

# --- 6. Results ---
match_df = pd.DataFrame(match_results)
match_df.to_csv("node_profile_matches.csv", index=False)
print("\n📁 Match results saved to: node_profile_matches.csv")

# --- 7. Analysis ---
print("\n📊 Node classification counts:")
print(match_df.groupby(["Island", "Best Match Profile"]).size().unstack(fill_value=0))

print("\n📈 Average similarity scores by profile:")
print(match_df.groupby("Best Match Profile")["Similarity Score"].mean())

# Boxplot
plt.figure(figsize=(10, 5))
sns.boxplot(data=match_df, x="Best Match Profile", y="Similarity Score")
plt.xticks(rotation=45)
plt.title("Similarity Score by Profile")
plt.tight_layout()
plt.show()

# Histogram
plt.figure(figsize=(10, 5))
plt.hist(match_df["Similarity Score"], bins=30, color="skyblue", edgecolor="k")
plt.axvline(0.9, color='r', linestyle='--', label='0.9 Threshold')
plt.title("Distribution of Similarity Scores")
plt.legend()
plt.tight_layout()
plt.show()

# Scatter: Similarity vs Avg Demand
avg_demand = demand_df.groupby("Region Code")[gwh_col].mean()
match_df = match_df.merge(avg_demand, on="Region Code")
plt.figure(figsize=(8, 5))
sns.scatterplot(data=match_df, x=gwh_col, y="Similarity Score", hue="Best Match Profile")
plt.title("Similarity vs Average Demand")
plt.tight_layout()
plt.show()

# Heatmap of all similarity scores
sim_matrix = cosine_similarity(demand_normalized, profile_curves[cols].T)
sns.heatmap(sim_matrix, cmap="viridis", xticklabels=cols)
plt.title("Similarity Scores: Nodes vs Profiles")
plt.xlabel("Profile")
plt.ylabel("Node Index")
plt.tight_layout()
plt.show()

# --- 8. Visualize a Day (low opacity lines) ---
summer_day = "2019-01-15"
day_start = pd.to_datetime(f"{summer_day} 00:00")
day_end = pd.to_datetime(f"{summer_day} 23:00")

for profile in ["NI Residential", "NI Heavy Ind.", "NI Ind & Comm."]:
    node_ids = match_df[match_df["Best Match Profile"] == profile]["Region Code"]
    day_data = demand_df[(demand_df["Region Code"].isin(node_ids)) &
                         (demand_df["Period start"] >= day_start) &
                         (demand_df["Period start"] <= day_end)]
    day_data = day_data.drop_duplicates(subset=["Region Code", "Period start"])
    
    pivot_day = day_data.pivot(index="Region Code", columns="Period start", values=gwh_col).dropna()
    norm_day = ((pivot_day.T - pivot_day.T.min()) / (pivot_day.T.max() - pivot_day.T.min())).T

    plt.figure(figsize=(12, 5))
    for row in norm_day.values:
        plt.plot(range(24), row, alpha=0.1, color="blue")

    profile_hourly = profile_curves[profile].values[14*24:15*24]
    plt.plot(range(24), profile_hourly, label=f"EMI {profile}", color="black", linewidth=2)
    plt.title(f"{profile} – Normalized Demand on {summer_day}")
    plt.xlabel("Hour")
    plt.ylabel("Normalized Demand")
    plt.legend()
    plt.tight_layout()
    plt.show()

for profile in ["SI Residential", "SI Heavy Ind.", "SI Ind & Comm."]:
    node_ids = match_df[match_df["Best Match Profile"] == profile]["Region Code"]
    day_data = demand_df[(demand_df["Region Code"].isin(node_ids)) &
                         (demand_df["Period start"] >= day_start) &
                         (demand_df["Period start"] <= day_end)]
    day_data = day_data.drop_duplicates(subset=["Region Code", "Period start"])
    
    pivot_day = day_data.pivot(index="Region Code", columns="Period start", values=gwh_col).dropna()
    norm_day = ((pivot_day.T - pivot_day.T.min()) / (pivot_day.T.max() - pivot_day.T.min())).T

    plt.figure(figsize=(12, 5))
    for row in norm_day.values:
        plt.plot(range(24), row, alpha=0.1, color="blue")

    profile_hourly = profile_curves[profile].values[14*24:15*24]
    plt.plot(range(24), profile_hourly, label=f"EMI {profile}", color="black", linewidth=2)
    plt.title(f"{profile} – Normalized Demand on {summer_day}")
    plt.xlabel("Hour")
    plt.ylabel("Normalized Demand")
    plt.legend()
    plt.tight_layout()
    plt.show()

# --- 9. Export normalized profiles ---
emi_export = profile_curves.copy()
emi_export["Hour"] = range(1, 8761)
emi_export.to_csv("normalized_emi_profiles.csv", index=False)
print("\n📁 Normalized EMI profiles exported to: normalized_emi_profiles.csv")

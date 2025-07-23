import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# === Load full demand and classification data ===
demand_df = pd.read_csv("hourly_demand_all_2019_with_regions.csv", parse_dates=["Period start"])
match_df = pd.read_csv("node_profile_matches.csv")

# === Classify into sectors ===
def classify_sector(profile):
    if "Residential" in profile:
        return "Residential"
    elif "Heavy Ind." in profile:
        return "Industrial"
    elif "Ind & Comm." in profile:
        return "Commercial"
    else:
        return "Unclassified"

match_df["Sector"] = match_df["Best Match Profile"].apply(classify_sector)

# === Print node sector classification ===
print("\n📋 Nodes per sector:")
for sector in ["Residential", "Industrial", "Commercial", "Unclassified"]:
    codes = match_df[match_df["Sector"] == sector]["Region Code"].tolist()
    print(f"\n🔹 {sector} ({len(codes)} nodes):\n{codes}")

# === Add date structure ===
demand_df = demand_df.rename(columns={"Period start": "Fecha/Hora"})
gwh_col = next(c for c in demand_df.columns if "gwh" in c.lower())
demand_df["yr"] = demand_df["Fecha/Hora"].dt.year
demand_df["mes"] = demand_df["Fecha/Hora"].dt.month
demand_df["dia"] = demand_df["Fecha/Hora"].dt.day
demand_df["hora"] = demand_df["Fecha/Hora"].dt.hour

# === Create 8760-row demand matrix for each sector ===
def make_sector_df(sector_name):
    codes = match_df[match_df["Sector"] == sector_name]["Region Code"].tolist()
    df = demand_df[demand_df["Region Code"].isin(codes)].copy()

    # Pivot using pivot_table to safely handle duplicates
    pivot = df.pivot_table(index="Fecha/Hora", columns="Region Code", values=gwh_col, aggfunc="sum")

    pivot = pivot.reset_index()
    pivot["yr"] = pivot["Fecha/Hora"].dt.year
    pivot["mes"] = pivot["Fecha/Hora"].dt.month
    pivot["dia"] = pivot["Fecha/Hora"].dt.day
    pivot["hora"] = pivot["Fecha/Hora"].dt.hour
    return pivot


df_res = make_sector_df("Residential")
df_com = make_sector_df("Commercial")
df_ind = make_sector_df("Industrial")
df_unc = make_sector_df("Unclassified")

# === Preview ===
print("\n🧾 Residential Demand DataFrame:")
print(df_res.head())

print("\n🧾 Industrial Demand DataFrame:")
print(df_ind.head())

# === 📊 Example plot using your format ===
# Choose 2 region codes to compare
node_1 = df_res.columns[1]  # 1st node in Residential
node_2 = df_com.columns[1]  # 1st node in Commercial

sns.set(style="whitegrid")

fig, axes = plt.subplots(1, 2, figsize=(16, 6), dpi=300, facecolor="w", sharey=True)
fig.subplots_adjust(wspace=0.03)

axes[0].grid(True, color='gainsboro', linestyle='dashed')
axes[1].grid(True, color='gainsboro', linestyle='dashed')

sns.lineplot(x=df_res["hora"], y=df_res[node_1], hue=df_res["mes"], palette='Blues', ax=axes[0])
sns.lineplot(x=df_com["hora"], y=df_com[node_2], hue=df_com["mes"], palette='summer_r', ax=axes[1])

sns.lineplot(x=df_res["hora"], y=df_res[node_1], color='k', ax=axes[0], label="Promedio", ci=None, legend=False)
sns.lineplot(x=df_com["hora"], y=df_com[node_2], color='k', ax=axes[1], label="Promedio", ci=None, legend=False)

axes[0].set_title(node_1)
axes[1].set_title(node_2)

axes[0].set_ylabel('Demanda (GWh)')
axes[1].set_ylabel('Demanda (GWh)')
axes[0].set_xlabel('\n Hora')
axes[1].set_xlabel('\n Hora')

xticks = range(0, 24, 2)
xticklabels = ['\n 0:00','\n 2:00','\n 4:00','\n 6:00','\n 8:00','\n 10:00',
               '\n 12:00','\n 14:00','\n 16:00','\n 18:00','\n 20:00','\n 22:00']
axes[0].set_xticks(xticks)
axes[0].set_xticklabels(xticklabels, rotation=45)
axes[1].set_xticks(xticks)
axes[1].set_xticklabels(xticklabels, rotation=45)

axes[0].legend(title='Mes')
axes[1].legend(title='Mes')

plt.suptitle("Demanda promedio por hora del día")
plt.tight_layout()
plt.show()

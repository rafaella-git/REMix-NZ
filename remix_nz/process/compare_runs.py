# compare_runs.py
import gdxpds
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import os

# Config
group_name = "hadi"
base_case = "pypsa"
scenario_case = "pypsa-low"
label1, label2 = "Base", "Dispatch"
target_year = "2030"
pairs = [("Inflow", "Water_in"), ("Wholesale", "Elec"), ("Slack", "Elec")]

sns.set_theme(style="whitegrid", context="talk")
os.chdir(Path(__file__).parent.resolve())

# Paths
gdx1_path = f"../project/{group_name}/{base_case}/result/{base_case}_opt.gdx"
gdx2_path = f"../project/{group_name}/{scenario_case}/result/{scenario_case}_dispatch.gdx"

# ------------------------------------------------------------------------------

def load_gdx_data(path):
    p = Path(path).resolve()
    if not p.exists():
        raise FileNotFoundError(f"GDX not found: {p}")
    data = gdxpds.to_dataframes(str(p))
    print(f"Loaded {len(data)} tables from {path}")
    return data


def compare_converter_caps(df1, df2, year, label1, label2):
    d1, d2 = df1.reset_index().copy(), df2.reset_index().copy()
    d1.columns = d1.columns.str.lower(); d2.columns = d2.columns.str.lower()
    d1["accyears"], d2["accyears"] = d1["accyears"].astype(str), d2["accyears"].astype(str)
    d1 = d1[(d1["accyears"] == year) & (d1["captype"].str.lower() == "total") &
            (d1["commodity"].str.lower() == "elec") &
            (~d1["accnodesmodel"].str.lower().eq("global"))]
    d2 = d2[(d2["accyears"] == year) & (d2["captype"].str.lower() == "total") &
            (d2["commodity"].str.lower() == "elec") &
            (~d2["accnodesmodel"].str.lower().eq("global"))]
    d1 = d1.groupby(["techs"], as_index=False)["value"].sum()
    d2 = d2.groupby(["techs"], as_index=False)["value"].sum()

    comp = pd.merge(d1, d2, on="techs", how="outer",
                    suffixes=(f"_{label1}", f"_{label2}"))
    comp["Difference"] = comp[f"value_{label2}"].fillna(0) - comp[f"value_{label1}"].fillna(0)
    comp = comp.round(2)

    print("\n--- Installed Capacity Comparison (2030) ---")
    print(comp)
    total1, total2 = comp[f"value_{label1}"].sum(), comp[f"value_{label2}"].sum()
    print(f"\nTotal installed capacity ({label1}): {total1:.2f} GW")
    print(f"Total installed capacity ({label2}): {total2:.2f} GW")
    if total1 != 0:
        print(f"Change: {(total2 - total1) / total1 * 100:+.2f}%")
    if comp["Difference"].abs().sum() < 0.01:
        print("All capacities identical between runs.")

    return comp


def compare_balance(df1, df2, year, label1, label2, pairs):
    d1, d2 = df1.reset_index().copy(), df2.reset_index().copy()
    d1.columns = d1.columns.str.lower(); d2.columns = d2.columns.str.lower()
    d1["accyears"], d2["accyears"] = d1["accyears"].astype(str), d2["accyears"].astype(str)
    d1 = d1[(d1["accyears"] == year) & (~d1["accnodesmodel"].str.lower().eq("global")) &
            (d1["balancetype"].str.lower() == "net")]
    d2 = d2[(d2["accyears"] == year) & (~d2["accnodesmodel"].str.lower().eq("global")) &
            (d2["balancetype"].str.lower() == "net")]

    mask1 = pd.Series(False, index=d1.index); mask2 = pd.Series(False, index=d2.index)
    for tech, comm in pairs:
        mask1 |= (d1["techs"] == tech) & (d1["commodity"] == comm)
        mask2 |= (d2["techs"] == tech) & (d2["commodity"] == comm)
    d1, d2 = d1[mask1], d2[mask2]
    d1 = d1.groupby(["techs"], as_index=False)["value"].sum()
    d2 = d2.groupby(["techs"], as_index=False)["value"].sum()

    comp = pd.merge(d1, d2, on="techs", how="outer",
                    suffixes=(f"_{label1}", f"_{label2}"))
    comp["Difference"] = comp[f"value_{label2}"].fillna(0) - comp[f"value_{label1}"].fillna(0)
    comp = comp.round(2)

    print("\n--- Commodity Balance Comparison (2030) ---")
    print(comp)
    total1, total2 = comp[f"value_{label1}"].sum(), comp[f"value_{label2}"].sum()
    print(f"\nTotal net flow ({label1}): {total1:.2f} TWh")
    print(f"Total net flow ({label2}): {total2:.2f} TWh")
    print(f"Change: {total2 - total1:+.2f} TWh")
    if comp["Difference"].abs().sum() < 0.01:
        print("All flows identical between runs.")

    return comp


def derive_dispatch(df1, df2, label1, label2):
    def extract_dispatch(df):
        d = df.reset_index().copy()
        d.columns = d.columns.str.lower()
        d = d[
            (d["accnodesmodel"].str.lower() == "global") &
            (d["commodity"].str.lower().isin(["elec", "h2"])) &
            (d["balancetype"].str.lower() == "net") &
            (d["accyears"].astype(str) == target_year)
        ]
        return d.groupby(["techs", "commodity"], as_index=False)["value"].sum()

    d1, d2 = extract_dispatch(df1), extract_dispatch(df2)
    comp = pd.merge(d1, d2, on=["techs", "commodity"], how="outer",
                    suffixes=(f"_{label1}", f"_{label2}"))
    comp["Difference"] = comp[f"value_{label2}"].fillna(0) - comp[f"value_{label1}"].fillna(0)
    comp = comp.round(2)

    print("\n--- Global Dispatch Comparison (2030, Elec + H₂) ---")
    print(comp)
    if comp["Difference"].abs().sum() < 0.01:
        print("Dispatch identical between runs.")

    return comp


def plot_combined(capacity_df, balance_df, dispatch_df, label1, label2):
    sns.set_theme(style="whitegrid", context="talk")
    palette = sns.color_palette("Set2", n_colors=2)
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    plt.subplots_adjust(wspace=0.3)
    l1, l2 = label1.lower(), label2.lower()

    # 1. Capacity
    cap = capacity_df.copy(); cap.columns = cap.columns.str.lower()
    cap_m = cap.melt(id_vars=["techs"], value_vars=[f"value_{l1}", f"value_{l2}"],
                     var_name="Scenario", value_name="Capacity (GW)")
    sns.barplot(data=cap_m, x="techs", y="Capacity (GW)", hue="Scenario",
                ax=axes[0], palette=palette, errorbar=None)
    axes[0].set_title("Installed Capacity (2030)")
    axes[0].tick_params(axis="x", rotation=45)

    # 2. Balance
    bal = balance_df.copy(); bal.columns = bal.columns.str.lower()
    bal_m = bal.melt(id_vars=["techs"], value_vars=[f"value_{l1}", f"value_{l2}"],
                     var_name="Scenario", value_name="Net Flow (TWh)")
    sns.barplot(data=bal_m, x="techs", y="Net Flow (TWh)", hue="Scenario",
                ax=axes[1], palette=palette, errorbar=None)
    axes[1].set_title("Net Energy Flow (2030)")
    axes[1].tick_params(axis="x", rotation=45)

    # 3. Dispatch
    disp = dispatch_df.copy(); disp.columns = disp.columns.str.lower()
    disp_m = disp.melt(id_vars=["techs", "commodity"],
                       value_vars=[f"value_{l1}", f"value_{l2}"],
                       var_name="Scenario", value_name="Output (TWh)")
    sns.barplot(data=disp_m, x="techs", y="Output (TWh)", hue="Scenario",
                ax=axes[2], palette=palette, errorbar=None)
    axes[2].set_title("Global Dispatch (Elec + H₂, 2030)")
    axes[2].tick_params(axis="x", rotation=45)

    for ax in axes:
        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)
        ax.grid(True, axis="y", linestyle="--", alpha=0.6)
        ax.legend(title="Scenario")

    plt.tight_layout()
    plt.show()

# ------------------------------------------------------------------------------

data1 = load_gdx_data(gdx1_path)
data2 = load_gdx_data(gdx2_path)
conv_caps1, conv_caps2 = data1.get("converter_caps"), data2.get("converter_caps")
bal1, bal2 = data1.get("commodity_balance_annual"), data2.get("commodity_balance_annual")

if any(x is None for x in [conv_caps1, conv_caps2, bal1, bal2]):
    raise KeyError("Missing required tables in one of the GDX files.")

cap_df = compare_converter_caps(conv_caps1, conv_caps2, target_year, label1, label2)
bal_df = compare_balance(bal1, bal2, target_year, label1, label2, pairs)
disp_df = derive_dispatch(bal1, bal2, label1, label2)
plot_combined(cap_df, bal_df, disp_df, label1, label2)

print("\n✅ Comparison complete.")

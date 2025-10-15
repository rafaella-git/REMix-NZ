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
years = [2020, 2030]
pairs = [("Inflow", "Water_in"), ("Wholesale", "Elec"), ("Slack", "Elec")]

sns.set_theme(style="whitegrid")
os.chdir(Path(__file__).parent.resolve())

# Paths (match outputs from run_remix.py)
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

def compare_converter_caps(df1, df2, years, label1, label2):
    d1, d2 = df1.reset_index().copy(), df2.reset_index().copy()
    d1.columns = d1.columns.str.lower(); d2.columns = d2.columns.str.lower()
    d1["accyears"], d2["accyears"] = d1["accyears"].astype(str), d2["accyears"].astype(str)
    years = [str(y) for y in years]
    d1 = d1[~d1["accnodesmodel"].str.lower().eq("global")]
    d2 = d2[~d2["accnodesmodel"].str.lower().eq("global")]
    d1 = d1[(d1["captype"].str.lower() == "total") & (d1["commodity"].str.lower() == "elec")]
    d2 = d2[(d2["captype"].str.lower() == "total") & (d2["commodity"].str.lower() == "elec")]
    d1, d2 = d1[d1["accyears"].isin(years)], d2[d2["accyears"].isin(years)]
    d1 = d1.groupby(["accyears", "techs"], as_index=False)["value"].sum()
    d2 = d2.groupby(["accyears", "techs"], as_index=False)["value"].sum()
    comp = pd.merge(d1, d2, on=["accyears", "techs"], how="outer",
                    suffixes=(f"_{label1}", f"_{label2}"))
    comp["Difference"] = comp[f"value_{label2}"].fillna(0) - comp[f"value_{label1}"].fillna(0)
    comp = comp.round(2)
    print("\nInstalled capacity comparison (Elec):")
    print(comp)
    return comp

def compare_balance(df1, df2, years, label1, label2, pairs):
    d1, d2 = df1.reset_index().copy(), df2.reset_index().copy()
    d1.columns = d1.columns.str.lower(); d2.columns = d2.columns.str.lower()
    d1["accyears"], d2["accyears"] = d1["accyears"].astype(str), d2["accyears"].astype(str)
    years = [str(y) for y in years]
    d1 = d1[~d1["accnodesmodel"].str.lower().eq("global")]
    d2 = d2[~d2["accnodesmodel"].str.lower().eq("global")]
    mask1 = pd.Series(False, index=d1.index); mask2 = pd.Series(False, index=d2.index)
    for tech, comm in pairs:
        mask1 |= (d1["techs"] == tech) & (d1["commodity"] == comm)
        mask2 |= (d2["techs"] == tech) & (d2["commodity"] == comm)
    d1, d2 = d1[mask1], d2[mask2]
    d1, d2 = d1[d1["accyears"].isin(years)], d2[d2["accyears"].isin(years)]
    d1, d2 = d1[d1["balancetype"].str.lower() == "net"], d2[d2["balancetype"].str.lower() == "net"]
    d1 = d1.groupby(["accyears", "techs", "commodity"], as_index=False)["value"].sum()
    d2 = d2.groupby(["accyears", "techs", "commodity"], as_index=False)["value"].sum()
    comp = pd.merge(d1, d2, on=["accyears", "techs", "commodity"], how="outer",
                    suffixes=(f"_{label1}", f"_{label2}"))
    comp["Difference"] = comp[f"value_{label2}"].fillna(0) - comp[f"value_{label1}"].fillna(0)
    comp = comp.round(2)
    print("\nCommodity balance comparison:")
    print(comp)
    return comp

def plot_combined(capacity_df, balance_df, label1, label2):
    sns.set_theme(style="whitegrid", context="talk")
    palette = sns.color_palette("Set2", n_colors=2)
    fig, axes = plt.subplots(2, 1, figsize=(12, 10))
    plt.subplots_adjust(hspace=0.35)

    # Capacity
    cap = capacity_df.copy(); cap.columns = cap.columns.str.lower()
    l1, l2 = label1.lower(), label2.lower()
    cap_m = cap.melt(id_vars=["accyears", "techs"],
                     value_vars=[f"value_{l1}", f"value_{l2}"],
                     var_name="Scenario", value_name="Capacity (GW)")
    sns.barplot(data=cap_m, x="techs", y="Capacity (GW)", hue="Scenario",
                ax=axes[0], palette=palette, errorbar=None)
    axes[0].set_title("Installed Electricity Capacity by Technology")
    axes[0].set_ylabel("Capacity (GW)")
    axes[0].tick_params(axis="x", rotation=45)
    axes[0].legend(title="Scenario")

    # Balance
    bal = balance_df.copy(); bal.columns = bal.columns.str.lower()
    bal_m = bal.melt(id_vars=["accyears", "techs"],
                     value_vars=[f"value_{l1}", f"value_{l2}"],
                     var_name="Scenario", value_name="Net Flow (TWh)")
    sns.barplot(data=bal_m, x="techs", y="Net Flow (TWh)", hue="Scenario",
                ax=axes[1], palette=palette, errorbar=None)
    axes[1].set_title("Net Energy Flow by Technology")
    axes[1].set_ylabel("Net Flow (TWh)")
    axes[1].set_xlabel("Technology")
    axes[1].tick_params(axis="x", rotation=45)
    axes[1].legend(title="Scenario")

    for ax in axes:
        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)
        ax.grid(True, axis="y", linestyle="--", alpha=0.6)

    plt.tight_layout()
    plt.show()

# ------------------------------------------------------------------------------

data1 = load_gdx_data(gdx1_path)
data2 = load_gdx_data(gdx2_path)
conv_caps1, conv_caps2 = data1.get("converter_caps"), data2.get("converter_caps")
bal1, bal2 = data1.get("commodity_balance_annual"), data2.get("commodity_balance_annual")

if any(x is None for x in [conv_caps1, conv_caps2, bal1, bal2]):
    raise KeyError("Missing required tables in one of the GDX files.")

cap_df = compare_converter_caps(conv_caps1, conv_caps2, years, label1, label2)
bal_df = compare_balance(bal1, bal2, years, label1, label2, pairs)
plot_combined(cap_df, bal_df, label1, label2)

print("\nComparison complete.")


# possibel improvement, add a number with % of bariation on top of the bars, change name of the x labels to something prettier
# have smaller font and figures be side to side
# do a geographical plot of the differences in capacity and energy output, also on slack
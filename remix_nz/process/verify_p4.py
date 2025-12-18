# # analyse_2050_all_scenarios.py

# import gdxpds  # must come before pandas on Windows

# import os
# from pathlib import Path
# import pandas as pd
# import matplotlib.pyplot as plt

# # ------------------------------------------------------------------------
# # CONFIG
# cases = [
#     ("GP-NT-ELEC-BIO-H2", "nz_case_GP_2050"),
#     ("GP-NT-ELEC-BIO-H2", "nz_case_NT_2050"),
#     ("GP-NT-ELEC-BIO-H2", "nz_case_ELEC+_2050"),
#     ("GP-NT-ELEC-BIO-H2", "nz_case_BIO+_2050"),
#     ("GP-NT-ELEC-BIO-H2", "nz_case_H2+_2050"),
# ]

# YEAR = "2050"

# base_dir = Path(__file__).parent.resolve() / ".." / "project"

# # ------------------------------------------------------------------------
# # HELPERS
# def scen_label(case_name: str) -> str:
#     return case_name.replace("nz_case_", "").replace("_2050", "")

# def load_gdx_tables(gdx_path: Path) -> dict:
#     data = gdxpds.to_dataframes(str(gdx_path))
#     return data

# def tech_group(tech: str) -> str:
#     # Map your builder’s tech names into groups for capacity comparison
#     t = tech
#     if t == "Hydro":
#         return "Hydro"
#     if t.lower().startswith("wind"):
#         return "Wind"
#     if t.lower().startswith("pv"):
#         return "Solar"
#     if t in ["Geoth", "geoth"]:
#         return "Geothermal"
#     if t in ["BIO", "COAL", "DIE"]:
#         return "Thermal_fossil_bio"
#     if t == "Electrolyser":
#         return "Electrolyser"
#     if t == "DAC":
#         return "DAC"
#     if t == "Methanizer":
#         return "Methanizer"
#     if t == "MethanolSyn":
#         return "MethanolSyn"
#     if t == "FTropschSyn":
#         return "FTropsch"
#     if t == "H2CCGT":
#         return "H2CCGT"
#     if t == "H2FC":
#         return "H2FC"
#     return "Other"

# def bar_table(df: pd.DataFrame, title: str, ylabel: str):
#     ax = df.plot(kind="bar", figsize=(10, 6))
#     ax.set_title(title)
#     ax.set_ylabel(ylabel)
#     ax.set_xlabel("")
#     ax.legend(title="Scenario")
#     ax.grid(axis="y", linestyle="--", alpha=0.4)
#     plt.tight_layout()
#     plt.show()

# # ------------------------------------------------------------------------
# # METRIC EXTRACTION FOR ONE SCENARIO (2050)
# # ------------------------------------------------------------------------
# def extract_metrics_2050(data: dict) -> dict:
#     """
#     data: dict from gdxpds.to_dataframes

#     Returns:
#       - capacity_by_group_GW: Series(index=tech_group, value=GW)
#       - final_energy_TWh:     Series(index=commodity, value=TWh demand)
#       - total_cost_MEUR:      float
#       - cost_by_indicator:    Series(index=indicator, value)
#       - storage_throughput_GWh: float (batteries + H2 storage)
#     """
#     m = {}

#     # 1) converter_caps → capacities in 2050
#     cap = data.get("converter_caps")
#     if cap is None or cap.empty:
#         raise KeyError("converter_caps not found or empty in GDX")

#     cap_2050 = cap[cap["accYears"] == YEAR].copy()
#     cap_2050["group"] = cap_2050["techs"].map(tech_group)
#     # Filter to total Elec capacity (GW)
#     mask_total = (cap_2050["capType"] == "total") & (cap_2050["commodity"] == "Elec")
#     cap_2050 = cap_2050[mask_total]
#     cap_by_group = cap_2050.groupby("group")["Value"].sum()
#     m["capacity_by_group_GW"] = cap_by_group

#     # 2) commodity_balance_annual → final energy demand by commodity (TWh)
#     cba = data.get("commodity_balance_annual")
#     if cba is None or cba.empty:
#         raise KeyError("commodity_balance_annual not found or empty in GDX")

#     cba_2050 = cba[(cba["accYears"] == YEAR) & (cba["balanceType"] == "net")].copy()
#     # Negative = demand; convert to positive demand
#     fe_by_comm = -cba_2050.groupby("commodity")["Value"].sum() / 1e3  # TWh
#     m["final_energy_TWh"] = fe_by_comm

#     # 3) indicator_accounting → cost metrics
#     ind = data.get("indicator_accounting")
#     if ind is None or ind.empty:
#         raise KeyError("indicator_accounting not found or empty in GDX")

#     # Use accYears == "horizon" for total cost (typical REMix convention), fallback to YEAR if missing
#     ia_hor = ind[ind["accYears"] == "horizon"]
#     if ia_hor.empty:
#         ia_use = ind[ind["accYears"] == YEAR]
#     else:
#         ia_use = ia_hor

#     total_cost = ia_use[ia_use["indicator"] == "SystemCost"]["Value"].sum()
#     m["total_cost_MEUR"] = total_cost

#     cost_by_indicator = ia_use.groupby("indicator")["Value"].sum()
#     m["cost_by_indicator"] = cost_by_indicator

#     # 4) storage_flows_annual → storage throughput (only batteries + H2 storage)
#     sf = data.get("storage_flows_annual")
#     if sf is None or sf.empty:
#         raise KeyError("storage_flows_annual not found or empty in GDX")

#     sf_2050 = sf[(sf["accYears"] == YEAR) & (sf["balanceType"].isin(["positive", "negative"]))].copy()
#     # keep only Battery and H2_storage
#     mask_batt = sf_2050["techs"].str.contains("Battery", case=False, na=False)
#     mask_h2st = sf_2050["techs"].str.contains("H2_storage", case=False, na=False)
#     sf_2050 = sf_2050[mask_batt | mask_h2st]

#     # Annual throughput ≈ sum of positive flows (GWh)
#     storage_throughput = sf_2050[sf_2050["balanceType"] == "positive"]["Value"].sum()
#     m["storage_throughput_GWh"] = storage_throughput

#     return m

# # ------------------------------------------------------------------------
# # MAIN: loop over scenarios, build tables, plot
# # ------------------------------------------------------------------------
# def main():
#     cap_rows = []
#     fe_rows = []
#     cost_rows = []
#     cost_bd_rows = []
#     ops_rows = []

#     for group_name, case_name in cases:
#         scen = scen_label(case_name)
#         gdx_path = base_dir / group_name / case_name / "result" / f"{case_name}.gdx"
#         if not gdx_path.exists():
#             print(f"[WARN] Missing GDX for {scen}: {gdx_path}")
#             continue

#         print(f"\n=== Scenario: {scen} ===")
#         data = load_gdx_tables(gdx_path)

#         metrics = extract_metrics_2050(data)

#         cap_rows.append(metrics["capacity_by_group_GW"].rename(scen))
#         fe_rows.append(metrics["final_energy_TWh"].rename(scen))

#         cost_rows.append(pd.Series({
#             "scenario": scen,
#             "total_cost_MEUR_2050": metrics["total_cost_MEUR"],
#         }))
#         cost_bd_rows.append(metrics["cost_by_indicator"].rename(scen))

#         ops_rows.append(pd.Series({
#             "scenario": scen,
#             "storage_throughput_GWh_2050": metrics["storage_throughput_GWh"],
#         }))

#     # ---------- Table 1: capacities ----------
#     tbl_cap = pd.concat(cap_rows, axis=1).fillna(0.0)
#     tbl_cap.index.name = "Tech_group"
#     print("\nTable 1: Installed capacity by tech group in 2050 (GW)")
#     print(tbl_cap)
#     bar_table(tbl_cap, "Installed capacity in 2050 by tech group", "GW")

#     # ---------- Table 2: final energy demand ----------
#     tbl_fe = pd.concat(fe_rows, axis=1).fillna(0.0)
#     tbl_fe.index.name = "Commodity"
#     print("\nTable 2: Final energy demand in 2050 (TWh)")
#     print(tbl_fe)
#     bar_table(tbl_fe, "Final energy demand in 2050", "TWh")

#     # ---------- Table 3: total cost ----------
#     tbl_cost = pd.DataFrame(cost_rows).set_index("scenario")
#     print("\nTable 3: Total SystemCost in 2050 (M€)")
#     print(tbl_cost)
#     bar_table(tbl_cost[["total_cost_MEUR_2050"]], "Total SystemCost in 2050", "M€")

#     # ---------- Table 4: cost breakdown ----------
#     tbl_cost_bd = pd.concat(cost_bd_rows, axis=1).fillna(0.0)
#     tbl_cost_bd.index.name = "Indicator"
#     print("\nTable 4: Cost breakdown by indicator in 2050 (M€)")
#     print(tbl_cost_bd)
#     bar_table(tbl_cost_bd, "Cost breakdown by indicator in 2050", "M€")

#     # ---------- Table 5: storage throughput ----------
#     tbl_ops = pd.DataFrame(ops_rows).set_index("scenario")
#     print("\nTable 5: Storage throughput (Battery + H2_storage) in 2050 (GWh)")
#     print(tbl_ops)
#     bar_table(tbl_ops[["storage_throughput_GWh_2050"]],
#               "Storage throughput in 2050 (Battery + H2_storage)", "GWh")

# if __name__ == "__main__":
#     main()

#!/usr/bin/env python3
"""
figure_installed_capacity_2020_2050_H2plus.py

Create Figure:
Installed capacity by technology group in 2020 and 2050 (GW)
for scenario: nz_case_H2+_2020-2050

Requirements:
    pip install gdxpds pandas matplotlib
"""

import gdxpds  # must come before pandas on Windows
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap

# ----------------------------------------------------------------------
# CONFIG
YEARS = ["2020", "2050"]

# Adjust these three to match your project layout
BASE_GROUP = "GP-NT-ELEC-BIO-H2"
CASE_NAME = "nz_case_H2+_2020-2050"
GDX_ROOT = Path(__file__).parent.resolve() / ".." / "project"

GDX_PATH = GDX_ROOT / BASE_GROUP / CASE_NAME / "result" / f"{CASE_NAME}.gdx"

FIGURE_OUTFILE = (
    Path(__file__).parent.resolve()
    / f"figure_installed_capacity_{CASE_NAME}_2020_2050.png"
)

# ----------------------------------------------------------------------
# TECHNOLOGY GROUPING
TECH_GROUP_DEFINITION = {
    "Hydropower": ["Hydro"],
    "Wind onshore": ["wind_onshore"],
    "Wind offshore": ["wind_offshore_foundation", "wind_offshore_floating"],
    "Solar PV": ["pv_decentral", "pv_central_fixed"],
    "Geothermal": ["geoth"],
    "Thermal": ["BIO", "COAL", "DIE"],
    "Electrolyser": ["Electrolyser"],
    "Gas turbines": ["GT", "CCGT", "OCGT"],
    "DAC": ["DAC"],
    "Methaniser": ["Methanizer"],
    "Fischer-Tropsch": ["FTropschSyn"],
    "Modified gas turbines (H2)": ["GT_H2", "CCGT_H2"],
    "H2 fuel cell": ["H2_FC"],
    "Battery": ["Battery"],
    "H2 storage": ["H2_storage"],
    "Hydro reservoir": ["Hydro_reservoir"],
}

# Flatten mapping: tech -> group
TECH_TO_GROUP = {}
for group_name, tech_list in TECH_GROUP_DEFINITION.items():
    for t in tech_list:
        TECH_TO_GROUP[t] = group_name


# ----------------------------------------------------------------------
# HELPERS
def load_gdx_tables(gdx_path: Path) -> dict:
    if not gdx_path.exists():
        raise FileNotFoundError(f"GDX file not found: {gdx_path}")
    data = gdxpds.to_dataframes(str(gdx_path))
    return data


def map_tech_to_group(tech: str) -> str:
    """
    Map a raw tech name to one of the predefined technology groups.
    Unmapped techs are returned as 'Other'.
    """
    if tech in TECH_TO_GROUP:
        return TECH_TO_GROUP[tech]
    return "Other"


def extract_capacities_by_group_for_year(data: dict, year: str) -> pd.Series:
    """
    From converter_caps, extract total electric capacity (GW) for a given year
    aggregated by technology group.
    """
    cap = data.get("converter_caps")
    if cap is None or cap.empty:
        raise KeyError("converter_caps not found or empty in GDX")

    cap_y = cap[
        (cap["accYears"] == year)
        & (cap["commodity"] == "Elec")
        & (cap["capType"] == "total")
        & (cap["accNodesModel"] == "global")
    ].copy()

    cap_y["Tech_group"] = cap_y["techs"].map(map_tech_to_group)

    cap_by_group = cap_y.groupby("Tech_group")["Value"].sum()

    ordered_groups = list(TECH_GROUP_DEFINITION.keys())
    cap_by_group = cap_by_group.reindex(ordered_groups, fill_value=0.0)

    return cap_by_group


def make_grouped_bar_figure(cap_df: pd.DataFrame) -> None:
    """
    Produce the grouped bar figure:
    X-axis: technology groups
    Two bars per group: 2020 and 2050.
    """
    tech_groups = cap_df.index.tolist()
    x = np.arange(len(tech_groups))
    width = 0.35

    fig, ax = plt.subplots(figsize=(12, 6))

    cmap = get_cmap("viridis")
    # Sample two distinct colors from viridis
    color_2020 = cmap(0.2)
    color_2050 = cmap(0.8)

    bars_2020 = ax.bar(
        x - width / 2,
        cap_df["2020"].values,
        width,
        label="2020",
        color=color_2020,
    )
    bars_2050 = ax.bar(
        x + width / 2,
        cap_df["2050"].values,
        width,
        label="2050",
        color=color_2050,
    )

    ax.set_xlabel("Technology groups")
    ax.set_ylabel("Installed capacity (GW) in 2050")
    ax.set_xticks(x)
    ax.set_xticklabels(tech_groups, rotation=45, ha="right")

    ax.set_ylim(bottom=0)
    ax.legend(title="Year")
    ax.grid(axis="y", linestyle="--", alpha=0.4)

    fig.tight_layout()
    fig.savefig(FIGURE_OUTFILE, dpi=300)
    print(f"Figure saved to: {FIGURE_OUTFILE}")


def print_stats(cap_df: pd.DataFrame, scenario_label: str) -> None:
    """
    Print some relevant stats to console for both years.
    """
    print("\n=== Installed capacity by tech group (GW) ===")
    print(f"Scenario: {scenario_label}")
    print(cap_df.to_string())

    for year in cap_df.columns:
        total_cap = cap_df[year].sum()
        top3 = cap_df[year].sort_values(ascending=False).head(3)

        print(f"\n--- Year {year} ---")
        print(f"Total installed capacity: {total_cap:,.2f} GW")
        print("Top 3 tech groups by capacity:")
        for name, val in top3.items():
            share = 100.0 * val / total_cap if total_cap > 0 else 0.0
            print(f"  {name:30s}: {val:10.2f} GW  ({share:5.1f}%)")


# ----------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------
def main():
    scenario_label = CASE_NAME.replace("nz_case_", "").replace("_2020-2050", "")
    print(f"Loading GDX for scenario {scenario_label} ...")
    print(f"  {GDX_PATH}")

    data = load_gdx_tables(GDX_PATH)

    cap_by_year = {}
    for year in YEARS:
        cap_by_year[year] = extract_capacities_by_group_for_year(data, year)

    # Combine into a single DataFrame with columns = years, rows = tech groups
    cap_df = pd.DataFrame(cap_by_year)

    print_stats(cap_df, scenario_label)
    make_grouped_bar_figure(cap_df)


if __name__ == "__main__":
    main()

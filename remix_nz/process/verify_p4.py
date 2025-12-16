# Verify that 2050 fuel/e-fuel use in the REMix result GDX respects
# the annual upper bounds implied by the carriers Excel.
#
# It reconstructs the 2050 annual "demand" per (node, sector, commodity)
# from datasummary.xlsx (using the same mappings as the build script),
# then compares it to sourcesink_flow (or sourcesink_flow_annual) in GDX.
# ----------------------------------------------------------------------

import os
from pathlib import Path

import gdxpds        
import numpy as np   
import pandas as pd 

# ----------------- USER SETTINGS --------------------------------------

# GDX result file
GDX_PATH = Path("C:/Local/REMix/remix_nz/project/GP-NT-ELEC-BIO-H2/nz_case_GP_MVP/result/nz_case_GP_MVP.gdx")
# Carriers Excel used by working_script.txt
CARRIERS_EXCEL = Path("C:/Local/REMix/remix_nz/input/demand/GP-NT-ELEC-BIO-H2/data_summary.xlsx")
BASE_SCENARIO = "GP"      # basescenario in working_script.txt 
YEAR_CHECK = 2050

# Commodities to check (as in your add_paid_fuel_consumption_from_excel + e-fuels). 
CHECK_COMMODITIES = [
    "LF_fossil",
    "Gas_fossil",
    "Coal_fossil",
    "LF_bio",
    "Gas_bio",
    "Wood",
    "REfuel",
    "CH4",
    "H2",       # include H2 if you want to check it
]



# ----------------- HELPERS FROM BUILD SCRIPT --------------------------

def map_region_to_remix(region_value: str) -> str:
    """Simplified version of mapRegionToRemix from working_script.txt."""  #
    region_to_remix = {
        "Auckland": "AKL",
        "Bay of Plenty": "BOP",
        "Canterbury": "CAN",
        "Gisborne": "HBY",
        "Hawkes Bay": "HBY",
        "Hawke's Bay": "HBY",
        "Manawatu-Whanganui": "CEN",
        "Manawatū-Whanganui": "CEN",
        "Marlborough": "NEL",
        "Nelson": "NEL",
        "Northland": "NIS",
        "Otago": "OTG",
        "Southland": "OTG",
        "Taranaki": "TRN",
        "Tasman": "NEL",
        "Waikato": "WTO",
        "Wellington": "WEL",
        "West Coast": "CAN",
    }
    s = str(region_value).strip()
    if " " in s and s.split()[-1].isupper():
        # e.g. "Auckland AKL" → "AKL" [file:1]
        return s.split()[-1].strip()
    return region_to_remix.get(s, None)

def read_carriers_excel_long(excel_path: Path, scenario_name: str) -> pd.DataFrame:
    """Replicate the carriers-long reader from working_script.txt for one scenario."""  # [file:1]
    df = pd.read_excel(excel_path, sheet_name="carriers")
    df.columns = [c.strip() for c in df.columns]

    df = df.loc[df["Scenario"].astype(str).str.strip() == str(scenario_name)].copy()
    if df.empty:
        raise ValueError(f"No rows found in carriers Excel for scenario {scenario_name}")

    df["REMixRegion"] = df["Region"].apply(map_region_to_remix)
    missing = df.loc[df["REMixRegion"].isna(), "Region"].unique().tolist()
    if missing:
        raise ValueError("Unmapped regions in carriers Excel: " + ", ".join(str(x) for x in missing))

    year_cols = [c for c in df.columns if str(c).isdigit()]
    if not year_cols:
        raise ValueError("No year columns found in carriers sheet (expected 2020, 2025, ...).")

    dflong = df.melt(
        id_vars=["Scenario", "Sector", "Carrier", "Region", "REMixRegion"],
        value_vars=year_cols,
        var_name="Year",
        value_name="Demand",
    )
    dflong["Year"] = dflong["Year"].astype(int)
    dflong["Demand"] = pd.to_numeric(dflong["Demand"], errors="coerce").fillna(0.0)

    dflong = dflong.loc[(dflong["Demand"] != 0.0)].copy()
    return dflong


def build_annual_demand_from_excel(df_long: pd.DataFrame, year: int) -> pd.DataFrame:
    """
    Build annual demand upper bounds for fuels & e-fuels as in the helpers:
    index: (node, commodity), column 'upper'.
    Sector is aggregated because commodity_balance_annual is already aggregated. [file:1]
    """
    df = df_long.loc[df_long["Year"] == year].copy()
    if df.empty:
        print(f"No carriers Excel demands for year {year}.")
        return pd.DataFrame()

    df["CarrierKey"] = df["Carrier"].astype(str).str.strip().str.lower()

    paid_carrier_to_commodity = {
        "fossil (lf)": "LF_fossil",
        "fossil (gas)": "Gas_fossil",
        "coal": "Coal_fossil",
        "biofuel (lf)": "LF_bio",
        "biofuel (gas)": "Gas_bio",
        "wood": "Wood",
    }

    efuel_carrier_to_commodity = {
        "e-fuel (lf)": "REfuel",
        "e-fuel (gas)": "CH4",
    }

    # H2 (if you want to check it here)
    h2_carrier_to_commodity = {
        "hydrogen": "H2",
        "h2": "H2",
        "h2-feedstock": "H2",
    }


    df_paid = df.loc[df["CarrierKey"].isin(paid_carrier_to_commodity)].copy()
    df_paid["commodity"] = df_paid["CarrierKey"].map(paid_carrier_to_commodity)

    df_efuel = df.loc[df["CarrierKey"].isin(efuel_carrier_to_commodity)].copy()
    df_efuel["commodity"] = df_efuel["CarrierKey"].map(efuel_carrier_to_commodity)

    df_h2 = df.loc[df["CarrierKey"].isin(h2_carrier_to_commodity)].copy()
    df_h2["commodity"] = df_h2["CarrierKey"].map(h2_carrier_to_commodity)

    df_all = pd.concat([df_paid, df_efuel, df_h2], axis=0, ignore_index=True)


    agg = (
        df_all.groupby(["node", "commodity"], as_index=True)["Demand"]
        .sum()
        .to_frame(name="upper")
    )

    if CHECK_COMMODITIES:
        agg = agg.loc[agg.index.get_level_values("commodity").isin(CHECK_COMMODITIES)]


    return agg


# ----------------- GDX HELPERS ----------------------------------------

def load_gdx(path: Path) -> dict:
    data = gdxpds.to_dataframes(str(path))  # [web:3]
    print(f"Loaded {len(data)} tables from {path.name}")
    return data


def normalize_cols(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = df.columns.str.lower()
    return df


def extract_value_column(df: pd.DataFrame) -> pd.Series:
    cols = df.columns.tolist()
    if len(cols) >= 6 and cols[-5] == "level":
        return df["level"].astype(float)
    return df.iloc[:, -1].astype(float)


def prepare_commodity_balance_annual(data: dict, year: int) -> pd.DataFrame:
    """
    Use commodity_balance_annual to get net 2050 fuel use per (node, commodity).
    Net > 0 means net supply; for demand-type carriers (fuels at sinks),
    the sign convention in your build leads to positive flows for consumption
    when looking at the commodity side. [file:1]
    """
    if "commodity_balance_annual" not in data:
        raise KeyError("No table 'commodity_balance_annual' in GDX result.")

    df = data["commodity_balance_annual"]
    d = normalize_cols(df).reset_index()

    needed = ["accnodesmodel", "accyears", "commodity", "balancetype"]
    missing = [c for c in needed if c not in d.columns]
    if missing:
        raise KeyError(
            f"commodity_balance_annual missing columns {missing}. "
            f"Available: {list(d.columns)}"
        )

    # Local node balances, not global. [file:1]
    d = d.loc[
        (d["balancetype"] == "net")
        & (~d["accnodesmodel"].eq("global"))
        & (d["accyears"] == str(year))
    ].copy()
    if d.empty:
        print(f"No annual commodity balances for year {year}.")
        return pd.DataFrame()

    d = d.loc[d["commodity"].isin(CHECK_COMMODITIES)]
    if d.empty:
        print(f"No annual balances for CHECK_COMMODITIES in year {year}.")
        return pd.DataFrame()

    val = extract_value_column(d)

    d["used"] = val
    out = (
        d[["accnodesmodel", "commodity", "used"]]
        .rename(columns={"accnodesmodel": "node"})
        .groupby(["node", "commodity"], as_index=True)["used"]
        .sum()
        .to_frame()
    )
    return out


# ----------------- COMPARISON & PRINTING ------------------------------

def compare_and_print(upper_df: pd.DataFrame, used_df: pd.DataFrame, year: int):
    comp = upper_df.join(used_df, how="outer")
    comp = comp.fillna(0.0)
    comp["gap"] = comp["upper"] - comp["used"]
    comp["rel_use"] = np.where(comp["upper"] > 0, comp["used"] / comp["upper"], np.nan)

    if comp.empty:
        print("No comparison rows (no overlapping fuels between Excel and balances).")
        return

    print(f"\n=== 2050 Fuel/E-fuel Demand vs Use (Node, Commodity) ===")
    print(f"Rows: {len(comp)}")

    total_upper = comp["upper"].sum()
    total_used = comp["used"].sum()
    print(f"Total 2050 Excel upper (all checked fuels): {total_upper:,.3f} GWh")
    print(f"Total 2050 net use     (from balances):     {total_used:,.3f} GWh")
    if total_upper > 0:
        print(f"Overall utilization: {100 * total_used / total_upper:5.1f} %")

    overshoot = (comp["used"] > comp["upper"] + 1e-6).sum()
    exactish = (np.isclose(comp["used"], comp["upper"], atol=1e-3)).sum()
    unused = ((comp["upper"] > 0) & (comp["used"] == 0)).sum()
    print(f"Entries where use > upper (violations?): {overshoot}")
    print(f"Entries ~exactly using upper:           {exactish}")
    print(f"Entries with upper>0 but used=0:        {unused}")

    big_gap = comp.reindex(
        comp["gap"].abs().sort_values(ascending=False).index
    ).reset_index()

    print("\nTop 15 largest absolute gaps (upper - used) in 2050:")
    print(
        big_gap.head(15).to_string(
            index=False,
            formatters={
                "upper":   lambda x: f"{x:8.3f}",
                "used":    lambda x: f"{x:8.3f}",
                "gap":     lambda x: f"{x:+8.3f}",
                "rel_use": lambda x: "   n/a" if pd.isna(x) else f"{100*x:6.1f} %",
            },
        )
    )

    comp_reset = comp.reset_index()
    per_comm = (
        comp_reset.groupby("commodity")[["upper", "used"]]
        .sum()
        .assign(utilization=lambda df: np.where(df["upper"] > 0, df["used"] / df["upper"], np.nan))
        .sort_values("upper", ascending=False)
    )

    print("\nPer-commodity totals (2050):")
    print(
        per_comm.to_string(
            formatters={
                "upper":       lambda x: f"{x:8.3f}",
                "used":        lambda x: f"{x:8.3f}",
                "utilization": lambda x: "   n/a" if pd.isna(x) else f"{100*x:6.1f} %",
            }
        )
    )

    out_csv = GDX_PATH.with_name(f"{GDX_PATH.stem}_fuel_check_cb_{year}.csv")
    comp_reset.to_csv(out_csv, index=False)
    print(f"\nWrote detailed comparison to: {out_csv}")


# ----------------- MAIN -----------------------------------------------

def main():
    print("\n--- Checking 2050 fuel/e-fuel demand vs use (Excel vs commodity_balance_annual) ---")
    print(f"GDX:    {GDX_PATH}")
    print(f"Excel:  {CARRIERS_EXCEL}")
    print(f"Year:   {YEAR_CHECK}")

    carriers_long = read_carriers_excel_long(CARRIERS_EXCEL, BASE_SCENARIO)
    upper_df = build_annual_demand_from_excel(carriers_long, YEAR_CHECK)

    if upper_df.empty:
        print("No Excel-based annual demand found for checked fuels; nothing to compare.")
        return

    data = load_gdx(GDX_PATH)
    used_df = prepare_commodity_balance_annual(data, YEAR_CHECK)

    if used_df.empty:
        print("No commodity_balance_annual entries for checked fuels; nothing to compare.")
        return

    compare_and_print(upper_df, used_df, YEAR_CHECK)


if __name__ == "__main__":
    main()
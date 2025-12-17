from pathlib import Path
import pandas as pd
import gdxpds
import numpy as np

# --- USER SETTINGS ---
EXCEL_PATH = Path(
    "C:/Local/REMix/remix_nz/input/demand/GP-NT-ELEC-BIO-H2/data_summary.xlsx"
)
DATA_DIR = Path(
    "C:/Local/REMix/remix_nz/project/GP-NT-ELEC-BIO-H2/nz_case_GP_2050/data"
)
GDX_PATH = Path(
    "C:/Local/REMix/remix_nz/project/GP-NT-ELEC-BIO-H2/nz_case_GP_2050/result/nz_case_GP_2050.gdx"
)

SCENARIO = "GP"
YEARS_CHECK = [2050]

# Mapping (same as in your build)
REGION_MAP = {
    "Auckland": "AKL",
    "Bay of Plenty": "BOP",
    "Canterbury": "CAN",
    "Gisborne": "HBY",
    "Hawkes Bay": "HBY",
    "Hawke's Bay": "HBY",
    "Manawatu-Whanganui": "CEN",
    "Manawatū-Whanganui": "CEN",
    "Marlborough": "NEL",
    " Marlborough": "NEL",
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


CARRIER_MAP = {
    "biofuel (gas)": "Gas_bio", "biofuel (lf)": "LF_bio", "wood": "Wood",
    "coal": "Coal_fossil", "fossil (gas)": "Gas_fossil", "fossil (lf)": "LF_fossil",
    "hydrogen": "H2", "h2": "H2",
}


def load_excel_demand(path, scenario, years):
    """Load Excel and aggregate to (REMixRegion, Year, Carrier_Clean)."""
    print("[EXCEL] Loading Excel carriers sheet...")
    df = pd.read_excel(path, sheet_name="carriers")
    df.columns = [c.strip() for c in df.columns]
    df = df[df["Scenario"].astype(str).str.strip() == scenario].copy()

    df["REMixRegion"] = df["Region"].map(REGION_MAP)
    if df["REMixRegion"].isna().any():
        unmapped = df[df["REMixRegion"].isna()]["Region"].unique()
        print(f"  WARNING: Unmapped regions: {list(unmapped)}")
        df = df.dropna(subset=["REMixRegion"])

    df["Carrier_Clean"] = df["Carrier"].str.strip().str.lower().map(CARRIER_MAP)
    df = df.dropna(subset=["Carrier_Clean"])

    val_cols = [c for c in df.columns if str(c) in [str(y) for y in years]]
    long = df.melt(
        id_vars=["REMixRegion", "Carrier_Clean", "Sector"],
        value_vars=val_cols, var_name="Year", value_name="Demand",
    )
    long["Year"] = long["Year"].astype(int)
    agg = long.groupby(["REMixRegion", "Year", "Carrier_Clean"])["Demand"].sum()
    agg.name = "Excel"
    print(f"[EXCEL] ✓ Loaded. Total: {agg.sum():,.2f} GWh")
    return agg


def load_csv_demand(data_dir, years):
    """Load sourcesink_annualsum.csv (build input data)."""
    csv_path = data_dir / "sourcesink_annualsum.csv"
    if not csv_path.exists():
        print(f"[CSV] ERROR: {csv_path} not found.")
        return pd.Series()

    print(f"[CSV] Loading {csv_path}...")
    df = pd.read_csv(csv_path, index_col=[0, 1, 2, 3])
    df.index.names = ["nodesdata", "years", "sector", "commodity"]
    df = df.reset_index()

    df = df[df["years"].astype(str).isin([str(y) for y in years])].copy()

    def map_csv_comm(c):
        c = str(c).strip()
        if c == "H2": return "H2"
        if c.startswith("FuelService_"): return c.replace("FuelService_", "")
        return None

    df["Carrier_Clean"] = df["commodity"].apply(map_csv_comm)
    df = df.dropna(subset=["Carrier_Clean"])

    df["Demand"] = -df["upper"].astype(float)

    agg = df.groupby(["nodesdata", "years", "Carrier_Clean"])["Demand"].sum()
    agg.index.names = ["REMixRegion", "Year", "Carrier_Clean"]
    agg.name = "CSV"
    agg.index = agg.index.set_levels(
        pd.to_numeric(agg.index.levels[1], errors="coerce"), level="Year"
    )
    print(f"[CSV] ✓ Loaded. Total: {agg.sum():,.2f} GWh")
    return agg


def load_gdx_model_demand(path, years):
    """Load model DEMAND (FuelService consumption) from commodity_balance_annual."""
    print(f"[GDX] Loading {path}...")
    if not path.exists():
        print(f"[GDX] ERROR: {path} not found.")
        return pd.Series()

    try:
        dfs = gdxpds.to_dataframes(str(path))
    except Exception as e:
        print(f"[GDX] ERROR reading GDX: {e}")
        return pd.Series()

    # Look for commodity_balance_annual (model results)
    if "commodity_balance_annual" not in dfs:
        print(f"[GDX] ERROR: commodity_balance_annual not found in GDX.")
        print(f"[GDX]   Available tables: {list(dfs.keys())[:15]}")
        return pd.Series()

    cb = dfs["commodity_balance_annual"]
    cb.columns = [c.lower() for c in cb.columns]

    print(f"[GDX]   commodity_balance_annual shape: {cb.shape}")
    print(f"[GDX]   Columns: {list(cb.columns)}")

    # Filter for years
    cb = cb[cb["accyears"].astype(str).isin([str(y) for y in years])].copy()

    # We want negative balance (demand): FuelService_* and H2
    def map_gdx_comm(c):
        c = str(c).strip()
        if c == "H2": return "H2"
        if c.startswith("FuelService_"): return c.replace("FuelService_", "")
        return None

    cb["Carrier_Clean"] = cb["commodity"].apply(map_gdx_comm)
    cb = cb.dropna(subset=["Carrier_Clean"])

    # Filter for negative balance (demand side)
    cb = cb[cb["balancetype"] == "negative"].copy()

    # Value column
    numeric_cols = cb.select_dtypes(include=[np.number]).columns.tolist()
    value_col = numeric_cols[0] if numeric_cols else None
    if value_col is None:
        print("[GDX] ERROR: No numeric column found.")
        return pd.Series()

    # Demand is absolute value (negative balance means demand)
    cb["Demand"] = cb[value_col].astype(float).abs()

    # Aggregate
    agg = cb.groupby(["accnodesmodel", "accyears", "Carrier_Clean"])["Demand"].sum()
    agg.index.names = ["REMixRegion", "Year", "Carrier_Clean"]
    agg.name = "GDX_Model"

    agg.index = agg.index.set_levels(
        pd.to_numeric(agg.index.levels[1], errors="coerce"), level="Year"
    )

    print(f"[GDX] ✓ Loaded. Total model demand: {agg.sum():,.2f} GWh")
    return agg


def main():
    print("=" * 80)
    print("VERIFICATION: Excel → CSV (build) → GDX (model result)")
    print("=" * 80)

    print("\n--- Step 1: Load Excel ---")
    excel = load_excel_demand(EXCEL_PATH, SCENARIO, YEARS_CHECK)

    print("\n--- Step 2: Load CSV (build input) ---")
    csv = load_csv_demand(DATA_DIR, YEARS_CHECK)

    print("\n--- Step 3: Load GDX (model demand from commodity_balance_annual) ---")
    gdx = load_gdx_model_demand(GDX_PATH, YEARS_CHECK)

    # Combine
    print("\n" + "=" * 80)
    print("COMPARISON")
    print("=" * 80 + "\n")

    df_compare = pd.concat([excel, csv, gdx], axis=1).fillna(0.0)

    df_compare["Excel-CSV"] = df_compare["Excel"] - df_compare["CSV"]
    df_compare["CSV-GDX"] = df_compare["CSV"] - df_compare["GDX_Model"]
    df_compare["Excel-GDX"] = df_compare["Excel"] - df_compare["GDX_Model"]

    tol = 1e-2
    df_compare["Excel=CSV"] = df_compare["Excel-CSV"].abs() < tol
    df_compare["CSV=GDX"] = df_compare["CSV-GDX"].abs() < tol

    print(df_compare.to_string())

    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    tot_excel = df_compare["Excel"].sum()
    tot_csv = df_compare["CSV"].sum()
    tot_gdx = df_compare["GDX_Model"].sum()

    print(f"\nTotal Demands:")
    print(f"  Excel:              {tot_excel:>15,.2f} GWh")
    print(f"  CSV (build input):  {tot_csv:>15,.2f} GWh")
    print(f"  GDX (model result): {tot_gdx:>15,.2f} GWh")

    print(f"\nAbsolute Differences:")
    print(f"  Excel - CSV:        {(tot_excel - tot_csv):>15,.2f} GWh")
    print(f"  CSV - GDX:          {(tot_csv - tot_gdx):>15,.2f} GWh")
    print(f"  Excel - GDX:        {(tot_excel - tot_gdx):>15,.2f} GWh")

    print(f"\nMatch Count (within {tol} GWh):")
    print(f"  Excel = CSV:        {df_compare['Excel=CSV'].sum()} / {len(df_compare)}")
    print(f"  CSV = GDX:          {df_compare['CSV=GDX'].sum()} / {len(df_compare)}")

    print("\n" + "=" * 80)
    print("MISMATCHES")
    print("=" * 80)

    ex_csv = df_compare[~df_compare["Excel=CSV"]]
    if ex_csv.empty:
        print("\n✓ Excel → CSV: Perfect match!")
    else:
        print("\n✗ Excel → CSV mismatches:")
        print(ex_csv[["Excel", "CSV", "Excel-CSV"]].to_string())

    csv_gdx = df_compare[~df_compare["CSV=GDX"]]
    if csv_gdx.empty:
        print("\n✓ CSV → GDX: Model fully used all demands!")
    else:
        print("\n✗ CSV → GDX (model may not use all demands):")
        print(csv_gdx[["CSV", "GDX_Model", "CSV-GDX"]].to_string())

    print("\n" + "=" * 80)


if __name__ == "__main__":
    main()
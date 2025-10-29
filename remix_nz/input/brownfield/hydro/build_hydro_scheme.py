"""
Process JADE hydrological inflow data for REMix.
This script:
  1. Reads the JADE weekly inflow file (cumecs)
  2. Maps each lake or inflow catchment to the corresponding REMix-NZ commodity
  3. Analyses annual inflows to identify dry, wet and average years
  4. Builds hourly inflow profiles for selected years
  5. Exports them to REMix demand format (node, year, sector, commodity, t0001..t8760)
"""

import pandas as pd
import numpy as np
import os
import re


# File paths and basic settings

# File paths and basic settings

INFLOWS_CSV = r"C:/Local/REMix/remix_nz/input/brownfield/hydro/JADE_hydro/inflows.csv"
OUTPUT_FOLDER = r"C:/Local/REMix/remix_nz/input/brownfield/hydro/inflows_remix-nz"

# Output filename automatically includes year mapping
SOURCE_YEARS = [2012, 2014]  # years from inflows.csv
TARGET_YEARS = [2020, 2030]  # names to use in REMix
YEAR_OF_INTEREST = 2022      # for hydrological analysis

output_tags = "_".join([f"{s}-to-{t}" for s, t in zip(SOURCE_YEARS, TARGET_YEARS)])
OUTPUT_FILE = f"inflow_{output_tags}.csv"


# Select which years to extract and what to rename them to
SOURCE_YEARS = [2012, 2014]  # years from inflows.csv
TARGET_YEARS = [2020, 2030]  # names to use in REMix
YEAR_OF_INTEREST = 2022      # year to analyse in detail / see how average it is


HOURS_PER_WEEK = 168
HOURS_PER_YEAR = 8760
SECTOR = "Inflow"          # sector name for REMix
NEGATIVE_VALUES = True     # inflows are negative so they become positive in add_demand()


# Mapping between JADE lakes and REMix-NZ commodities / turbines / regions
COMMODITY_MAP = {
    "Lake_Tekapo": ("Tekapo_in", "CAN", "Tekapo_A", 0.2322),
    "Lake_Pukaki": ("Pukaki_in", "CAN", "Ohau_A", 0.5005),
    "Lake_Ohau": ("Ohau_in", "CAN", "Ohau_A", 0.5005),
    "Lake_Benmore": ("Benmore_in", "CAN", "Benmore", 0.8177),
    "Lake_Aviemore": ("Aviemore_in", "CAN", "Aviemore", 0.3101),
    "Lake_Waitaki": ("Waitaki_in", "CAN", "Waitaki", 0.1622),
    "Lake_Coleridge": ("Coleridge_in", "CAN", "Coleridge", 1.009),
    "Lake_Taupo": ("Aratiatia_in", "WTO", "Aratiatia", 0.2841),
    "Lake_Ohakuri": ("Ohakuri_in", "WTO", "Ohakuri", 0.2841),
    "Lake_Atiamuri": ("Atiamuri_in", "WTO", "Ātiamuri", 0.1957),
    "Lake_Whakamaru": ("Whakamaru_in", "WTO", "Whakamaru", 0.3165),
    "Lake_Maraetai": ("Maraetai_in", "WTO", "Maraetai", 0.5263),
    "Lake_Waipapa": ("Waipapa_in", "WTO", "Waipapa", 0.1385),
    "Lake_Arapuni": ("Arapuni_in", "WTO", "Arapuni", 0.4619),
    "Lake_Karapiro": ("Karapiro_in", "WTO", "Karapiro", 0.2639),
    "Lake_Hawea": ("Dunstan_in", "OTG", "Clyde_220kV", 0.5181),
    "Lake_Wanaka": ("Dunstan_in", "OTG", "Clyde_220kV", 0.5181),
    "Lake_Dunstan": ("Dunstan_in", "OTG", "Clyde_220kV", 0.5181),
    "Lake_Roxburgh": ("Roxburgh_in", "OTG", "Roxburgh", 0.4016),
    "Lakes_Manapouri_Te_Anau": ("Manapouri_in", "OTG", "Manapouri", 1.5314),
    "Lake_Waikaremoana": ("Waikaremoana_in", "HBY", "Waikaremoana", 3.535),
    "Lake_Rotoaira": ("Tokaanu_in", "CEN", "Tokaanu", 1.75),
    "Lake_Moawhango": ("Rangipo_in", "CEN", "Rangipo", 1.96),
    "Lake_Matahina": ("Matahina_in", "CEN", "Matahina", 0.595),
    "Mangahao_head": ("Mangahao_in", "CEN", "Mangahao", 2.53),
    "Lake_Cobb": ("Cobb_in", "NEL", "Cobb", 4.405),
}


def normalize_name(text):
    """Clean column names for consistent matching."""
    t = str(text).strip().replace(" ", "_").replace("-", "_")
    return re.sub(r"__+", "_", t)

def read_jade_inflows(path):
    """
    Read JADE inflows.csv correctly for files structured as:

    % comment
    <blank>
    CATCHMENT,,Lake_Arapuni,Lake_Aratiatia,...
    INFLOW_REGION,,NI,NI,NI,SI,SI,...
    YEAR,WEEK,<data...>
    1932,1,<numbers>...

    Returns:
        DataFrame with numeric YEAR, WEEK, and inflows for each lake.
    """
    import pandas as pd
    from io import StringIO

    # Read all lines, skipping comments and blank lines
    with open(path, "r", encoding="utf-8-sig") as f:
        lines = [line.strip() for line in f if line.strip() and not line.startswith("%")]

    # Find the important rows
    catch_line = next(i for i, l in enumerate(lines) if l.startswith("CATCHMENT"))
    region_line = next(i for i, l in enumerate(lines) if l.startswith("INFLOW_REGION"))
    data_start = region_line + 2  # skip region line + one "YEAR,WEEK" line

    # Extract headers from CATCHMENT row
    catchments = [c.strip() for c in lines[catch_line].split(",")]
    # Fix double commas (blank columns)
    catchments = [c if c else None for c in catchments]

    # The data starts two lines after the INFLOW_REGION line
    data_text = "\n".join(lines[data_start:])

    # Read data
    df = pd.read_csv(StringIO(data_text), header=None)

    # Assign headers: first two are YEAR, WEEK; rest come from catchments[2:]
    headers = ["YEAR", "WEEK"] + [normalize_name(c) for c in catchments[2:]]
    df.columns = headers

    # Clean numeric columns
    df["YEAR"] = pd.to_numeric(df["YEAR"], errors="coerce")
    df["WEEK"] = pd.to_numeric(df["WEEK"], errors="coerce")

    df = df.dropna(subset=["YEAR", "WEEK"]).reset_index(drop=True)

    print(f"Loaded inflow file: {len(df)} rows, {len(df.columns)} columns")
    print(f"Years range: {int(df['YEAR'].min())}–{int(df['YEAR'].max())}")
    return df


def years_with_full_weeks(df):
    """Return only years with at least 52 weeks of data."""
    weeks = df.groupby("YEAR")["WEEK"].nunique()
    return weeks[weeks >= 52].index

def weekly_to_hourly(weekly_values):
    """Repeat each weekly value for 168 hours (1 week)."""
    hourly = np.repeat(weekly_values, HOURS_PER_WEEK)
    if hourly.size < HOURS_PER_YEAR:
        hourly = np.append(hourly, np.repeat(weekly_values[-1], HOURS_PER_YEAR - hourly.size))
    return hourly

def flow_to_power(q_cumecs, specific_power):
    """Convert average flow (cumecs) to power (GW) using the turbine's specific power (MW/cumec)."""
    return (specific_power / 1000.0) * q_cumecs

def analyse_water_years(df):
    """Print driest, wettest and average water years."""
    lakes = [normalize_name(k) for k in COMMODITY_MAP.keys() if normalize_name(k) in df.columns]
    df = df[df["YEAR"].isin(years_with_full_weeks(df))]
    annual = df.groupby("YEAR")[lakes].sum().sum(axis=1)
    mean = annual.mean()
    driest = annual.nsmallest(5)
    wettest = annual.nlargest(5)
    avg_like = (annual - mean).abs().nsmallest(5)

    print(f"\nAverage annual inflow (sum over mapped lakes): {mean:.2f}")
    print("\nFive driest years:")
    for y, v in driest.items(): print(f"  {int(y)} → {v:.2f}")
    print("\nFive wettest years:")
    for y, v in wettest.items(): print(f"  {int(y)} → {v:.2f}")
    print("\nFive most average years:")
    for y, v in avg_like.items(): print(f"  {int(y)} → {annual.loc[y]:.2f} (Δ={v:.2f})")

    # Recent decade statistics (2012–2022)
    # Compute annual inflow totals across all mapped catchments 
    # Sum weekly inflows for each year, across all catchments that were successfully mapped
    inflow_cols = [col for col in df.columns if col not in ("YEAR", "WEEK")]
    yearly_totals = (
        df.groupby("YEAR")[inflow_cols].sum().sum(axis=1).reset_index()
    )
    yearly_totals.columns = ["YEAR", "TOTAL"]

    print(f"\nAverage annual inflow (sum over mapped lakes): {yearly_totals['TOTAL'].mean():.2f}")

    recent_years = list(range(2012, 2023))
    recent = yearly_totals.loc[yearly_totals["YEAR"].isin(recent_years)]

    long_term_mean = yearly_totals["TOTAL"].mean()
    yearly_totals["rel_to_mean_%"] = (yearly_totals["TOTAL"] / long_term_mean - 1) * 100
    # --- Year-to-year change in inflow (% difference from previous year) ---
    yearly_totals["change_from_prev_%"] = yearly_totals["TOTAL"].pct_change() * 100
    # --- Rolling 10-year mean comparison ---
    yearly_totals["rolling10_mean"] = yearly_totals["TOTAL"].rolling(window=10, min_periods=5).mean()
    yearly_totals["rel_to_rolling10_%"] = (yearly_totals["TOTAL"] / yearly_totals["rolling10_mean"] - 1) * 100



    if not recent.empty:
        avg_recent = recent["TOTAL"].mean()
        min_recent = recent["TOTAL"].min()
        max_recent = recent["TOTAL"].max()
        driest_recent = recent.loc[recent["TOTAL"].idxmin()]
        wettest_recent = recent.loc[recent["TOTAL"].idxmax()]

        print("\nHydrological summary for the last decade (2012–2022):")
        print(f"  Average inflow across years: {avg_recent:,.2f}")
        print(f"  Driest year: {int(driest_recent['YEAR'])} → {driest_recent['TOTAL']:,.2f}")
        print(f"  Wettest year: {int(wettest_recent['YEAR'])} → {wettest_recent['TOTAL']:,.2f}")
        print(f"  Range (wet–dry): {max_recent - min_recent:,.2f}")
    else:
        print("\nNo data found for 2012–2022 in inflows.csv.")

def describe_hydro_year(df, year):
    """
    Analyse inflows for a given year in table form.
    Includes annual, winter (Jun–Sep), and non-winter comparisons
    across long-term, 30-year, 10-year, previous, and next periods.
    """

    inflow_cols = [c for c in df.columns if c not in ("YEAR", "WEEK")]
    df["PERIOD"] = np.where(df["WEEK"].between(22, 39), "Winter", "NonWinter")

    # Aggregate inflows by year and season
    annual = df.groupby("YEAR")[inflow_cols].sum().sum(axis=1).rename("Annual")
    winter = df.loc[df["PERIOD"] == "Winter"].groupby("YEAR")[inflow_cols].sum().sum(axis=1).rename("Winter")
    nonwinter = df.loc[df["PERIOD"] == "NonWinter"].groupby("YEAR")[inflow_cols].sum().sum(axis=1).rename("NonWinter")

    yearly_totals = pd.concat([annual, winter, nonwinter], axis=1).reset_index()
    yearly_totals["Winter_share"] = yearly_totals["Winter"] / yearly_totals["Annual"]

    if year not in yearly_totals["YEAR"].values:
        print(f"Year {year} not found in inflow records.")
        return

    def pct_diff(value, ref):
        return (value / ref - 1) * 100 if ref != 0 else np.nan

    # Reference periods for comparison
    ranges = {
        "All years": yearly_totals["YEAR"].between(yearly_totals["YEAR"].min(), yearly_totals["YEAR"].max()),
        "1992–2022": yearly_totals["YEAR"].between(1992, 2022),
        "2012–2022": yearly_totals["YEAR"].between(2012, 2022),
        "Previous year": yearly_totals["YEAR"] == year - 1,
        "Next year": yearly_totals["YEAR"] == year + 1,
    }

    # Summary table
    data, labels = [], []
    row = yearly_totals.set_index("YEAR").loc[year]

    for label, mask in ranges.items():
        ref = yearly_totals.loc[mask]
        if ref.empty:
            continue
        ref_row = ref.mean()
        data.append([
            pct_diff(row["Annual"], ref_row["Annual"]),
            pct_diff(row["Winter"], ref_row["Winter"]),
            pct_diff(row["NonWinter"], ref_row["NonWinter"])
        ])
        labels.append(label)

    summary = pd.DataFrame(data, columns=["Annual Δ%", "Winter Δ%", "Non-winter Δ%"], index=labels)
    print(f"\nHydrological comparison for {year}")
    print(summary.round(2).to_string())

    # Helper for concise comments
    def seasonal_comment(row, ref_df):
        wdiff = pct_diff(row["Winter"], ref_df["Winter"].mean())
        nwdiff = pct_diff(row["NonWinter"], ref_df["NonWinter"].mean())
        wsign = "above" if wdiff > 0 else "below"
        nwsign = "above" if nwdiff > 0 else "below"
        return f"(Winter inflows were {abs(wdiff):.1f}% {wsign} average; rest of year {abs(nwdiff):.1f}% {nwsign} average)"

    # Reference year summaries by period
    def print_reference_years(title, mask):
        ref_df = yearly_totals.loc[mask].copy()
        if ref_df.empty:
            return
        mean_annual = ref_df["Annual"].mean()
        ref_df["DiffFromMean"] = (ref_df["Annual"] - mean_annual).abs()
        driest = ref_df.loc[ref_df["Annual"].idxmin()]
        wettest = ref_df.loc[ref_df["Annual"].idxmax()]
        avg = ref_df.loc[ref_df["DiffFromMean"].idxmin()]

        print(f"\nReference years for {title}:")
        print(f"  Driest year:  {int(driest['YEAR'])}  → {driest['Annual']:,.1f}  {seasonal_comment(driest, ref_df)}")
        print(f"  Wettest year: {int(wettest['YEAR'])} → {wettest['Annual']:,.1f}  {seasonal_comment(wettest, ref_df)}")
        print(f"  Most average: {int(avg['YEAR'])}    → {avg['Annual']:,.1f}  {seasonal_comment(avg, ref_df)}")

    print_reference_years("entire record", yearly_totals["YEAR"].between(yearly_totals["YEAR"].min(), yearly_totals["YEAR"].max()))
    print_reference_years("1992–2022", yearly_totals["YEAR"].between(1992, 2022))
    print_reference_years("2012–2022", yearly_totals["YEAR"].between(2012, 2022))

    print("--------------------------------------------------")


def build_inflow_profiles(df, source_years, target_years):
    """Build hourly inflow profiles for each selected year and rename them."""
    if len(source_years) != len(target_years):
        raise ValueError("Source and target year lists must have the same length.")

    df = df[df["YEAR"].isin(years_with_full_weeks(df))]
    lake_cols = [normalize_name(k) for k in COMMODITY_MAP.keys() if normalize_name(k) in df.columns]

    results = []
    for src, tgt in zip(source_years, target_years):
        df_y = df[df["YEAR"] == src].sort_values("WEEK")
        rows, data = [], []

        for lake_col in lake_cols:
            lake_name = [k for k in COMMODITY_MAP.keys() if normalize_name(k) == lake_col][0]
            commodity, region, turbine, sp = COMMODITY_MAP[lake_name]

            weekly_q = df_y[lake_col].to_numpy(float)
            weekly_gw = flow_to_power(weekly_q, sp)
            hourly_gw = weekly_to_hourly(weekly_gw)
            if NEGATIVE_VALUES:
                hourly_gw *= -1

            rows.append((region, tgt, SECTOR, commodity))
            data.append(hourly_gw)

        tcols = [f"t{t:04d}" for t in range(1, HOURS_PER_YEAR + 1)]
        df_out = pd.DataFrame(
            np.vstack(data),
            index=pd.MultiIndex.from_tuples(rows, names=["node", "year", "sector", "commodity"]),
            columns=tcols
        )
        results.append(df_out)

    inflow_profiles = pd.concat(results).sort_index()
    return inflow_profiles





if __name__ == "__main__":

    print("\n--- Processing JADE inflow data ---")
    df = read_jade_inflows(INFLOWS_CSV)
    print(f"Loaded inflow data: {df.shape[0]} rows, {df.shape[1]} columns.")

    # Check and print mappings
    available = [c for c in df.columns if c not in ("YEAR", "WEEK")]
    print("\n--- Mapping inflow catchments to REMix commodities ---")
    matched, missing = [], []
    for lake in COMMODITY_MAP.keys():
        lake_n = normalize_name(lake)
        if lake_n in available:
            c, r, t, sp = COMMODITY_MAP[lake]
            print(f"  {lake} → {c}  (region {r}, turbine {t}, SP={sp} MW/cumec)")
            matched.append(lake)
        else:
            missing.append(lake)
    if missing:
        print("\nWarning: the following catchments were not found in inflows.csv:")
        for m in missing: print("  -", m)

    # Analyse and display annual inflow characteristics
    analyse_water_years(df)

    # Build hourly inflow profiles for the selected years
    print("\n--- Building inflow profiles ---")
    inflows_hourly = build_inflow_profiles(df, SOURCE_YEARS, TARGET_YEARS)

    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    output_path = os.path.join(OUTPUT_FOLDER, OUTPUT_FILE)
    inflows_hourly.to_csv(output_path, float_format="%.9f")

    print(f"\nHydro inflow profiles saved to:\n  {output_path}")
    print(f"Rows: {inflows_hourly.shape[0]} | Columns: {inflows_hourly.shape[1]}")


    describe_hydro_year(df, YEAR_OF_INTEREST)
    


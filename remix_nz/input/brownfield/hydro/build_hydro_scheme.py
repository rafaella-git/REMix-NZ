"""
Process JADE hydrological inflow data for REMix.

Creates:
  • Hourly inflows in both cumecs (water) and GW (energy)
  • Full hydrological summaries for selected SOURCE_YEARS
  • Historical stats (driest, wettest, average, decade range)
"""

import pandas as pd
import numpy as np
import os, re
from io import StringIO


INFLOWS_CSV = r"C:/Local/REMix/remix_nz/input/brownfield/hydro/JADE_hydro/inflows.csv"
OUTPUT_FOLDER = r"C:/Local/REMix/remix_nz/input/brownfield/hydro/inflows_remix-nz"
SOURCE_YEARS = [2021, 2012]   # JADE years to extract
TARGET_YEARS = [2020, 2030]   # Names to assign for REMix
HOURS_PER_WEEK = 168
HOURS_PER_YEAR = 8760
SECTOR = "Inflow"
NEGATIVE_VALUES = True

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
    "Lake_Atiamuri": ("Atiamuri_in", "WTO", "Atiamuri", 0.1957),
    "Lake_Whakamaru": ("Whakamaru_in", "WTO", "Whakamaru", 0.3165),
    "Lake_Maraetai": ("Maraetai_in", "WTO", "Maraetai", 0.5263),
    "Lake_Waipapa": ("Waipapa_in", "WTO", "Waipapa", 0.1385),
    "Lake_Arapuni": ("Arapuni_in", "WTO", "Arapuni", 0.4619),
    "Lake_Karapiro": ("Karapiro_in", "WTO", "Karapiro", 0.2639),
    "Lake_Hawea": ("Dunstan_in", "OTG", "Clyde_220kV", 0.5181),
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
    return re.sub(r"__+", "_", str(text).strip().replace(" ", "_").replace("-", "_"))

def read_jade_inflows(path):
    """Read JADE inflows.csv correctly."""
    with open(path, "r", encoding="utf-8-sig") as f:
        lines = [l.strip() for l in f if l.strip() and not l.startswith("%")]
    catch_line = next(i for i, l in enumerate(lines) if l.startswith("CATCHMENT"))
    region_line = next(i for i, l in enumerate(lines) if l.startswith("INFLOW_REGION"))
    data_start = region_line + 2
    cols = [normalize_name(c) for c in lines[catch_line].split(",")[2:]]
    df = pd.read_csv(StringIO("\n".join(lines[data_start:])), header=None)
    df.columns = ["YEAR", "WEEK"] + cols
    df = df.dropna(subset=["YEAR", "WEEK"]).reset_index(drop=True)
    print(f"Loaded inflow data: {df.shape[0]} rows, {df.shape[1]} columns.")
    return df

def years_with_full_weeks(df):
    w = df.groupby("YEAR")["WEEK"].nunique()
    return w[w >= 52].index

def weekly_to_hourly(v):
    h = np.repeat(v, HOURS_PER_WEEK)
    if len(h) < HOURS_PER_YEAR:
        h = np.append(h, np.repeat(v[-1], HOURS_PER_YEAR - len(h)))
    return h

def flow_to_power(flow, sp):
    """Convert water flow (cumecs) to electrical power (GW)."""
    return (sp / 1000.0) * flow

def analyse_water_years(df):
    """Print driest, wettest, and average years."""
    lakes = [normalize_name(k) for k in COMMODITY_MAP if normalize_name(k) in df.columns]
    df = df[df["YEAR"].isin(years_with_full_weeks(df))]
    annual = df.groupby("YEAR")[lakes].sum().sum(axis=1)
    mean = annual.mean()
    driest = annual.nsmallest(5)
    wettest = annual.nlargest(5)
    avg_like = (annual - mean).abs().nsmallest(5)

    print("\n--- Overall inflow statistics ---")
    print(f"Average annual inflow (sum over all mapped lakes): {mean:,.2f}")
    print("\nFive driest years:")
    for y,v in driest.items(): print(f"  {int(y)} → {v:,.2f}")
    print("\nFive wettest years:")
    for y,v in wettest.items(): print(f"  {int(y)} → {v:,.2f}")
    print("\nFive most average years:")
    for y,v in avg_like.items(): print(f"  {int(y)} → {annual.loc[y]:,.2f} (Δ={v:,.2f})")

    inflow_cols = [c for c in df.columns if c not in ("YEAR","WEEK")]
    yearly_totals = df.groupby("YEAR")[inflow_cols].sum().sum(axis=1).reset_index()
    yearly_totals.columns = ["YEAR","TOTAL"]
    recent = yearly_totals[yearly_totals["YEAR"].between(2012,2022)]
    if not recent.empty:
        avg_recent = recent["TOTAL"].mean()
        driest_recent = recent.loc[recent["TOTAL"].idxmin()]
        wettest_recent = recent.loc[recent["TOTAL"].idxmax()]
        print("\n--- Hydrological summary for 2012–2022 ---")
        print(f"Average inflow: {avg_recent:,.2f}")
        print(f"Driest year: {int(driest_recent['YEAR'])} → {driest_recent['TOTAL']:,.2f}")
        print(f"Wettest year: {int(wettest_recent['YEAR'])} → {wettest_recent['TOTAL']:,.2f}")
        print(f"Range (wet–dry): {wettest_recent['TOTAL'] - driest_recent['TOTAL']:,.2f}")


def print_hydro_summary(df, year):
    """Print annual, winter, and non-winter flow summaries with Δ and energy equiv."""
    df = df[df["YEAR"].isin(years_with_full_weeks(df))].copy()
    cols = [normalize_name(k) for k in COMMODITY_MAP if normalize_name(k) in df.columns]
    df["PERIOD"] = np.where(df["WEEK"].between(22,39),"Winter","NonWinter")

    annual = df.groupby("YEAR")[cols].sum().sum(axis=1)
    winter = df[df["PERIOD"]=="Winter"].groupby("YEAR")[cols].sum().sum(axis=1)
    nonwinter = df[df["PERIOD"]=="NonWinter"].groupby("YEAR")[cols].sum().sum(axis=1)
    data = pd.DataFrame({"Annual":annual,"Winter":winter,"NonWinter":nonwinter}).sort_index()

    if year not in data.index:
        print(f"Year {year} not found.")
        return

    ref10 = data.loc[(data.index>=year-10)&(data.index<year)].mean()
    ref30 = data.loc[(data.index>=year-30)&(data.index<year)].mean()
    row = data.loc[year]
    pct = lambda v,r: (v/r - 1)*100 if r else np.nan

    avg_sp = np.mean([v[3] for v in COMMODITY_MAP.values()])
    energy = row * avg_sp * HOURS_PER_WEEK * 52 / 1000

    print(f"\nHydrological summary for {year}")
    print("-----------------------------------------------------------------------------")
    print(f"{'Metric':<12}{'Flow (cumecs)':>18}{'Δ10yr %':>10}{'Δ30yr %':>10}{'Energy-equiv (GWh)':>20}")
    print("-----------------------------------------------------------------------------")
    for k in data.columns:
        print(f"{k:<12}{row[k]:>18.1f}{pct(row[k],ref10[k]):>10.2f}{pct(row[k],ref30[k]):>10.2f}{energy[k]:>20,.1f}")
    print("-----------------------------------------------------------------------------")
    wet = "wetter" if pct(row["Annual"],ref10["Annual"])>0 else "drier"
    print(f"{year} was {abs(pct(row['Annual'],ref10['Annual'])):.1f}% {wet} than the 10-year mean.\n")


def build_inflow_profiles(df, src, tgt, energy=False):
    df = df[df["YEAR"].isin(years_with_full_weeks(df))]
    lakes = [normalize_name(k) for k in COMMODITY_MAP if normalize_name(k) in df.columns]
    results = []
    for s,t in zip(src,tgt):
        df_y = df[df["YEAR"]==s].sort_values("WEEK")
        rows, data = [], []
        for lake in lakes:
            lake_name = [k for k in COMMODITY_MAP if normalize_name(k)==lake][0]
            commodity, region, _, sp = COMMODITY_MAP[lake_name]
            weekly = df_y[lake].to_numpy(float)
            hourly = weekly_to_hourly(weekly)
            if energy: hourly = flow_to_power(hourly, sp)
            if NEGATIVE_VALUES: hourly *= -1
            rows.append((region,t,SECTOR,commodity))
            data.append(hourly)
        tcols = [f"t{t:04d}" for t in range(1,HOURS_PER_YEAR+1)]
        results.append(pd.DataFrame(np.vstack(data),
            index=pd.MultiIndex.from_tuples(rows,names=["node","year","sector","commodity"]),
            columns=tcols))
    return pd.concat(results).sort_index()

if __name__ == "__main__":
    print("\n--- Processing JADE inflow data ---")
    df = read_jade_inflows(INFLOWS_CSV)
    analyse_water_years(df)

    inflow_cumecs = build_inflow_profiles(df, SOURCE_YEARS, TARGET_YEARS, energy=False)
    inflow_energy = build_inflow_profiles(df, SOURCE_YEARS, TARGET_YEARS, energy=True)

    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    base = "_".join([f"{s}-to-{t}" for s,t in zip(SOURCE_YEARS,TARGET_YEARS)])
    f1 = os.path.join(OUTPUT_FOLDER, f"inflow_{base}_INFLOWS_cumecs.csv")
    f2 = os.path.join(OUTPUT_FOLDER, f"inflow_{base}_ENERGY_GW.csv")
    inflow_cumecs.to_csv(f1, float_format="%.9f")
    inflow_energy.to_csv(f2, float_format="%.9f")

    print(f"\nHourly inflow profiles saved to:\n  {f1}\n  {f2}")
    print(f"Rows: {inflow_cumecs.shape[0]} | Columns: {inflow_cumecs.shape[1]}")

    # Print hydrological summaries for both extracted years
    for y in SOURCE_YEARS:
        print_hydro_summary(df, y)

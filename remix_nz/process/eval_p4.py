# compare_results_paper4.py
# Multi-scenario comparison plotting + printing for REMix NZ (Paper 4)
#
# Run from VSCode (Python), not bash.
# - Loads multiple scenario result GDX files using remix.framework.tools.gdx.GDXEval
# - Robust to different year selections and different tech sets across scenarios
# - Prints schema/debug to console so you can show me what your GDX actually contains
# - Writes CSV and Excel tables
# - Saves PNG + SVG plots into comparison folder
#
# Notes:
# - Uses only matplotlib (no seaborn) for maximum portability
# - Sign convention:
#   commodity_balance: positive = production, negative = consumption (kept negative)
#   For "demand plots" we plot abs(negative) as positive bars
#
# Author: ChatGPT (tailored to your REMix NZ project structure)

from __future__ import annotations

import os
import sys
import math
import traceback
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Sequence

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from remix.framework.tools.gdx import GDXEval

# Optional (maps). If missing, script will continue without maps.
try:
    import geopandas as gpd  # noqa
    from remix.framework.tools.plots import plot_network, plot_choropleth  # noqa
    HAS_GEO = True
except Exception:
    HAS_GEO = False


# -----------------------------
# User settings (EDIT HERE)
# -----------------------------
BASE_PROJECT = Path(r"C:\Local\REMix\remix_nz\project\GP-NT-ELEC-BIO-H2")

SCENARIOS = [
    "nz_case_GP_2020-2050",
    "nz_case_NT_2020-2050",
    "nz_case_ELEC+_2020-2050",
    "nz_case_BIO+_2020-2050",
    "nz_case_H2+_2020-2050",
]

# Output folder (auto-created)
OUT_DIR = BASE_PROJECT / "comparison"

# Canonical strings
GLOBAL_NODE = "global"
ELEC_COMMODITY = "Elec"

# Thresholds / robustness
CAP_NEG_TOL = 0.01  # GW: if negative capacity abs > this, print debug
CAP_POS_TOL = 0.0001  # GW: plot only capacities > this (after fixing tiny negatives)
ENERGY_TOL = 0.01  # GWh: filter tiny annual/hours for cleanliness

# Time slices (hour index based, assuming timeModel is 1..8760 or similar sortable)
# We will use "hour numbers" after sorting timeModel.
WINTER_CENTER_HOUR = 500   # adjust if you want
SUMMER_CENTER_HOUR = 4500  # adjust if you want
SLICE_WEEKS = 2
HOURS_PER_WEEK = 168
HOURS_PER_SLICE = SLICE_WEEKS * HOURS_PER_WEEK  # 336

# If you want to compare a specific model year across hourly plots, keep None to auto-pick max year per scenario,
# or set e.g. "2050".
HOURLY_YEAR_OVERRIDE: Optional[str] = None

# Plot appearance
plt.rcParams.update({"figure.autolayout": True})
plt.rcParams.update({"savefig.dpi": 300})
plt.rcParams.update({"figure.dpi": 120})
plt.rcParams.update({"font.size": 11})
plt.rcParams.update({"axes.titlesize": 12})
plt.rcParams.update({"axes.labelsize": 11})


# -----------------------------
# Helpers
# -----------------------------
def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def safe_get_symbol(results: GDXEval, name: str) -> Optional[pd.Series]:
    """Try to fetch a symbol; return None if missing."""
    try:
        x = results[name]
        # Many REMix symbols return a Series-like object with MultiIndex and a value column,
        # but in practice it's usually a pandas Series with name 'value'.
        return x
    except Exception:
        return None


def series_to_df(s: pd.Series) -> pd.DataFrame:
    """Normalize to DataFrame with a 'value' column and preserved index levels."""
    if s is None:
        return pd.DataFrame()
    if isinstance(s, pd.DataFrame):
        if "value" in s.columns:
            return s.copy()
        # If it's already a df but no value column, try first column
        df = s.copy()
        df = df.rename(columns={df.columns[0]: "value"})
        return df
    # Series
    df = s.to_frame(name="value")
    return df


def print_schema(name: str, s: Optional[pd.Series], max_levels_preview: int = 12) -> None:
    if s is None:
        print(f"  - {name}: MISSING")
        return
    idx = s.index
    print(f"  - {name}: type={type(s).__name__}  rows={len(s):,}")
    if hasattr(idx, "names"):
        print(f"    index.names = {idx.names}")
        # Show unique counts (capped)
        for lvl in idx.names:
            try:
                nuniq = idx.get_level_values(lvl).nunique()
                print(f"    - {lvl}: n_unique={nuniq}")
            except Exception:
                pass
        # Preview a few items per level
        for lvl in idx.names[:max_levels_preview]:
            try:
                vals = list(pd.Index(idx.get_level_values(lvl)).unique()[:8])
                print(f"    - {lvl} sample: {vals}")
            except Exception:
                pass
    else:
        print("    (no MultiIndex names found)")


def infer_years(s: pd.Series) -> List[str]:
    """Return sorted list of available accYears as strings, if present."""
    if s is None:
        return []
    if not hasattr(s.index, "names") or "accYears" not in s.index.names:
        return []
    years = pd.Index(s.index.get_level_values("accYears")).unique().astype(str)
    # sort numerically if possible
    try:
        years_sorted = sorted(years, key=lambda x: int(str(x)))
    except Exception:
        years_sorted = sorted(map(str, years))
    return list(years_sorted)


def pick_hourly_year(scen_years: Sequence[str]) -> Optional[str]:
    if not scen_years:
        return None
    if HOURLY_YEAR_OVERRIDE is not None:
        if HOURLY_YEAR_OVERRIDE in scen_years:
            return HOURLY_YEAR_OVERRIDE
        # If override not present, fallback to max year but print warning
        print(f"[warn] HOURLY_YEAR_OVERRIDE={HOURLY_YEAR_OVERRIDE} not in {scen_years}; using max year.")
    # Default: max year
    try:
        return sorted(scen_years, key=lambda x: int(str(x)))[-1]
    except Exception:
        return scen_years[-1]


def scenario_gdx_path(case_name: str) -> Path:
    # Expected structure: ...\{case}\result\{case}.gdx
    return BASE_PROJECT / case_name / "result" / f"{case_name}.gdx"


def save_fig(fig: plt.Figure, outdir: Path, stem: str) -> None:
    ensure_dir(outdir)
    png = outdir / f"{stem}.png"
    svg = outdir / f"{stem}.svg"
    fig.savefig(png, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)


def pivot_year_tech(
    df: pd.DataFrame,
    year_col: str = "accYears",
    tech_col: str = "techs",
    value_col: str = "value",
    agg: str = "sum"
) -> pd.DataFrame:
    """Return DataFrame indexed by year with tech columns."""
    if df.empty:
        return pd.DataFrame()
    if year_col not in df.columns or tech_col not in df.columns or value_col not in df.columns:
        return pd.DataFrame()
    if agg == "sum":
        g = df.groupby([year_col, tech_col], as_index=False)[value_col].sum()
    else:
        g = df.groupby([year_col, tech_col], as_index=False)[value_col].mean()
    piv = g.pivot(index=year_col, columns=tech_col, values=value_col).fillna(0.0)
    # sort years numerically if possible
    try:
        piv = piv.reindex(sorted(piv.index, key=lambda x: int(str(x))))
    except Exception:
        piv = piv.sort_index()
    return piv


def sort_time_index(time_idx: pd.Index) -> pd.Index:
    """Try to sort timeModel robustly (supports t0001, 1, etc.)."""
    # Convert to string; extract digits if present
    s = pd.Index(time_idx).astype(str)

    def key(x: str) -> Tuple[int, str]:
        digits = "".join([c for c in x if c.isdigit()])
        if digits != "":
            try:
                return (int(digits), x)
            except Exception:
                return (10**18, x)
        return (10**18, x)

    order = sorted(list(s.unique()), key=key)
    return pd.Index(order)


def duration_curve(series: pd.Series) -> pd.Series:
    """Sort descending, index as rank (1..N)."""
    vals = np.asarray(series.values, dtype=float)
    vals = np.sort(vals)[::-1]
    return pd.Series(vals, index=np.arange(1, len(vals) + 1))


# -----------------------------
# Core Scenario container
# -----------------------------
@dataclass
class ScenarioData:
    name: str
    gdx_path: Path
    results: GDXEval

    converter_caps: Optional[pd.Series] = None
    commodity_balance_annual: Optional[pd.Series] = None
    commodity_balance_hourly: Optional[pd.Series] = None
    storage_flows: Optional[pd.Series] = None
    indicator_accounting: Optional[pd.Series] = None
    indicator_accounting_detailed: Optional[pd.Series] = None

    years_caps: List[str] = None
    years_annual: List[str] = None
    years_hourly: List[str] = None

    def load_symbols(self) -> None:
        self.converter_caps = safe_get_symbol(self.results, "converter_caps")
        self.commodity_balance_annual = safe_get_symbol(self.results, "commodity_balance_annual")
        self.commodity_balance_hourly = safe_get_symbol(self.results, "commodity_balance")
        self.storage_flows = safe_get_symbol(self.results, "storage_flows")
        self.indicator_accounting = safe_get_symbol(self.results, "indicator_accounting")
        self.indicator_accounting_detailed = safe_get_symbol(self.results, "indicator_accounting_detailed")

        self.years_caps = infer_years(self.converter_caps)
        self.years_annual = infer_years(self.commodity_balance_annual)
        self.years_hourly = infer_years(self.commodity_balance_hourly)

    def print_debug_overview(self) -> None:
        print(f"\n=== Scenario: {self.name} ===")
        print(f"gdx: {self.gdx_path}")
        print("Symbols/schema preview:")
        print_schema("converter_caps", self.converter_caps)
        print_schema("commodity_balance_annual", self.commodity_balance_annual)
        print_schema("commodity_balance (hourly)", self.commodity_balance_hourly)
        print_schema("storage_flows (hourly)", self.storage_flows)
        print_schema("indicator_accounting", self.indicator_accounting)
        print_schema("indicator_accounting_detailed", self.indicator_accounting_detailed)
        print(f"Detected years (caps):   {self.years_caps}")
        print(f"Detected years (annual): {self.years_annual}")
        print(f"Detected years (hourly): {self.years_hourly}")

    # ------------
    # Capacity
    # ------------
    def capacities_total_gw(self) -> pd.DataFrame:
        """
        Return total installed capacities (GW) by (accYears, techs) for GLOBAL_NODE and capType='total'.
        - prints debug if significant negatives exist.
        - returns pivot: index=year, columns=tech, values=GW
        """
        s = self.converter_caps
        if s is None:
            return pd.DataFrame()

        df = series_to_df(s).reset_index()

        required = {"accNodesModel", "accYears", "techs", "commodity", "capType", "value"}
        if not required.issubset(set(df.columns)):
            print(f"[warn] {self.name}: converter_caps missing expected columns {required - set(df.columns)}")
            return pd.DataFrame()

        # Filter to global node + total
        df = df[df["accNodesModel"].astype(str) == GLOBAL_NODE]
        df = df[df["capType"].astype(str) == "total"]

        # Units are expected GW already (per your note)
        # Handle negative capacities: ignore tiny negative, debug large negative
        neg = df[df["value"] < -CAP_NEG_TOL]
        if not neg.empty:
            for _, r in neg.iterrows():
                print(
                    f"[remix-framework debug] negative caps for {self.name}: "
                    f"{r['techs']} {r['commodity']} {r['accYears']} = {r['value']:.6g} GW"
                )

        # Clamp tiny negatives to 0, then drop negative values for plotting
        df.loc[df["value"].between(-CAP_NEG_TOL, 0.0), "value"] = 0.0

        # For plotting/reporting: keep only positive above a tiny threshold
        df = df[df["value"] > CAP_POS_TOL].copy()

        piv = pivot_year_tech(df, year_col="accYears", tech_col="techs", value_col="value", agg="sum")
        return piv

    # ------------
    # Annual Elec mix
    # ------------
    def annual_elec_net_by_tech(self) -> pd.DataFrame:
        """
        Annual net Elec balance by tech (GWh) at GLOBAL_NODE.
        Returns pivot: index=year, columns=tech, values=GWh (net, can be + or -)
        """
        s = self.commodity_balance_annual
        if s is None:
            return pd.DataFrame()

        df = series_to_df(s).reset_index()
        required = {"accNodesModel", "accYears", "techs", "commodity", "balanceType", "value"}
        if not required.issubset(set(df.columns)):
            print(f"[warn] {self.name}: commodity_balance_annual missing expected columns {required - set(df.columns)}")
            return pd.DataFrame()

        df = df[df["accNodesModel"].astype(str) == GLOBAL_NODE]
        df = df[df["commodity"].astype(str) == ELEC_COMMODITY]
        df = df[df["balanceType"].astype(str) == "net"]

        # Remove tiny noise
        df = df[df["value"].abs() > ENERGY_TOL].copy()

        piv = pivot_year_tech(df, year_col="accYears", tech_col="techs", value_col="value", agg="sum")
        return piv

    # ------------
    # Indicators
    # ------------
    def indicators_table(self, indicators: Sequence[str]) -> pd.DataFrame:
        """
        Attempt to extract indicator values by year for requested indicators.
        Works with indicator_accounting or indicator_accounting_detailed as available.
        Returns table indexed by year with indicator columns, values rounded later.
        """
        sources = []
        if self.indicator_accounting is not None:
            sources.append(("indicator_accounting", self.indicator_accounting))
        if self.indicator_accounting_detailed is not None:
            sources.append(("indicator_accounting_detailed", self.indicator_accounting_detailed))

        if not sources:
            return pd.DataFrame()

        best_df = pd.DataFrame()

        for src_name, s in sources:
            df = series_to_df(s).reset_index()

            # Heuristic: find which column holds indicator names
            # Common patterns in REMix: "indicators" or "indicator" or similar.
            possible_cols = [c for c in df.columns if c.lower() in {"indicator", "indicators", "accounting", "accIndicator".lower()}]
            if not possible_cols:
                # fallback: any col containing 'indicator'
                possible_cols = [c for c in df.columns if "indicator" in c.lower()]

            if not possible_cols:
                print(f"[warn] {self.name}: can't find indicator column in {src_name}. columns={list(df.columns)}")
                continue

            ind_col = possible_cols[0]

            # Heuristic for years
            if "accYears" not in df.columns:
                # maybe "year"?
                year_candidates = [c for c in df.columns if "year" in c.lower()]
                if not year_candidates:
                    print(f"[warn] {self.name}: can't find year column in {src_name}.")
                    continue
                year_col = year_candidates[0]
            else:
                year_col = "accYears"

            # Optional filter: global node if present
            if "accNodesModel" in df.columns:
                df = df[df["accNodesModel"].astype(str) == GLOBAL_NODE]

            df[ind_col] = df[ind_col].astype(str)
            df[year_col] = df[year_col].astype(str)

            df = df[df[ind_col].isin(list(indicators))].copy()
            if df.empty:
                continue

            # If multiple dims remain, sum over them
            group_cols = [year_col, ind_col]
            df2 = df.groupby(group_cols, as_index=False)["value"].sum()
            piv = df2.pivot(index=year_col, columns=ind_col, values="value").fillna(0.0)

            # Prefer the table with the most indicator columns found
            if piv.shape[1] > best_df.shape[1]:
                best_df = piv.copy()

        # Sort years
        if not best_df.empty:
            try:
                best_df = best_df.reindex(sorted(best_df.index, key=lambda x: int(str(x))))
            except Exception:
                best_df = best_df.sort_index()

        return best_df

    # ------------
    # Hourly Elec series
    # ------------
    def hourly_elec_by_tech(self, year: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Returns (generation_df, demand_df) for hourly Elec net by tech.
        generation_df: index=timeModel_sorted, columns=tech, values=positive GWh (>=0)
        demand_df:      index=timeModel_sorted, columns=tech, values=negative GWh (<=0)
        """
        s = self.commodity_balance_hourly
        if s is None:
            return pd.DataFrame(), pd.DataFrame()

        df = series_to_df(s).reset_index()

        # Required columns (timeModel plus key dims)
        if "timeModel" not in df.columns:
            # Try to find a time-like column
            time_candidates = [c for c in df.columns if "time" in c.lower()]
            if not time_candidates:
                print(f"[warn] {self.name}: hourly commodity_balance has no timeModel.")
                return pd.DataFrame(), pd.DataFrame()
            time_col = time_candidates[0]
        else:
            time_col = "timeModel"

        # Some REMix variants include balanceType here; you said possibly.
        # We prefer net if balanceType exists; otherwise we assume already net.
        if "balanceType" in df.columns:
            df = df[df["balanceType"].astype(str) == "net"]

        required = {"accNodesModel", "accYears", "techs", "commodity", "value"}
        missing = required - set(df.columns)
        if missing:
            print(f"[warn] {self.name}: hourly commodity_balance missing cols: {missing}")
            return pd.DataFrame(), pd.DataFrame()

        df = df[df["accNodesModel"].astype(str) == GLOBAL_NODE]
        df = df[df["accYears"].astype(str) == str(year)]
        df = df[df["commodity"].astype(str) == ELEC_COMMODITY]
        df = df[df["value"].abs() > ENERGY_TOL].copy()

        if df.empty:
            return pd.DataFrame(), pd.DataFrame()

        # Group to (time, tech)
        g = df.groupby([time_col, "techs"], as_index=False)["value"].sum()
        piv = g.pivot(index=time_col, columns="techs", values="value").fillna(0.0)

        # Sort time
        ordered_time = sort_time_index(piv.index)
        piv = piv.reindex(ordered_time)

        gen = piv.clip(lower=0.0)
        dem = piv.clip(upper=0.0)  # negative or 0
        return gen, dem

    def hourly_storage_by_tech(self, year: str) -> pd.DataFrame:
        """
        Try to read storage_flows (hourly). Returns pivot index=time, columns=techs (and maybe commodity),
        values=sum flows. If commodity exists, we keep separate columns tech|commodity.
        """
        s = self.storage_flows
        if s is None:
            return pd.DataFrame()

        df = series_to_df(s).reset_index()

        # Find time column
        if "timeModel" not in df.columns:
            time_candidates = [c for c in df.columns if "time" in c.lower()]
            if not time_candidates:
                return pd.DataFrame()
            time_col = time_candidates[0]
        else:
            time_col = "timeModel"

        if "accYears" not in df.columns:
            year_candidates = [c for c in df.columns if "year" in c.lower()]
            if not year_candidates:
                return pd.DataFrame()
            year_col = year_candidates[0]
        else:
            year_col = "accYears"

        # Optional filter global node if present
        if "accNodesModel" in df.columns:
            df = df[df["accNodesModel"].astype(str) == GLOBAL_NODE]

        df = df[df[year_col].astype(str) == str(year)].copy()
        if df.empty:
            return pd.DataFrame()

        # Build column key
        if "commodity" in df.columns:
            df["tech_comm"] = df["techs"].astype(str) + "|" + df["commodity"].astype(str)
            col = "tech_comm"
        else:
            col = "techs"

        g = df.groupby([time_col, col], as_index=False)["value"].sum()
        piv = g.pivot(index=time_col, columns=col, values="value").fillna(0.0)

        ordered_time = sort_time_index(piv.index)
        piv = piv.reindex(ordered_time)

        return piv


# -----------------------------
# Load scenarios
# -----------------------------
def load_all_scenarios() -> List[ScenarioData]:
    scenarios: List[ScenarioData] = []
    for name in SCENARIOS:
        gdx = scenario_gdx_path(name)
        if not gdx.exists():
            print(f"[ERROR] Missing GDX for {name}: {gdx}")
            continue
        try:
            res = GDXEval(str(gdx))
            sd = ScenarioData(name=name, gdx_path=gdx, results=res)
            sd.load_symbols()
            scenarios.append(sd)
        except Exception as e:
            print(f"[ERROR] Failed loading {name}: {e}")
            traceback.print_exc()
    return scenarios


# -----------------------------
# Tables export
# -----------------------------
def export_tables(
    outdir: Path,
    caps_tables: Dict[str, pd.DataFrame],
    gen_tables: Dict[str, pd.DataFrame],
    dem_tables: Dict[str, pd.DataFrame],
    ind_tables: Dict[str, pd.DataFrame],
) -> None:
    ensure_dir(outdir)

    # Save per-scenario CSV
    for scen, df in caps_tables.items():
        df.round(6).to_csv(outdir / f"table_caps_total_GW__{scen}.csv")
    for scen, df in gen_tables.items():
        df.round(6).to_csv(outdir / f"table_elec_generation_TWh__{scen}.csv")
    for scen, df in dem_tables.items():
        df.round(6).to_csv(outdir / f"table_elec_demand_TWh__{scen}.csv")
    for scen, df in ind_tables.items():
        df.round(6).to_csv(outdir / f"table_indicators__{scen}.csv")

    # Also write a single Excel workbook
    xlsx_path = outdir / "comparison_tables.xlsx"
    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        for scen, df in caps_tables.items():
            df.round(6).to_excel(writer, sheet_name=f"caps_{scen[:24]}")
        for scen, df in gen_tables.items():
            df.round(6).to_excel(writer, sheet_name=f"gen_{scen[:25]}")
        for scen, df in dem_tables.items():
            df.round(6).to_excel(writer, sheet_name=f"dem_{scen[:25]}")
        for scen, df in ind_tables.items():
            df.round(6).to_excel(writer, sheet_name=f"ind_{scen[:25]}")

    print(f"\n[ok] Wrote tables to: {outdir}")
    print(f"[ok] Excel workbook: {xlsx_path}")


# -----------------------------
# Plotting
# -----------------------------
def plot_multi_panel_stacked_bars(
    tables: Dict[str, pd.DataFrame],
    title: str,
    ylabel: str,
    outdir: Path,
    stem: str,
    convert_to_twh: bool = False,
    force_abs: bool = False,
    common_ylim: bool = True,
) -> None:
    """
    Create one subplot per scenario, stacked bars by year with tech columns.
    - convert_to_twh: if True divides values by 1000 (GWh->TWh)
    - force_abs: if True uses abs(values) (used for demand)
    """
    scen_names = list(tables.keys())
    n = len(scen_names)
    if n == 0:
        return

    # Determine global y max for consistent axis limits
    ymax = 0.0
    for scen, df in tables.items():
        if df.empty:
            continue
        vals = df.copy()
        if force_abs:
            vals = vals.abs()
        if convert_to_twh:
            vals = vals / 1000.0
        # stacked sum per year
        s = vals.sum(axis=1)
        if len(s) > 0:
            ymax = max(ymax, float(s.max()))
    if ymax <= 0:
        ymax = 1.0
    ymax *= 1.10

    cols = 2
    rows = int(math.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(13, 4.2 * rows), squeeze=False)
    fig.suptitle(title)

    for i, scen in enumerate(scen_names):
        r, c = divmod(i, cols)
        ax = axes[r][c]
        df = tables[scen]
        ax.set_title(scen)

        if df.empty:
            ax.text(0.5, 0.5, "No data", ha="center", va="center")
            ax.set_axis_off()
            continue

        vals = df.copy()
        if force_abs:
            vals = vals.abs()
        if convert_to_twh:
            vals = vals / 1000.0

        # Remove all-zero columns for cleaner legend
        vals = vals.loc[:, (vals.sum(axis=0) != 0)]
        # Plot
        vals.plot(kind="bar", stacked=True, ax=ax, width=0.85, legend=False)

        ax.set_xlabel("Year")
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", alpha=0.35)
        if common_ylim:
            ax.set_ylim(0, ymax)

    # Turn off unused axes
    for j in range(n, rows * cols):
        r, c = divmod(j, cols)
        axes[r][c].set_axis_off()

    # Create one global legend (from last non-empty plot)
    handles, labels = None, None
    for scen in scen_names:
        df = tables[scen]
        if not df.empty:
            tmp = df.copy()
            if force_abs:
                tmp = tmp.abs()
            if convert_to_twh:
                tmp = tmp / 1000.0
            tmp = tmp.loc[:, (tmp.sum(axis=0) != 0)]
            ax = plt.gca()
            # get from a dummy plot if needed
            break

    # Build legend from first non-empty axes
    for ax in fig.axes:
        h, l = ax.get_legend_handles_labels()
        if l:
            handles, labels = h, l
            break
    if labels:
        fig.legend(handles, labels, loc="upper center", ncol=5, bbox_to_anchor=(0.5, 0.98))

    save_fig(fig, outdir, stem)
    print(f"[ok] Saved plot: {stem}.png/.svg")


def plot_hourly_slices_and_duration(
    scenarios: List[ScenarioData],
    outdir: Path,
) -> None:
    """
    For each scenario:
      - Two-week winter slice: stacked generation + demand line (abs)
      - Two-week summer slice: stacked generation + demand line (abs)
      - Duration curves: total load (abs demand), total slack usage if found, and total generation
    Uses the *same hour windows* across scenarios in terms of sorted time index positions.
    """
    ensure_dir(outdir)

    # Precompute slice indices
    def slice_bounds(center: int) -> Tuple[int, int]:
        start = max(0, center - HOURS_PER_SLICE // 2)
        end = start + HOURS_PER_SLICE
        return start, end

    w0, w1 = slice_bounds(WINTER_CENTER_HOUR)
    s0, s1 = slice_bounds(SUMMER_CENTER_HOUR)

    for sd in scenarios:
        years = sd.years_hourly or []
        year = pick_hourly_year(years)
        if year is None:
            print(f"[warn] {sd.name}: no hourly years detected; skipping hourly plots.")
            continue

        gen, dem = sd.hourly_elec_by_tech(year)
        if gen.empty and dem.empty:
            print(f"[warn] {sd.name}: no hourly elec data for year={year}; skipping.")
            continue

        # Total series
        total_gen = gen.sum(axis=1)
        total_dem_abs = dem.sum(axis=1).abs()

        # Identify slack usage (if Slack techs exist in columns)
        slack_cols = [c for c in gen.columns if str(c).lower().startswith("slack")] + \
                     [c for c in dem.columns if str(c).lower().startswith("slack")]
        slack_series = None
        if slack_cols:
            # Slack might appear either positive or negative depending on modeling; take abs total
            slack_series = (gen[slack_cols].sum(axis=1) + dem[slack_cols].sum(axis=1)).abs()

        # ---- Two-week slice plots
        def plot_slice(start: int, end: int, season: str) -> None:
            fig, ax1 = plt.subplots(figsize=(13, 4.8))
            ax1.set_title(f"{sd.name} | Elec hourly operation | {season} slice | year={year} | hours {start}:{end}")
            # stacked generation
            gg = gen.iloc[start:end].copy()
            gg = gg.loc[:, (gg.sum(axis=0) != 0)]
            gg.plot.area(stacked=True, ax=ax1, linewidth=0.0)
            ax1.set_ylabel("Generation (GWh per timestep)")
            ax1.grid(alpha=0.25)

            ax2 = ax1.twinx()
            dd = dem.iloc[start:end].sum(axis=1).abs()
            ax2.plot(dd.index.astype(str), dd.values, linewidth=1.8)
            ax2.set_ylabel("Demand abs (GWh per timestep)")
            ax2.set_ylim(ax1.get_ylim())

            # reduce x tick clutter
            xt = np.linspace(0, max(1, len(dd) - 1), 8).astype(int)
            ax1.set_xticks(dd.index[xt])
            ax1.set_xticklabels([str(x) for x in dd.index[xt]], rotation=0)

            save_fig(fig, outdir, f"hourly_{sd.name}__{season}_2weeks__year_{year}")

        plot_slice(w0, w1, "winter")
        plot_slice(s0, s1, "summer")

        # ---- Duration curves
        fig, ax = plt.subplots(figsize=(10.5, 5))
        ax.set_title(f"{sd.name} | Duration curves | year={year}")

        dc_load = duration_curve(total_dem_abs)
        dc_gen = duration_curve(total_gen)
        ax.plot(dc_load.index, dc_load.values, label="Load (abs demand)")
        ax.plot(dc_gen.index, dc_gen.values, label="Total generation")

        if slack_series is not None:
            dc_slack = duration_curve(slack_series)
            ax.plot(dc_slack.index, dc_slack.values, label="Slack usage (abs)")

        ax.set_xlabel("Hour rank (descending)")
        ax.set_ylabel("GWh per timestep")
        ax.grid(alpha=0.3)
        ax.legend()

        save_fig(fig, outdir, f"duration_{sd.name}__year_{year}")
        print(f"[ok] Saved hourly slice + duration plots for {sd.name} (year={year}).")


def plot_storage_overview(
    scenarios: List[ScenarioData],
    outdir: Path,
) -> None:
    """
    If storage_flows exists, plots total storage flow by tech|commodity for the chosen hourly year.
    This is a fallback when SOC isn't available.
    """
    ensure_dir(outdir)

    for sd in scenarios:
        years = sd.years_hourly or []
        year = pick_hourly_year(years)
        if year is None:
            continue
        stor = sd.hourly_storage_by_tech(year)
        if stor.empty:
            print(f"[info] {sd.name}: no storage_flows data for year={year}.")
            continue

        # Plot a 2-week winter slice of storage flows for readability
        start = max(0, WINTER_CENTER_HOUR - HOURS_PER_SLICE // 2)
        end = start + HOURS_PER_SLICE
        ss = stor.iloc[start:end].copy()
        ss = ss.loc[:, (ss.sum(axis=0).abs() != 0)]

        # If too many columns, keep top 12 by absolute sum
        if ss.shape[1] > 12:
            top = ss.abs().sum(axis=0).sort_values(ascending=False).head(12).index
            ss = ss[top]

        fig, ax = plt.subplots(figsize=(13, 4.8))
        ax.set_title(f"{sd.name} | storage_flows (fallback) | winter 2-week slice | year={year}")
        ss.plot(ax=ax, linewidth=1.3)
        ax.set_xlabel("timeModel")
        ax.set_ylabel("Storage flow (model units per timestep)")
        ax.grid(alpha=0.25)
        ax.legend(loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.15))
        save_fig(fig, outdir, f"storage_flows_{sd.name}__winter_2weeks__year_{year}")
        print(f"[ok] Saved storage plot for {sd.name} (year={year}).")


def try_maps_optional(
    scenarios: List[ScenarioData],
    outdir: Path,
) -> None:
    """
    Optional: maps/network plots if shapefiles are present and geopandas + remix plot utils import works.
    Uses path_geo = C:/Local/REMix/remix_nz/input/shapefiles and geofile="11regionsNZ.geojson" (as you suggested).
    """
    if not HAS_GEO:
        print("[info] Maps skipped: geopandas/remix plot utils not available in this environment.")
        return

    path_base = Path("C:/Local/REMix")
    path_geo = path_base / "remix_nz" / "input" / "shapefiles"
    geojson = path_geo / "11regionsNZ.geojson"
    shp = path_geo / "11regionsNZ.shp"
    shp_attrcol = "id"

    if not geojson.exists() and not shp.exists():
        print(f"[info] Maps skipped: no geojson/shp found in {path_geo}")
        return

    ensure_dir(outdir / "maps")

    # Example: map annual PV generation by region in the latest year (per scenario),
    # only if commodity_balance_annual includes regional nodes.
    for sd in scenarios:
        s = sd.commodity_balance_annual
        if s is None:
            continue
        years = sd.years_annual or []
        if not years:
            continue
        year = years[-1]

        df = series_to_df(s).reset_index()
        required = {"accNodesModel", "accYears", "techs", "commodity", "balanceType", "value"}
        if not required.issubset(df.columns):
            continue

        # skip if only global node exists
        if df["accNodesModel"].astype(str).nunique() <= 1:
            continue

        df = df[(df["accYears"].astype(str) == str(year)) &
                (df["commodity"].astype(str) == ELEC_COMMODITY) &
                (df["balanceType"].astype(str) == "net")].copy()
        if df.empty:
            continue

        # pick PV-like techs by name contains 'pv'
        df["techs_l"] = df["techs"].astype(str).str.lower()
        df_pv = df[df["techs_l"].str.contains("pv")].copy()
        if df_pv.empty:
            continue

        # Aggregate by region node
        map_data = df_pv.groupby("accNodesModel", as_index=True)["value"].sum() / 1000.0  # TWh
        map_data = map_data.to_frame(name="value")

        try:
            fig = plt.figure(figsize=(10.5, 5.5))
            plot_choropleth(
                data=map_data,
                shp=str(shp) if shp.exists() else str(geojson),
                shp_attr=shp_attrcol,
                title=f"{sd.name} | PV generation | {year}",
                clabel="TWh",
            )
            save_fig(fig, outdir / "maps", f"map_pv_gen_{sd.name}__{year}")
            print(f"[ok] Saved map for {sd.name} (year={year}).")
        except Exception as e:
            print(f"[warn] Map failed for {sd.name}: {e}")


# -----------------------------
# Main
# -----------------------------
def main() -> None:
    ensure_dir(OUT_DIR)
    print(f"[info] Output folder: {OUT_DIR}")

    scenarios = load_all_scenarios()
    if not scenarios:
        print("[ERROR] No scenarios loaded. Check your BASE_PROJECT path and scenario folder names.")
        sys.exit(1)

    # 1) Print schema/debug for each scenario so you can show me structure
    for sd in scenarios:
        sd.print_debug_overview()

    # 2) Build key tables
    caps_tables: Dict[str, pd.DataFrame] = {}
    elec_net_tables: Dict[str, pd.DataFrame] = {}
    elec_gen_tables: Dict[str, pd.DataFrame] = {}
    elec_dem_tables: Dict[str, pd.DataFrame] = {}
    ind_tables: Dict[str, pd.DataFrame] = {}

    indicators_needed = ["CO2_emission", "FuelCost", "SlackCost", "SystemCost"]

    for sd in scenarios:
        caps = sd.capacities_total_gw()
        caps_tables[sd.name] = caps

        elec_net = sd.annual_elec_net_by_tech()
        elec_net_tables[sd.name] = elec_net

        # generation = positive net elec
        gen = elec_net.clip(lower=0.0)
        # demand = negative net elec (keep negative in table, but also store abs version as "demand table" below)
        dem = elec_net.clip(upper=0.0)

        # For annual stacks in your paper, you wanted:
        # - generation mix = positive net elec (stack)
        # - demand/consumption = negative net elec plotted as abs (stack)
        elec_gen_tables[sd.name] = gen
        elec_dem_tables[sd.name] = dem  # still negative here; plotting will abs()

        inds = sd.indicators_table(indicators_needed)
        ind_tables[sd.name] = inds

        # Console prints (rounded to 2 decimals)
        print(f"\n--- TABLES (rounded) | {sd.name} ---")
        if not caps.empty:
            print("Installed capacities (GW) capType=total, global:")
            print(caps.round(2))
        else:
            print("Installed capacities: (no data)")

        if not elec_net.empty:
            print("\nAnnual Elec net balance by tech (GWh), global:")
            print(elec_net.round(2))
            print("\nAnnual Elec generation mix (positive net, GWh), global:")
            print(gen.round(2))
            print("\nAnnual Elec demand/consumption (negative net, GWh), global (kept negative):")
            print(dem.round(2))
        else:
            print("Annual Elec balance: (no data)")

        if not inds.empty:
            print("\nIndicators (raw units as in GDX), global:")
            print(inds.round(2))
        else:
            print("\nIndicators: (no data)")

    # 3) Export tables (CSV + Excel)
    export_tables(OUT_DIR, caps_tables, elec_gen_tables, elec_dem_tables, ind_tables)

    # 4) Plots: multi-panel per scenario
    plot_multi_panel_stacked_bars(
        tables=caps_tables,
        title="Installed capacities (capType=total, global) by year",
        ylabel="GW",
        outdir=OUT_DIR,
        stem="caps_total_GW__stacked_by_year__panels",
        convert_to_twh=False,
        force_abs=False,
        common_ylim=True,
    )

    plot_multi_panel_stacked_bars(
        tables=elec_gen_tables,
        title="Electricity generation mix (positive net Elec, global) by year",
        ylabel="TWh",
        outdir=OUT_DIR,
        stem="elec_generation_TWh__stacked_by_year__panels",
        convert_to_twh=True,   # your annual balance is in GWh -> show TWh
        force_abs=False,
        common_ylim=True,
    )

    plot_multi_panel_stacked_bars(
        tables=elec_dem_tables,
        title="Electricity demand/consumption (abs(negative net Elec), global) by year",
        ylabel="TWh",
        outdir=OUT_DIR,
        stem="elec_demand_TWh__stacked_by_year__panels",
        convert_to_twh=True,
        force_abs=True,  # plot abs value
        common_ylim=True,
    )

    # 5) Hourly comparisons (2-week slices + duration curves)
    plot_hourly_slices_and_duration(scenarios, OUT_DIR / "hourly")

    # 6) Storage overview (fallback: storage_flows)
    plot_storage_overview(scenarios, OUT_DIR / "storage")

    # 7) Optional maps (only if deps + files exist)
    try_maps_optional(scenarios, OUT_DIR)

    print("\n[done] Comparison outputs written to:")
    print(f"  {OUT_DIR}\n")


if __name__ == "__main__":
    main()

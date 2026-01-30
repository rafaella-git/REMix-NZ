from __future__ import annotations

import warnings
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib import colormaps, cm
from remix.framework.tools.gdx import GDXEval


# ---------------------------------------------------------------------------
# BASIC CONFIG
# ---------------------------------------------------------------------------

warnings.filterwarnings(
    "ignore",
    message=r".*GAMS version.*differs from the API version.*",
    category=UserWarning,
)

BASE_DIR = Path(r"C:\Local\REMix\remix_nz\project\GP-NT-ELEC-BIO-H2")

CASE_DIRS = [
    "nz_case_GP_2020-2025-2030-2035-2040-2045-2050",
    "nz_case_NT_2020-2025-2030-2035-2040-2045-2050",
    "nz_case_ELEC+_2020-2025-2030-2035-2040-2045-2050",
    "nz_case_BIO+_2020-2025-2030-2035-2040-2045-2050",
    "nz_case_H2+_2020-2025-2030-2035-2040-2045-2050",
]

YEAR_INTS = [2020, 2025, 2030, 2035, 2040, 2045, 2050]
SCEN_ORDER = ["GP", "NT", "ELEC+", "BIO+", "H2+"]

DX_IN = 1.0
YEAR_GAP = 1.5
YEAR_PAD = 55
BAR_WIDTH = 0.85

GEN_TECH_MIN_TWH = 0.01
CAP_TECH_MIN_GW = 0.01  # 10 MW

EXPORT_FIGURES = False  # toggled by user only; if False, only plt.show() is used

FUEL_CONV_TECHS = ["Methanizer", "FTropschSyn", "Electrolyser", "DAC"]
FUEL_CONV_STACK_ORDER = ["Electrolyser", "Methanizer", "FTropschSyn", "DAC"]

COST_COMPONENTS = ["FuelCost", "Invest", "OMFix", "SlackCost", "SpillPenalty"]

technology_colors = {
    "Hydropower": "#298c81",
    "Geothermal": "#ba91b1",
    "Solar PV": "#f9d002",
    "Offshore wind": "#6895dd",
    "Onshore wind": "#4F7EC9",
    "Battery": "#708090",
    "H2 Storage": "#bf13a0",
    "Fuel Cell (H2)": "#c251ae",
    "Electrolyser": "#AF9968",
    "CCGT": "#ee8340",
    "OCGT": "#FAA460",
    "GT": "#ee8340",
    "Gas turbines": "#ee8340",
    "Diesel": "#B5A642",
    "Coal": "#505050",
    "Biomass": "#008000",
}
DEFAULT_COLOR = "#9E9E9E"


# ---------------------------------------------------------------------------
# UTILS / LOADING
# ---------------------------------------------------------------------------

def scenario_name(case_dir: str) -> str:
    if "_GP_" in case_dir:
        return "GP"
    if "_NT_" in case_dir:
        return "NT"
    if "_ELEC+_" in case_dir:
        return "ELEC+"
    if "_BIO+_" in case_dir:
        return "BIO+"
    if "_H2+_" in case_dir:
        return "H2+"
    return case_dir


def load_results(case_dir: str) -> GDXEval | None:
    gdx_path = BASE_DIR / case_dir / "result" / f"{case_dir}.gdx"
    if not gdx_path.is_file():
        print("Missing:", gdx_path)
        return None
    return GDXEval(str(gdx_path))


def read_symbol(results: GDXEval, name: str) -> pd.DataFrame | None:
    try:
        df = results[name]
    except KeyError:
        return None
    if isinstance(df, pd.Series):
        df = df.to_frame("value")
    if "value" not in df.columns:
        df = df.copy()
        df.columns = ["value"]
    return df


def ensure_year_col(df: pd.DataFrame, year_col_candidates=("accYears", "years")) -> pd.DataFrame:
    df = df.copy()
    year_col = None
    for c in year_col_candidates:
        if c in df.columns:
            year_col = c
            break
    if year_col is None:
        raise KeyError(
            f"Could not find a year column. Tried {year_col_candidates}. "
            f"Columns: {list(df.columns)}"
        )
    df["year"] = pd.to_numeric(df[year_col], errors="coerce")
    df = df[df["year"].isin(YEAR_INTS)].copy()
    df["year"] = df["year"].astype(int)
    return df


def load_needed_data():
    cba_list, caps_list, ind_list, ind_det_list = [], [], [], []

    for case_dir in CASE_DIRS:
        scen = scenario_name(case_dir)
        res = load_results(case_dir)
        if res is None:
            continue

        cba = read_symbol(res, "commodity_balance_annual")
        if cba is not None:
            df = ensure_year_col(cba.reset_index(), year_col_candidates=("accYears",))
            df["scenario"] = scen
            cba_list.append(df)

        caps = read_symbol(res, "converter_caps")
        if caps is not None:
            df = ensure_year_col(caps.reset_index(), year_col_candidates=("accYears",))
            df["scenario"] = scen
            caps_list.append(df)

        ind = read_symbol(res, "indicator_accounting")
        if ind is not None:
            df = ensure_year_col(ind.reset_index(), year_col_candidates=("accYears",))
            df["scenario"] = scen
            ind_list.append(df)

        ind_det = read_symbol(res, "indicator_accounting_detailed")
        if ind_det is not None:
            df = ensure_year_col(ind_det.reset_index(), year_col_candidates=("years", "accYears"))
            df["scenario"] = scen
            ind_det_list.append(df)

    if not cba_list or not caps_list:
        raise RuntimeError("No data loaded. Check BASE_DIR / CASE_DIRS and GDX files.")

    cba_long = pd.concat(cba_list, ignore_index=True)
    caps_long = pd.concat(caps_list, ignore_index=True)
    ind_long = pd.concat(ind_list, ignore_index=True) if ind_list else pd.DataFrame()
    ind_det_long = pd.concat(ind_det_list, ignore_index=True) if ind_det_list else pd.DataFrame()

    cba_long["scenario"] = pd.Categorical(cba_long["scenario"], categories=SCEN_ORDER, ordered=True)
    caps_long["scenario"] = pd.Categorical(caps_long["scenario"], categories=SCEN_ORDER, ordered=True)
    if not ind_long.empty:
        ind_long["scenario"] = pd.Categorical(ind_long["scenario"], categories=SCEN_ORDER, ordered=True)
    if not ind_det_long.empty:
        ind_det_long["scenario"] = pd.Categorical(
            ind_det_long["scenario"], categories=SCEN_ORDER, ordered=True
        )

    return cba_long, caps_long, ind_long, ind_det_long


# ---------------------------------------------------------------------------
# X-AXIS MULTIINDEX HELPERS
# ---------------------------------------------------------------------------

def build_grouped_x(
    years: Iterable[int],
    scen_order: list[str],
    dx_in: float = 1.0,
    year_gap: float = 2.0,
):
    years = list(years)
    x, year_centers = [], []
    pos = 0.0
    for _y in years:
        xs = [pos + i * dx_in for i in range(len(scen_order))]
        x.extend(xs)
        year_centers.append(float(np.mean(xs)))
        pos = xs[-1] + dx_in + year_gap
    return np.array(x), np.array(year_centers)


def apply_two_level_x(
    ax: plt.Axes,
    x: np.ndarray,
    years: list[int],
    year_centers: np.ndarray,
    scen_order: list[str],
    year_pad: int = 74,
    scen_pad: int = 2,
    show_scen: bool = True,
):
    if show_scen:
        scen_labels = [s for _y in years for s in scen_order]
        ax.set_xticks(x)
        ax.set_xticklabels(scen_labels, rotation=90, va="top")
        ax.tick_params(axis="x", pad=scen_pad)
    else:
        ax.set_xticks([])
        ax.set_xticklabels([])

    sec = ax.secondary_xaxis("bottom")
    sec.set_xlim(ax.get_xlim())
    sec.set_xticks(year_centers)
    sec.set_xticklabels([str(y) for y in years], rotation=0)
    sec.spines["bottom"].set_visible(False)
    sec.tick_params(axis="x", pad=year_pad, length=0)


# ---------------------------------------------------------------------------
# MAPPINGS / COLORS
# ---------------------------------------------------------------------------

def tech_group_pretty(tech: str) -> str:
    if tech in {"pv_central_fixed", "pv_decentral"}:
        return "Solar PV"
    if tech.startswith("wind_onshore_"):
        return "Onshore wind"
    if tech.startswith("wind_offshore_"):
        return "Offshore wind"
    if tech == "Hydro":
        return "Hydropower"
    if tech in {"GT", "OCGT", "CCGT"}:
        return "Gas turbines"
    if tech == "Thermal_Bio":
        return "Biomass"
    if tech == "Thermal_Coal":
        return "Coal"
    if tech == "Thermal_Diesel":
        return "Diesel"
    if tech == "H2_FC":
        return "Fuel Cell (H2)"
    return tech.replace("_", " ")


def color_map_for(categories: list[str]) -> dict[str, str]:
    return {c: technology_colors.get(c, DEFAULT_COLOR) for c in categories}


def evenly_spaced_colors_from_cmap(cmap_name: str, n: int, start: float = 0.0, stop: float = 1.0):
    cmap = cm.get_cmap(cmap_name)
    return [cmap(x) for x in np.linspace(start, stop, n)]


# ---------------------------------------------------------------------------
# GENERIC STACKED GROUPED BAR
# ---------------------------------------------------------------------------

def stacked_grouped_bars(
    ax: plt.Axes,
    df_long: pd.DataFrame,
    value_col: str,
    category_col: str,
    y_label: str,
    years: list[int],
    scen_order: list[str],
    cmap: dict[str, str],
    cat_order: list[str],
    dx_in: float = 1.0,
    year_gap: float = 2.0,
    year_pad: int = 74,
    show_scen_labels: bool = True,
):
    x, year_centers = build_grouped_x(years, scen_order, dx_in=dx_in, year_gap=year_gap)

    piv = (
        df_long.pivot_table(
            index=["year", "scenario"],
            columns=category_col,
            values=value_col,
            aggfunc="sum",
            observed=False,
        )
        .reindex(pd.MultiIndex.from_product([years, scen_order], names=["year", "scenario"]))
        .fillna(0.0)
    )
    cat_order = [c for c in cat_order if c in piv.columns]
    piv = piv[cat_order]

    bottom = np.zeros(len(x))
    for c in piv.columns:
        vals = piv[c].values
        ax.bar(x, vals, bottom=bottom, width=BAR_WIDTH, color=cmap.get(c, DEFAULT_COLOR), edgecolor="none")
        bottom += vals

    ax.grid(axis="y", alpha=0.25)
    ax.set_ylabel(y_label)

    apply_two_level_x(
        ax=ax,
        x=x,
        years=years,
        year_centers=year_centers,
        scen_order=scen_order,
        year_pad=year_pad,
        show_scen=show_scen_labels,
    )

    ymax = float((piv.sum(axis=1)).max())
    ax.set_ylim(0.0, ymax * 1.05 if ymax > 0 else 1.0)

    return piv.columns.tolist()


def shared_legend(
    fig: plt.Figure,
    labels: list[str],
    cmap: dict[str, str],
    title: str,
    x: float,
    y: float,
):
    handles = [plt.Line2D([0], [0], color=cmap.get(l, DEFAULT_COLOR), lw=10) for l in labels]
    leg = fig.legend(
        handles,
        labels,
        title=title,
        loc="center left",
        bbox_to_anchor=(x, y),
        frameon=True,
        borderpad=1.2,
        labelspacing=0.8,
        handlelength=1.8,
        handletextpad=0.8,
    )
    leg.get_frame().set_linewidth(1.2)
    return leg


# ---------------------------------------------------------------------------
# DATA BUILDERS FOR RQ1
# ---------------------------------------------------------------------------

def build_capacity_elec(caps_long: pd.DataFrame) -> pd.DataFrame:
    df = caps_long[
        (caps_long["accNodesModel"] == "global")
        & (caps_long["capType"] == "total")
        & (caps_long["commodity"] == "Elec")
        & (caps_long["year"].isin(YEAR_INTS))
    ].copy()
    df = df[df["value"].abs() > CAP_TECH_MIN_GW].copy()
    df["tech"] = df["techs"].astype(str).map(tech_group_pretty)
    out = df.groupby(["year", "scenario", "tech"], as_index=False, observed=False)["value"].sum()
    return out


def build_generation_elec(cba_long: pd.DataFrame) -> pd.DataFrame:
    df = cba_long[
        (cba_long["accNodesModel"] == "global")
        & (cba_long["commodity"] == "Elec")
        & (cba_long["balanceType"] == "net")
        & (cba_long["year"].isin(YEAR_INTS))
        & (cba_long["value"] > 0)
    ].copy()
    df["tech"] = df["techs"].astype(str).map(tech_group_pretty)
    df["value_twh"] = df["value"] / 1000.0
    df = df[df["value_twh"].abs() > GEN_TECH_MIN_TWH].copy()
    out = df.groupby(["year", "scenario", "tech"], as_index=False, observed=False)["value_twh"].sum()
    return out


def build_fuelconv_caps(caps_long: pd.DataFrame) -> pd.DataFrame:
    df = caps_long[
        (caps_long["accNodesModel"] == "global")
        & (caps_long["capType"] == "total")
        & (caps_long["year"].isin(YEAR_INTS))
        & (caps_long["techs"].astype(str).isin(FUEL_CONV_TECHS))
    ].copy()
    df = df[df["value"].abs() > 0.01].copy()
    df["tech"] = df["techs"].astype(str)
    return df[["year", "scenario", "tech", "value"]]


def build_fuel_supply_groups(cba_long: pd.DataFrame) -> pd.DataFrame:
    df = cba_long[
        (cba_long["year"].isin(YEAR_INTS))
        & (cba_long["balanceType"] == "net")
        & (cba_long["value"] > 0)
    ].copy()

    # Slack warnings
    slack = df[df["techs"].astype(str).str.startswith("SlackFuel_")].copy()
    if not slack.empty:
        check = (
            slack.groupby(["scenario", "year"], observed=False)["value"]
            .sum()
            .abs()
            .reset_index()
        )
        warn = check[check["value"] > 1.0]
        for _, row in warn.iterrows():
            print(
                f"Warning: SlackFuel_* absolute value > 1 GWh for scenario={row['scenario']} year={row['year']}"
            )
    df = df[~df["techs"].astype(str).str.startswith("SlackFuel_")].copy()

    def fuel_group(row) -> str | None:
        t = str(row["techs"])
        comm = str(row["commodity"])

        if comm == "e-CH4":
            return "e-CH4"
        if comm == "REfuel":
            return "e-FTL fuels"
        if comm == "H2":
            return "e-hydrogen"

        if t in {"FuelImport_Bio_LF", "FuelImport_Biofuel"}:
            return "Biofuels liquid"
        if t == "FuelImport_Bio_gas":
            return "Biogas"
        if t in {"FuelImport_CH4", "FuelImport_Fossil_CH4"}:
            return "Fossil methane"
        if t == "FuelImport_Coal":
            return "Fossil coal"
        if t in {"FuelImport_Diesel", "FuelImport_Fossil_LF"}:
            return "Fossil oil"
        return None

    df["fuel_group"] = df.apply(fuel_group, axis=1)
    df = df[df["fuel_group"].notna()].copy()
    df["value_twh"] = df["value"] / 1000.0
    df = df[df["value_twh"].abs() > 0.01].copy()

    out = df.groupby(["year", "scenario", "fuel_group"], as_index=False, observed=False)["value_twh"].sum()
    return out


# ---------------------------------------------------------------------------
# DATA BUILDERS FOR RQ3
# ---------------------------------------------------------------------------

def drop_doublecounted_fuel_import_fuelcost(df_cost: pd.DataFrame) -> pd.DataFrame:
    df = df_cost.copy()
    if "indicator" not in df.columns or "techs" not in df.columns:
        raise KeyError(f"Need columns ['indicator','techs']. Got: {list(df.columns)}")
    is_fuelcost = df["indicator"].astype(str) == "FuelCost"
    is_import = df["techs"].astype(str).str.startswith("FuelImport_", na=False)
    return df[~(is_fuelcost & is_import)].copy()


def build_cost_long_from_detailed(ind_det_long: pd.DataFrame) -> pd.DataFrame:
    if ind_det_long.empty:
        return pd.DataFrame(columns=["year", "scenario", "component", "value_bEUR"])

    df = ind_det_long.copy()
    df = df[(df["year"].isin(YEAR_INTS)) & (df["indicator"].isin(COST_COMPONENTS))].copy()
    df = drop_doublecounted_fuel_import_fuelcost(df)

    label_map = {
        "FuelCost": "FuelCost",
        "Invest": "Invest",
        "OMFix": "OMFix",
        "SlackCost": "SlackCost",
        "SpillPenalty": "SpillPenalty",
    }
    df["component"] = df["indicator"].map(lambda x: label_map.get(str(x), str(x)))

    agg = df.groupby(["year", "scenario", "component"], as_index=False, observed=False)["value"].sum()

    agg["value_bEUR"] = agg["value"] / 1000.0
    agg = agg[agg["value_bEUR"].abs() > 0.01].copy()

    return agg[["year", "scenario", "component", "value_bEUR"]]


def build_lcoe(ind_det_long: pd.DataFrame, cba_long: pd.DataFrame) -> pd.DataFrame:
    if ind_det_long.empty:
        return pd.DataFrame(columns=["year", "scenario", "LCOE_EUR_MWh"])

    df = ind_det_long.copy()
    df = df[df["year"].isin(YEAR_INTS)].copy()
    df = df[df["indicator"].isin(["FuelCost", "Invest", "OMFix"])].copy()

    include_techs = {
        "GT",
        "Thermal_Coal",
        "CCGT",
        "OCGT",
        "Thermal_Diesel",
        "Thermal_Bio",
        "Hydro",
        "H2_CCGT",
        "H2_FC",
        "pv_central_fixed",
        "pv_decentral",
        "wind_offshore_1",
        "wind_offshore_2",
        "wind_offshore_3",
        "wind_offshore_4",
        "wind_onshore_1",
        "wind_onshore_2",
        "wind_onshore_3",
        "wind_onshore_4",
        "FuelImport_CH4",
        "FuelImport_Coal",
        "FuelImport_Biofuel",
        "FuelImport_Diesel",
        "HV",
        "Battery",
    }
    exclude_techs = {
        "HeatDemand_Fossil_CH4",
        "HeatDemand_Fossil_Coal",
        "HeatDemand_Fossil_LF",
        "TranspDemand_Fossil_LF",
        "HeatDemand_e-CH4",
        "TranspDemand_REfuel",
        "TranspDemand_e-CH4",
        "FuelImport_Fossil_CH4",
        "FuelImport_Fossil_Coal",
        "FuelImport_Fossil_LF",
        "FuelImport_Bio_LF",
        "FuelImport_Bio_gas",
        "FuelImport_Biomass",
        "Elec_slack",
        "SlackFuel_H2",
        "SlackFuel_REfuel",
        "H2_storage",
        "Electrolyser",
        "FTropschSyn",
        "Methanizer",
    }

    tech = df["techs"].astype(str)
    df = df[tech.isin(include_techs) & ~tech.isin(exclude_techs)].copy()

    df = drop_doublecounted_fuel_import_fuelcost(df)

    df["cost_eur"] = df["value"] * 1e6
    cost_agg = df.groupby(["year", "scenario"], as_index=False, observed=False)["cost_eur"].sum()

    gen = cba_long[
        (cba_long["accNodesModel"] == "global")
        & (cba_long["commodity"] == "Elec")
        & (cba_long["balanceType"] == "net")
        & (cba_long["year"].isin(YEAR_INTS))
        & (cba_long["value"] > 0)
    ].copy()
    gen_gwh = gen.groupby(["year", "scenario"], as_index=False, observed=False)["value"].sum()
    gen_gwh = gen_gwh.rename(columns={"value": "gen_gwh"})
    gen_gwh["gen_mwh"] = gen_gwh["gen_gwh"] * 1000.0

    merged = cost_agg.merge(gen_gwh, on=["year", "scenario"], how="left")
    merged["LCOE_EUR_MWh"] = (merged["cost_eur"] / merged["gen_mwh"]).replace(
        [np.inf, -np.inf], np.nan
    )
    merged = merged.dropna(subset=["LCOE_EUR_MWh"])
    merged = merged[merged["LCOE_EUR_MWh"].abs() > 0.01].copy()

    return merged[["year", "scenario", "LCOE_EUR_MWh"]]


def build_total_co2(ind_long: pd.DataFrame) -> pd.DataFrame:
    df = ind_long[
        (ind_long["accNodesModel"] == "global")
        & (ind_long["indicator"] == "CO2_emission")
        & (ind_long["year"].isin(YEAR_INTS))
    ].copy()
    df["value_mt"] = df["value"] / 1000.0
    df = df[df["value_mt"].abs() > 0.01].copy()
    out = df.groupby(["year", "scenario"], as_index=False, observed=False)["value_mt"].sum()
    return out


def build_co2_contributors_2050(ind_det_long: pd.DataFrame) -> pd.DataFrame:
    if ind_det_long.empty:
        return pd.DataFrame(columns=["scenario", "category", "value_mt"])

    df = ind_det_long.copy()
    df = df[(df["indicator"] == "CO2_emission") & (df["year"] == 2050)].copy()
    df["value_mt"] = df["value"] / 1000.0

    tech = df["techs"].astype(str)

    def cat_map(t: str) -> str | None:
        if t in {
            "HeatDemand_Fossil_LF",
            "HeatDemand_Fossil_CH4",
            "HeatDemand_REfuel",
            "HeatDemand_e-CH4",
        }:
            return "Heat"
        if t in {
            "TranspDemand_Fossil_LF",
            "TranspDemand_REfuel",
            "TranspDemand_e-CH4",
        }:
            return "Transport"
        if t in {"CCGT", "GT", "OCGT"}:
            return "Electricity generation (gas)"
        if t == "Thermal_Coal":
            return "Electricity generation (coal)"
        if t == "Thermal_Diesel":
            return "Electricity generation (diesel)"
        if t == "DAC":
            return "Direct air capture"
        return t

    df["category"] = tech.map(cat_map)
    df = df[df["category"].notna()].copy()
    df = df[df["value_mt"].abs() > 0.01].copy()

    out = df.groupby(["scenario", "category"], as_index=False, observed=False)["value_mt"].sum()
    return out


# ---------------------------------------------------------------------------
# PLOTTING FUNCTIONS
# ---------------------------------------------------------------------------

def plot_rq1_a_capacity_lines(caps_long: pd.DataFrame, outdir: Path | None = None):
    sns.set_theme(style="whitegrid")

    cap = build_capacity_elec(caps_long)
    if cap.empty:
        print("RQ1.A: no capacity data.")
        return

    years = sorted(cap["year"].unique())
    techs = sorted(cap["tech"].unique())
    cmap = color_map_for(techs)

    fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    for tech in techs:
        df_t = cap[cap["tech"] == tech]
        series = df_t.groupby("year", as_index=False, observed=False)["value"].sum()
        ax.plot(
            series["year"],
            series["value"],
            label=tech,
            color=cmap.get(tech, DEFAULT_COLOR),
            marker="o",
        )

    ax.set_xlabel("Year")
    ax.set_ylabel("Installed electricity generation capacity (GW)")
    ax.grid(True, axis="y", alpha=0.25)

    handles, labels = ax.get_legend_handles_labels()
    leg = ax.legend(
        handles,
        labels,
        title="Technology",
        frameon=True,
    )
    leg.get_frame().set_linewidth(1.0)

    if EXPORT_FIGURES and outdir is not None:
        outdir.mkdir(parents=True, exist_ok=True)
        fig.savefig(outdir / "RQ1A_installed_capacity_lines.png", dpi=300, bbox_inches="tight")

    plt.show()


def plot_rq1_b_capacity_and_generation(caps_long: pd.DataFrame, cba_long: pd.DataFrame, outdir: Path | None = None):
    sns.set_theme(style="whitegrid")

    cap = build_capacity_elec(caps_long)
    gen = build_generation_elec(cba_long)
    if cap.empty or gen.empty:
        print("RQ1.B: capacity or generation data missing.")
        return

    techs = sorted(set(cap["tech"].unique()) | set(gen["tech"].unique()))
    cmap = color_map_for(techs)

    stack_order = [t for t in ["Hydropower", "Geothermal", "Onshore wind", "Solar PV"] if t in techs] + [
        t for t in techs if t not in {"Hydropower", "Geothermal", "Onshore wind", "Solar PV"}
    ]

    fig, axes = plt.subplots(2, 1, figsize=(18, 11), sharex=True)
    fig.subplots_adjust(left=0.08, right=0.75, top=0.98, bottom=0.12, hspace=0.06)

    year_gap = YEAR_GAP * 0.5

    stacked_grouped_bars(
        ax=axes[0],
        df_long=cap,
        value_col="value",
        category_col="tech",
        y_label="Installed electricity generation capacity (GW)",
        years=YEAR_INTS,
        scen_order=SCEN_ORDER,
        cmap=cmap,
        cat_order=stack_order,
        dx_in=DX_IN,
        year_gap=year_gap,
        year_pad=YEAR_PAD,
        show_scen_labels=False,
    )

    stacked_grouped_bars(
        ax=axes[1],
        df_long=gen,
        value_col="value_twh",
        category_col="tech",
        y_label="Electricity generation (TWh)",
        years=YEAR_INTS,
        scen_order=SCEN_ORDER,
        cmap=cmap,
        cat_order=stack_order,
        dx_in=DX_IN,
        year_gap=year_gap,
        year_pad=YEAR_PAD,
        show_scen_labels=True,
    )

    shared_legend(
        fig=fig,
        labels=stack_order,
        cmap=cmap,
        title="Installed electricity generation capacity (GW)",
        x=0.78,
        y=0.5,
    )

    if EXPORT_FIGURES and outdir is not None:
        outdir.mkdir(parents=True, exist_ok=True)
        fig.savefig(outdir / "RQ1B_capacity_and_generation.png", dpi=300, bbox_inches="tight")

    plt.show()


def plot_rq1_c_fuelconv_and_supply(caps_long: pd.DataFrame, cba_long: pd.DataFrame, outdir: Path | None = None):
    sns.set_theme(style="whitegrid")

    cap = build_fuelconv_caps(caps_long)
    sup = build_fuel_supply_groups(cba_long)

    if cap.empty or sup.empty:
        print("RQ1.C: fuel conversion capacity or supply data missing.")
        return

    cap_order = [t for t in FUEL_CONV_STACK_ORDER if t in cap["tech"].unique().tolist()]

    carrier_present = sorted(sup["fuel_group"].unique().tolist())
    infer = colormaps.get_cmap("inferno")
    n = max(1, len(carrier_present))
    base_colors = [infer(i / max(1, (n - 1))) for i in range(n)]
    carrier_cmap = {c: base_colors[i] for i, c in enumerate(carrier_present)}

    cap_cmap = {
        "Methanizer": carrier_cmap.get("e-CH4", DEFAULT_COLOR),
        "FTropschSyn": carrier_cmap.get("e-FTL fuels", DEFAULT_COLOR),
        "Electrolyser": carrier_cmap.get("e-hydrogen", DEFAULT_COLOR),
        "DAC": "#8B5A2B",
    }

    fig, axes = plt.subplots(2, 1, figsize=(18, 11), sharex=True)
    fig.subplots_adjust(left=0.08, right=0.75, top=0.98, bottom=0.12, hspace=0.06)

    year_gap = YEAR_GAP * 0.5

    stacked_grouped_bars(
        ax=axes[0],
        df_long=cap,
        value_col="value",
        category_col="tech",
        y_label="Total installed capacity for fuel conversion (GW_output)",
        years=YEAR_INTS,
        scen_order=SCEN_ORDER,
        cmap=cap_cmap,
        cat_order=cap_order,
        dx_in=DX_IN,
        year_gap=year_gap,
        year_pad=YEAR_PAD,
        show_scen_labels=False,
    )

    stacked_grouped_bars(
        ax=axes[1],
        df_long=sup,
        value_col="value_twh",
        category_col="fuel_group",
        y_label="Supply of fuels and chemicals (TWh)",
        years=YEAR_INTS,
        scen_order=SCEN_ORDER,
        cmap=carrier_cmap,
        cat_order=carrier_present,
        dx_in=DX_IN,
        year_gap=year_gap,
        year_pad=YEAR_PAD,
        show_scen_labels=True,
    )

    handles0 = [plt.Line2D([0], [0], color=cap_cmap.get(k, DEFAULT_COLOR), lw=10) for k in cap_order]
    leg0 = fig.legend(
        handles0,
        cap_order,
        title="Total installed capacity for fuel conversion (GW_output)",
        loc="center left",
        bbox_to_anchor=(0.78, 0.78),
        frameon=True,
        borderpad=1.2,
        labelspacing=0.8,
        handlelength=1.8,
        handletextpad=0.8,
    )
    leg0.get_frame().set_linewidth(1.2)

    handles1 = [plt.Line2D([0], [0], color=carrier_cmap.get(k, DEFAULT_COLOR), lw=10) for k in carrier_present]
    leg1 = fig.legend(
        handles1,
        carrier_present,
        title="Supply of fuels and chemicals (TWh)",
        loc="center left",
        bbox_to_anchor=(0.78, 0.22),
        frameon=True,
        borderpad=1.2,
        labelspacing=0.8,
        handlelength=1.8,
        handletextpad=0.8,
    )
    leg1.get_frame().set_linewidth(1.2)

    if EXPORT_FIGURES and outdir is not None:
        outdir.mkdir(parents=True, exist_ok=True)
        fig.savefig(outdir / "RQ1C_fuelconv_and_supply.png", dpi=300, bbox_inches="tight")

    plt.show()


def plot_rq3_a_costs_and_lcoe(ind_long: pd.DataFrame, ind_det_long: pd.DataFrame, cba_long: pd.DataFrame, outdir: Path | None = None):
    sns.set_theme(style="whitegrid")

    costs = build_cost_long_from_detailed(ind_det_long)
    lcoe_df = build_lcoe(ind_det_long, cba_long)

    if costs.empty or lcoe_df.empty:
        print("RQ3-A: cost or LCOE data missing.")
        return

    comp_order = ["FuelCost", "Invest", "OMFix", "SlackCost", "SpillPenalty"]
    comp_order = [c for c in comp_order if c in costs["component"].unique().tolist()]

    colors = evenly_spaced_colors_from_cmap("viridis", n=len(comp_order), start=0.0, stop=1.0)
    cmap = {c: col for c, col in zip(comp_order, colors)}

    fig, axes = plt.subplots(2, 1, figsize=(18, 11), sharex=True)
    fig.subplots_adjust(left=0.08, right=0.75, top=0.98, bottom=0.12, hspace=0.06)

    year_gap = YEAR_GAP * 0.5

    stacked_grouped_bars(
        ax=axes[0],
        df_long=costs,
        value_col="value_bEUR",
        category_col="component",
        y_label="Total annualised cost (billion EUR)",
        years=YEAR_INTS,
        scen_order=SCEN_ORDER,
        cmap=cmap,
        cat_order=comp_order,
        dx_in=DX_IN,
        year_gap=year_gap,
        year_pad=YEAR_PAD,
        show_scen_labels=False,
    )

    for scen in SCEN_ORDER:
        df_s = lcoe_df[lcoe_df["scenario"] == scen]
        if df_s.empty:
            continue
        df_s = df_s.sort_values("year")
        axes[1].plot(df_s["year"], df_s["LCOE_EUR_MWh"], marker="o", label=scen)

    axes[1].set_xlabel("Year")
    axes[1].set_ylabel("LCOE (EUR/MWh)")
    axes[1].grid(axis="y", alpha=0.25)

    scen_handles, scen_labels = axes[1].get_legend_handles_labels()
    axes[1].legend(
        scen_handles,
        scen_labels,
        title="LCOE (EUR/MWh)",
        frameon=True,
    )

    leg = shared_legend(
        fig=fig,
        labels=comp_order,
        cmap=cmap,
        title="Total annualised cost (billion EUR)",
        x=0.78,
        y=0.5,
    )

    if EXPORT_FIGURES and outdir is not None:
        outdir.mkdir(parents=True, exist_ok=True)
        fig.savefig(outdir / "RQ3A_costs_and_lcoe.png", dpi=300, bbox_inches="tight")

    plt.show()


def plot_rq3_b_co2(ind_long: pd.DataFrame, ind_det_long: pd.DataFrame, outdir: Path | None = None):
    sns.set_theme(style="whitegrid")

    total = build_total_co2(ind_long)
    contrib = build_co2_contributors_2050(ind_det_long)

    if total.empty or contrib.empty:
        print("RQ3-B: CO2 data missing.")
        return

    fig, axes = plt.subplots(2, 1, figsize=(10, 10))
    fig.subplots_adjust(left=0.10, right=0.95, top=0.95, bottom=0.10, hspace=0.25)

    for scen in SCEN_ORDER:
        df_s = total[total["scenario"] == scen]
        if df_s.empty:
            continue
        df_s = df_s.sort_values("year")
        axes[0].plot(df_s["year"], df_s["value_mt"], marker="o", label=scen)

    axes[0].set_xlabel("Year")
    axes[0].set_ylabel("Total CO2 emissions (Mt CO2)")
    axes[0].grid(axis="y", alpha=0.25)
    axes[0].legend(
        title="Total CO2 emissions (Mt CO2)",
        frameon=True,
    )

    contrib = contrib.copy()
    cats = sorted(contrib["category"].unique().tolist())
    colors = evenly_spaced_colors_from_cmap("plasma", n=len(SCEN_ORDER), start=0.0, stop=1.0)
    cmap_scen = {s: col for s, col in zip(SCEN_ORDER, colors)}

    x = np.arange(len(cats))
    width = 0.15

    for i, scen in enumerate(SCEN_ORDER):
        df_s = contrib[contrib["scenario"] == scen]
        if df_s.empty:
            continue
        vals = [df_s[df_s["category"] == c]["value_mt"].sum() for c in cats]
        axes[1].bar(x + i * width, vals, width=width, label=scen, color=cmap_scen[scen])

    axes[1].set_xticks(x + width * (len(SCEN_ORDER) - 1) / 2.0)
    axes[1].set_xticklabels(cats, rotation=45, ha="right")
    axes[1].set_ylabel("CO2 emissions in 2050 (Mt CO2)")
    axes[1].grid(axis="y", alpha=0.25)
    axes[1].legend(
        title="CO2 emissions in 2050 (Mt CO2)",
        frameon=True,
    )

    if EXPORT_FIGURES and outdir is not None:
        outdir.mkdir(parents=True, exist_ok=True)
        fig.savefig(outdir / "RQ3B_co2.png", dpi=300, bbox_inches="tight")

    plt.show()


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main():
    outdir = BASE_DIR / "comparison_outputs"

    cba_long, caps_long, ind_long, ind_det_long = load_needed_data()

    plot_rq1_a_capacity_lines(caps_long=caps_long, outdir=outdir)
    plot_rq1_b_capacity_and_generation(caps_long=caps_long, cba_long=cba_long, outdir=outdir)
    plot_rq1_c_fuelconv_and_supply(caps_long=caps_long, cba_long=cba_long, outdir=outdir)
    plot_rq3_a_costs_and_lcoe(ind_long=ind_long, ind_det_long=ind_det_long, cba_long=cba_long, outdir=outdir)
    plot_rq3_b_co2(ind_long=ind_long, ind_det_long=ind_det_long, outdir=outdir)


if __name__ == "__main__":
    main()

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib import colormaps
from remix.framework.tools.gdx import GDXEval


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

INCLUDE_2020 = True
USE_SEPARATE_FOLDERS = True
SINGLE_GDX_FOLDER = Path(r"C:\Local\REMix\remix_nz\project\GP-NT-ELEC-BIO-H2\gdxs")

BASE_YEARS = [2020, 2025, 2030, 2035, 2040, 2045, 2050]
YEAR_INTS = BASE_YEARS if INCLUDE_2020 else [y for y in BASE_YEARS if y != 2020]

SCEN_ORDER = ["GP", "NT", "ELEC+", "BIO+", "H2+"]

DX_IN = 1.0
YEAR_GAP = 1.5
YEAR_PAD = 55
BAR_WIDTH = 0.85

GEN_TECH_MIN_TWH = 0.01
CAP_TECH_MIN_GW = 0.01

FUEL_CONV_TECHS = ["Methanizer", "FTropschSyn", "Electrolyser", "DAC"]
FUEL_CONV_STACK_ORDER = [ "DAC",  "FTropschSyn","Methanizer", "Electrolyser", ]


COST_COMPONENTS = ["FuelCost", "Invest", "OMFix", "SlackCost", "SpillPenalty", "Slack_CO2"]

LEGEND_BORDER_PAD = 0.9
LEGEND_LABEL_SPACING = 0.55
LEGEND_FRAME_ALPHA = 1.0

EXPORT_FIGURES = False

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

def load_results_flexible(scenario: str) -> GDXEval | None:
    """Load GDX based on USE_SEPARATE_FOLDERS toggle."""
    if USE_SEPARATE_FOLDERS:
        # Original structure
        case_dir = next((d for d in CASE_DIRS if scenario_name(d) == scenario), None)
        if case_dir is None:
            print(f"No case_dir found for scenario: {scenario}")
            return None
        gdx_path = BASE_DIR / case_dir / "result" / f"{case_dir}.gdx"
    else:
        # Single folder - find matching GDX file
        candidate_files = list(SINGLE_GDX_FOLDER.glob(f"*{scenario}*.gdx"))
        if not candidate_files:
            print(f"No GDX file found for scenario '{scenario}' in {SINGLE_GDX_FOLDER}")
            return None
        gdx_path = candidate_files[0]  # Take first match
    
    if not gdx_path.is_file():
        print(f"Missing: {gdx_path}")
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
            f"Could not find a year column. Tried {year_col_candidates}. Columns: {list(df.columns)}"
        )
    df["year"] = pd.to_numeric(df[year_col], errors="coerce")
    df = df[df["year"].isin(BASE_YEARS)].copy()
    df["year"] = df["year"].astype(int)
    return df


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
        ax.set_xticklabels(
            scen_labels,
            rotation=90,
            va="top",
            fontsize=9,
            fontstretch="condensed",
        )
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


def legend_math_subscripts(s: str) -> str:
    s = s.replace("CO_2", r"$\mathrm{CO_2}$")
    s = s.replace("CH_4", r"$\mathrm{CH_4}$")
    s = s.replace("H_2", r"$\mathrm{H_2}$")
    s = s.replace("GW_{output}", r"$\mathrm{GW_{output}}$")
    return s


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
    # fuel conversion naming tweaks for legends
    if tech == "Electrolyser":
        return "Electrolyser"
    if tech == "Methanizer":
        return "Methaniser"
    if tech == "FTropschSyn":
        return "Fischer-Tropsch syn."
    if tech == "DAC":
        return "Direct air capture"
    return tech.replace("_", " ")


def color_map_for(categories: list[str]) -> dict[str, str]:
    return {c: technology_colors.get(c, DEFAULT_COLOR) for c in categories}


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
    years = list(years)
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
        ax.bar(
            x,
            vals,
            bottom=bottom,
            width=BAR_WIDTH,
            color=cmap.get(c, DEFAULT_COLOR),
            edgecolor="none",
        )
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

    # tighten x-limits to reduce gap after last bar
    if len(x) > 0:
        halfw = BAR_WIDTH / 2.0
        ax.set_xlim(x[0] - halfw, x[-1] + halfw)

    return piv.columns.tolist()

def filter_to_global(df: pd.DataFrame, node_cols=("accNodesModel", "nodesModel")) -> pd.DataFrame:
    df = df.copy()
    node_col = next((c for c in node_cols if c in df.columns), None)
    if node_col is None:
        return df
    s = df[node_col].astype(str).str.lower()
    if (s == "global").any():
        return df[s == "global"].copy()
    print(f"WARNING: no {node_col}=='global' rows; aggregating across {node_col}")
    return df

def shared_legend(fig, labels: list[str], cmap: dict[str, str], title: str, x: float, y: float):
    handles = [plt.Line2D([0], [0], color=cmap.get(l, DEFAULT_COLOR), lw=10) for l in labels]
    leg = fig.legend(
        handles,
        labels,
        title=title,
        loc="center left",
        bbox_to_anchor=(x, y),
        frameon=True,
        borderpad=LEGEND_BORDER_PAD,
        labelspacing=LEGEND_LABEL_SPACING,
        handlelength=1.8,
        handletextpad=0.8,
        framealpha=LEGEND_FRAME_ALPHA,
    )
    leg.get_frame().set_linewidth(1.2)
    return leg

def load_needed_data():
    """Load data using folder toggle."""
    if not USE_SEPARATE_FOLDERS and not SINGLE_GDX_FOLDER.exists():
        raise FileNotFoundError(f"Single GDX folder not found: {SINGLE_GDX_FOLDER}")
    
    cba_list, caps_list, ind_list, ind_det_list = [], [], [], []
    
    # Load for each scenario (handles both folder structures)
    for scenario in SCEN_ORDER:
        res = load_results_flexible(scenario)
        if res is None:
            continue

        cba = read_symbol(res, "commodity_balance_annual")
        if cba is not None:
            df = ensure_year_col(cba.reset_index(), year_col_candidates=("accYears",))
            df["scenario"] = scenario
            cba_list.append(df)

        caps = read_symbol(res, "converter_caps")
        if caps is not None:
            df = ensure_year_col(caps.reset_index(), year_col_candidates=("accYears",))
            df["scenario"] = scenario
            caps_list.append(df)

        ind = read_symbol(res, "indicator_accounting")
        if ind is not None:
            df = ensure_year_col(ind.reset_index(), year_col_candidates=("accYears",))
            df["scenario"] = scenario
            ind_list.append(df)

        ind_det = read_symbol(res, "indicator_accounting_detailed")
        if ind_det is not None:
            df = ensure_year_col(ind_det.reset_index(), year_col_candidates=("years", "accYears"))
            df["scenario"] = scenario
            ind_det_list.append(df)
    
    if not cba_list or not caps_list:
        raise RuntimeError("No data loaded. Check folder paths.")
    
    # Rest unchanged from your original function...
    cba_long = pd.concat(cba_list, ignore_index=True)
    caps_long = pd.concat(caps_list, ignore_index=True)
    ind_long = pd.concat(ind_list, ignore_index=True) if ind_list else pd.DataFrame()
    ind_det_long = pd.concat(ind_det_list, ignore_index=True) if ind_det_list else pd.DataFrame()
    
    cba_long["scenario"] = pd.Categorical(cba_long["scenario"], categories=SCEN_ORDER, ordered=True)
    caps_long["scenario"] = pd.Categorical(caps_long["scenario"], categories=SCEN_ORDER, ordered=True)
    if not ind_long.empty:
        ind_long["scenario"] = pd.Categorical(ind_long["scenario"], categories=SCEN_ORDER, ordered=True)
    if not ind_det_long.empty:
        ind_det_long["scenario"] = pd.Categorical(ind_det_long["scenario"], categories=SCEN_ORDER, ordered=True)
    
    return cba_long, caps_long, ind_long, ind_det_long


def build_capacity_grouped(caps_long: pd.DataFrame) -> pd.DataFrame:
    df = caps_long.copy()
    df = filter_to_global(df)
    df = df[
        (df["capType"] == "total")
        & (df["commodity"] == "Elec")
        & (df["year"].isin(YEAR_INTS))
        & (df["value"] > 0)
    ].copy()

    df["tech"] = df["techs"].astype(str)
    df = df[df["tech"] != "Battery"].copy()
    df = df[df["value"].abs() >= CAP_TECH_MIN_GW].copy()
    df["tech"] = df["tech"].map(tech_group_pretty)

    out = df.groupby(["year", "scenario", "tech"], as_index=False, observed=False)["value"].sum()
    return out


def build_generation_grouped(cba_long: pd.DataFrame) -> pd.DataFrame:
    df = cba_long.copy()
    df = filter_to_global(df)

    df = df[
        (df["commodity"] == "Elec")
        & (df["balanceType"] == "net")
        & (df["year"].isin(YEAR_INTS))
        & (df["value"] > 0)
    ].copy()

    df["tech"] = df["techs"].astype(str)
    df["tech"] = df["tech"].map(tech_group_pretty)
    df["value_twh"] = df["value"] / 1000.0

    totals = df.groupby("tech", observed=False)["value_twh"].sum()
    keep = totals[totals >= GEN_TECH_MIN_TWH].index.tolist()
    df = df[df["tech"].isin(keep)].copy()

    out = df.groupby(["year", "scenario", "tech"], as_index=False, observed=False)["value_twh"].sum()
    return out


def build_fuelconv_caps(caps_long: pd.DataFrame) -> pd.DataFrame:
    df = caps_long.copy()
    df = filter_to_global(df)
    df = df[
        (df["capType"] == "total")
        & (df["year"].isin(YEAR_INTS))
        & (df["value"] > 0)
        & (df["techs"].astype(str).isin(FUEL_CONV_TECHS))
    ].copy()

    df = df[df["value"].abs() > 0.01].copy()
    df["tech"] = df["techs"].astype(str)
    df["tech"] = df["tech"].map(tech_group_pretty)
    return df[["year", "scenario", "tech", "value"]]


def build_fuel_supply(cba_long: pd.DataFrame) -> pd.DataFrame:
    df = cba_long.copy()
    df = filter_to_global(df)

    df = df[
        (df["year"].isin(YEAR_INTS))
        & (df["balanceType"] == "net")
        & (df["value"] > 0)
    ].copy()

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

    def supply_group(row) -> str | None:
        t = str(row["techs"])
        comm = str(row["commodity"])

        if comm == "e-CH4":
            return "Gas (e-methane)"
        if comm == "REfuel":
            return "Liquid fuel (e-FTL)"
        if comm == "H2":
            return "e-Hydrogen"

        if t in ["FuelImport_Bio_LF", "FuelImport_Biofuel"]:
            return "Liquid fuel (bio)"
        if t == "FuelImport_Bio_gas":
            return "Gas (bio)"
        if t in ["FuelImport_CH4", "FuelImport_Fossil_CH4"]:
            return "Gas (fossil)"
        if t == "FuelImport_Coal":
            return "Solid (coal)"
        if t in ["FuelImport_Diesel", "FuelImport_Fossil_LF"]:
            return "Liquid fuel (fossil)"

        return None

    df["carrier"] = df.apply(supply_group, axis=1)
    df = df[df["carrier"].notna()].copy()
    df["value_twh"] = df["value"] / 1000.0
    df = df[df["value_twh"].abs() > 0.01].copy()

    out = df.groupby(["year", "scenario", "carrier"], as_index=False, observed=False)["value_twh"].sum()
    return out


def drop_doublecounted_fuel_import_fuelcost(df_cost: pd.DataFrame) -> pd.DataFrame:
    df = df_cost.copy()
    if "indicator" not in df.columns or "techs" not in df.columns:
        raise KeyError(f"Need columns ['indicator','techs']. Got: {list(df.columns)}")
    is_fuelcost = df["indicator"].astype(str) == "FuelCost"
    is_import = df["techs"].astype(str).str.startswith("FuelImport_", na=False)
    return df[~(is_fuelcost & is_import)].copy()


def build_cost_long_from_detailed_full(ind_det_long: pd.DataFrame) -> pd.DataFrame:
    """
    Full-system annualised cost:
    - include: FuelCost, Invest, OMFix, SlackCost, SpillPenalty (all sectors, all fuels)
    - drop slack techs (Elec_slack, SlackFuel_*)
    - aggregate across all nodesModel/techs.
    """
    if ind_det_long.empty:
        return pd.DataFrame(columns=["year", "scenario", "component", "value_bEUR"])

    df = ind_det_long.copy()
    df = filter_to_global(df)

    df = df[df["indicator"].isin(["FuelCost", "Invest", "OMFix", "SlackCost", "SpillPenalty"])].copy()

    tech = df["techs"].astype(str)
    slack_mask = (tech == "Elec_slack") | tech.str.startswith("SlackFuel_", na=False)
    df = df[~slack_mask].copy()

    df = df[df["year"].isin(YEAR_INTS)].copy()

    #df = drop_doublecounted_fuel_import_fuelcost(df)

    label_map = {
        "FuelCost": "Fuel cost",
        "Invest": "Investment",
        "OMFix": "OMFix",
        "SlackCost": "Slack cost",
        "SpillPenalty": "Spill penalty",
    }
    df["component"] = df["indicator"].map(lambda x: label_map.get(str(x), str(x)))

    agg = df.groupby(["year", "scenario", "component"], as_index=False, observed=False)["value"].sum()
    agg["value_bEUR"] = agg["value"] / 1000.0
    agg = agg[agg["value_bEUR"].abs() > 0.01].copy()
    return agg[["year", "scenario", "component", "value_bEUR"]]


def evenly_spaced_colors_from_cmap(cmap_name: str, n: int, start: float = 0.0, stop: float = 1.0):
    cmap = colormaps.get_cmap(cmap_name)
    return [cmap(x) for x in np.linspace(start, stop, n)]


def build_lcoe_components(ind_det_long: pd.DataFrame, cba_long: pd.DataFrame) -> pd.DataFrame:
    if ind_det_long.empty:
        return pd.DataFrame(columns=["year", "scenario", "component", "eur_per_mwh"])

    df = ind_det_long.copy()
    df = df[df["year"].isin(YEAR_INTS)].copy()
    df = df[df["indicator"].isin(["FuelCost", "Invest", "OMFix"])].copy()

    tech = df["techs"].astype(str)
    ind = df["indicator"].astype(str)

    elec_gen_techs = {
        "CCGT",
        "GT",
        "H2_CCGT",
        "Hydro",
        "OCGT",
        "Thermal_Bio",
        "Thermal_Coal",
        "Thermal_Diesel",
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
        "HV",
        "Battery",
    }

    demand_prefixes = ("HeatDemand_", "TranspDemand_", "Elec_demand")
    slack_mask = (tech == "Elec_slack") | tech.str.startswith("SlackFuel_", na=False)

    df = df[~tech.str.startswith(demand_prefixes, na=False) & ~slack_mask].copy()

    elec_fuel_imports = {
        "FuelImport_CH4",
        "FuelImport_Coal",
        "FuelImport_Biofuel",
        "FuelImport_Diesel",
    }

    gen_cost_mask = tech.isin(elec_gen_techs) & ind.isin(["Invest", "OMFix"])
    fuel_cost_mask = tech.isin(elec_fuel_imports) & (ind == "FuelCost")

    df_num = df[gen_cost_mask | fuel_cost_mask].copy()

    label_map = {
        "FuelCost": "Fuel cost",
        "Invest": "Investment",
        "OMFix": "OMFix",
    }
    df_num["component"] = df_num["indicator"].map(lambda x: label_map.get(str(x), str(x)))

    df_num["cost_eur"] = df_num["value"] * 1e6
    cost_agg = df_num.groupby(
        ["year", "scenario", "component"], as_index=False, observed=False
    )["cost_eur"].sum()

    gen = cba_long.copy()
    if "accNodesModel" in gen.columns and (gen["accNodesModel"] == "global").any():
        gen = gen[gen["accNodesModel"] == "global"].copy()

    gen = gen[
        (gen["commodity"] == "Elec")
        & (gen["balanceType"] == "net")
        & (gen["year"].isin(YEAR_INTS))
        & (gen["value"] > 0)
    ].copy()

    gen_tech = gen["techs"].astype(str)
    gen = gen[
        ~gen_tech.str.startswith(demand_prefixes, na=False)
        & ~gen_tech.str.contains("Battery", case=False, regex=False)
    ].copy()

    gen_gwh = gen.groupby(["year", "scenario"], as_index=False, observed=False)["value"].sum()
    gen_gwh = gen_gwh.rename(columns={"value": "gen_gwh"})
    gen_gwh["gen_mwh"] = gen_gwh["gen_gwh"] * 1000.0

    merged = cost_agg.merge(
        gen_gwh[["year", "scenario", "gen_mwh"]], on=["year", "scenario"], how="left"
    )
    merged["eur_per_mwh"] = (merged["cost_eur"] / merged["gen_mwh"]).replace(
        [np.inf, -np.inf], np.nan
    )

    out = merged.groupby(
        ["year", "scenario", "component"], as_index=False, observed=False
    )["eur_per_mwh"].sum()
    out = out[out["eur_per_mwh"].abs() > 0.01].copy()

    return out


def plot_capacity_and_generation_grouped(
    caps_long: pd.DataFrame, cba_long: pd.DataFrame, outdir: Path | None = None
):
    sns.set_theme(style="whitegrid")

    cap = build_capacity_grouped(caps_long)

    tmp = caps_long.copy()
    
    #print("converter_caps columns:", tmp.columns.tolist())
    if "nodesModel" in tmp.columns:
        print(tmp["nodesModel"].astype(str).value_counts().head(20))

    gen = build_generation_grouped(cba_long)

    if cap.empty or gen.empty:
        print("RQ1.B: capacity or generation data missing.")
        return

    techs = sorted(set(cap["tech"].unique()) | set(gen["tech"].unique()))
    cmap = color_map_for(techs)

    stack_order = (
        [t for t in ["Hydropower", "Geothermal", "Onshore wind", "Solar PV"] if t in techs]
        + [t for t in techs if t not in {"Hydropower", "Geothermal", "Onshore wind", "Solar PV"}]
    )

    fig, axes = plt.subplots(2, 1, figsize=(7.5, 7.5), sharex=True)
    fig.subplots_adjust(left=0.10, right=0.74, top=0.97, bottom=0.13, hspace=0.05)

    year_gap = YEAR_GAP * 0.5

    stacked_grouped_bars(
        ax=axes[0],
        df_long=cap,
        value_col="value",
        category_col="tech",
        y_label="Installed electricity generation\ncapacity (GW)",
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

    shared_legend(fig, stack_order, cmap, title="Technology", x=0.77, y=0.5)

    if EXPORT_FIGURES and outdir is not None:
        outdir.mkdir(parents=True, exist_ok=True)
        fig.savefig(outdir / "RQ1B_capacity_and_generation.png", dpi=300, bbox_inches="tight")

    plt.show()


def plot_fuelconv_and_supply(
    caps_long: pd.DataFrame, cba_long: pd.DataFrame, outdir: Path | None = None
):
    sns.set_theme(style="whitegrid")

    cap = build_fuelconv_caps(caps_long)
    sup = build_fuel_supply(cba_long)

    if cap.empty or sup.empty:
        print("RQ1.C: fuel conversion capacity or supply data missing.")
        return

    cap_order = [t for t in FUEL_CONV_STACK_ORDER if t in cap["tech"].unique().tolist()]

    # YlOrBr-style sequence for fuel conversion techs
    n_cap = len(cap_order) if cap_order else 1
    ylbr_colors = ylbr_colors[::-1]
    ylbr_colors = evenly_spaced_colors_from_cmap("YlOrBr", n_cap, start=0.2, stop=0.9)
    
    cap_cmap = {tech: col for tech, col in zip(cap_order, ylbr_colors)}

    all_carriers = sorted(sup["carrier"].unique().tolist())
    carrier_order = [
        "Solid (coal)",
        "Liquid fuel (fossil)",
        "Gas (fossil)",
        "Liquid fuel (bio)",
        "Gas (bio)",
        "Liquid fuel (e-FTL)",
        "Gas (e-methane)",
        "e-Hydrogen",
    ]
    carrier_order = [c for c in carrier_order if c in all_carriers]

    infer = colormaps.get_cmap("inferno")
    frac_map = {
        "Solid (coal)": 0.05,
        "Liquid fuel (fossil)": 0.15,
        "Gas (fossil)": 0.25,
        "Liquid fuel (bio)": 0.40,
        "Gas (bio)": 0.50,
        "Liquid fuel (e-FTL)": 0.65,
        "Gas (e-methane)": 0.80,
        "e-Hydrogen": 0.95,
    }
    carrier_colors = {c: infer(frac_map.get(c, 0.5)) for c in carrier_order}
    mid_col = infer(0.5)
    for c in all_carriers:
        if c not in carrier_colors:
            carrier_colors[c] = mid_col

    fig, axes = plt.subplots(2, 1, figsize=(7.5, 7.5), sharex=True)
    fig.subplots_adjust(left=0.10, right=0.74, top=0.97, bottom=0.13, hspace=0.05)

    year_gap = YEAR_GAP * 0.5

    stacked_grouped_bars(
        ax=axes[0],
        df_long=cap,
        value_col="value",
        category_col="tech",
        y_label="Total installed capacity for fuel conversion\n"
        r"$\mathrm{GW_{output}}$",
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
        category_col="carrier",
        y_label="Supply of fuels and chemicals (TWh)",
        years=YEAR_INTS,
        scen_order=SCEN_ORDER,
        cmap=carrier_colors,
        cat_order=carrier_order,
        dx_in=DX_IN,
        year_gap=year_gap,
        year_pad=YEAR_PAD,
        show_scen_labels=True,
    )

    LEGEND_X = 0.75
    LEGEND_Y_TOP = 0.80
    LEGEND_Y_BOT = 0.23

    handles0 = [plt.Line2D([0], [0], color=cap_cmap.get(k, DEFAULT_COLOR), lw=10) for k in cap_order]
    leg0 = fig.legend(
        handles0,
        cap_order,
        title="Technology",
        loc="center left",
        bbox_to_anchor=(LEGEND_X, LEGEND_Y_TOP),
        frameon=True,
        borderpad=1.0,
        labelspacing=0.6,
        handlelength=1.6,
        handletextpad=0.7,
    )
    leg0.get_frame().set_linewidth(1.2)

    carrierlabels = [legend_math_subscripts(x).replace(" for ", "\nfor ") for x in carrier_order]
    handles1 = [
        plt.Line2D([0], [0], color=carrier_colors.get(k, DEFAULT_COLOR), lw=10)
        for k in carrier_order
    ]
    leg1 = fig.legend(
        handles1,
        carrierlabels,
        title="Carrier",
        loc="center left",
        bbox_to_anchor=(LEGEND_X, LEGEND_Y_BOT),
        frameon=True,
        borderpad=1.0,
        labelspacing=0.6,
        handlelength=1.6,
        handletextpad=0.7,
    )
    leg1.get_frame().set_linewidth(1.2)

    if EXPORT_FIGURES and outdir is not None:
        outdir.mkdir(parents=True, exist_ok=True)
        fig.savefig(outdir / "RQ1C_fuelconv_and_supply.png", dpi=300, bbox_inches="tight")

    plt.show()


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
    # fuel conversion naming tweaks for legends
    if tech == "Electrolyser":
        return "Electrolyser"
    if tech == "Methanizer":
        return "Methaniser"
    if tech == "FTropschSyn":
        return "Fischer-Tropsch syn."
    if tech == "DAC":
        return "Direct air capture"
    return tech.replace("_", " ")


def build_fuelconv_caps(caps_long: pd.DataFrame) -> pd.DataFrame:
    df = caps_long.copy()
    
    if "accNodesModel" in df.columns:
        if (df["accNodesModel"].astype(str) == "global").any():
            df = df[df["accNodesModel"].astype(str) == "global"].copy()
        else:
            print("WARNING: no accNodesModel=='global' rows; aggregating across accNodesModel")
              




    df = df[
        (df["capType"] == "total")
        & (df["year"].isin(YEAR_INTS))
        & (df["value"] > 0)
        & (df["techs"].astype(str).isin(FUEL_CONV_TECHS))
    ].copy()

    df = df[df["value"].abs() > 0.01].copy()
    df["tech"] = df["techs"].astype(str)
    df["tech"] = df["tech"].map(tech_group_pretty)
    return df[["year", "scenario", "tech", "value"]]


def plot_fuelconv_and_supply(
    caps_long: pd.DataFrame, cba_long: pd.DataFrame, outdir: Path | None = None
):
    sns.set_theme(style="whitegrid")

    cap = build_fuelconv_caps(caps_long)
    sup = build_fuel_supply(cba_long)

    if cap.empty or sup.empty:
        print("RQ1.C: fuel conversion capacity or supply data missing.")
        return

    cap_order = [t for t in FUEL_CONV_STACK_ORDER if t in {"Electrolyser", "Methanizer", "FTropschSyn", "DAC"}]
    present = cap["tech"].unique().tolist()
    pretty_map = {
        "Electrolyser": "Electrolyser",
        "Methanizer": "Methaniser",
        "FTropschSyn": "Fischer-Tropsch syn.",
        "DAC": "Direct air capture",
    }
    cap_order = [pretty_map[t] for t in cap_order if pretty_map[t] in present]

    n_cap = len(cap_order) if cap_order else 1
    ylbr_colors = evenly_spaced_colors_from_cmap("YlOrBr", n_cap, start=0.2, stop=0.9)
    ylbr_colors = ylbr_colors[::-1]
    cap_cmap = {tech: col for tech, col in zip(cap_order, ylbr_colors)}

    all_carriers = sorted(sup["carrier"].unique().tolist())
    carrier_order = [
        "Solid (coal)",
        "Liquid fuel (fossil)",
        "Gas (fossil)",
        "Liquid fuel (bio)",
        "Gas (bio)",
        "Liquid fuel (e-FTL)",
        "Gas (e-methane)",
        "e-Hydrogen",
    ]
    carrier_order = [c for c in carrier_order if c in all_carriers]

    infer = colormaps.get_cmap("inferno")
    frac_map = {
        "Solid (coal)": 0.05,
        "Liquid fuel (fossil)": 0.15,
        "Gas (fossil)": 0.25,
        "Liquid fuel (bio)": 0.40,
        "Gas (bio)": 0.50,
        "Liquid fuel (e-FTL)": 0.65,
        "Gas (e-methane)": 0.80,
        "e-Hydrogen": 0.95,
    }
    carrier_colors = {c: infer(frac_map.get(c, 0.5)) for c in carrier_order}
    mid_col = infer(0.5)
    for c in all_carriers:
        if c not in carrier_colors:
            carrier_colors[c] = mid_col

    fig, axes = plt.subplots(2, 1, figsize=(7.5, 7.5), sharex=True)
    fig.subplots_adjust(left=0.10, right=0.74, top=0.97, bottom=0.13, hspace=0.05)

    year_gap = YEAR_GAP * 0.5

    stacked_grouped_bars(
        ax=axes[0],
        df_long=cap,
        value_col="value",
        category_col="tech",
        y_label="Total installed capacity for fuel conversion\n"
        r"(GW$_{output}$)",
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
        category_col="carrier",
        y_label="Supply of fuels and chemicals (TWh)",
        years=YEAR_INTS,
        scen_order=SCEN_ORDER,
        cmap=carrier_colors,
        cat_order=carrier_order,
        dx_in=DX_IN,
        year_gap=year_gap,
        year_pad=YEAR_PAD,
        show_scen_labels=True,
    )

    LEGEND_X = 0.72
    LEGEND_Y_TOP = 0.80
    LEGEND_Y_BOT = 0.23

    handles0 = [plt.Line2D([0], [0], color=cap_cmap.get(k, DEFAULT_COLOR), lw=10) for k in cap_order]
    leg0 = fig.legend(
        handles0,
        cap_order,
        title="Technology",
        loc="center left",
        bbox_to_anchor=(LEGEND_X, LEGEND_Y_TOP),
        frameon=True,
        borderpad=1.0,
        labelspacing=0.6,
        handlelength=1.6,
        handletextpad=0.7,
    )
    leg0.get_frame().set_linewidth(1.2)

    carrierlabels = [legend_math_subscripts(x).replace(" for ", "\nfor ") for x in carrier_order]
    handles1 = [
        plt.Line2D([0], [0], color=carrier_colors.get(k, DEFAULT_COLOR), lw=10)
        for k in carrier_order
    ]
    leg1 = fig.legend(
        handles1,
        carrierlabels,
        title="Carrier",
        loc="center left",
        bbox_to_anchor=(LEGEND_X, LEGEND_Y_BOT),
        frameon=True,
        borderpad=1.0,
        labelspacing=0.6,
        handlelength=1.6,
        handletextpad=0.7,
    )
    leg1.get_frame().set_linewidth(1.2)

    if EXPORT_FIGURES and outdir is not None:
        outdir.mkdir(parents=True, exist_ok=True)
        fig.savefig(outdir / "RQ1C_fuelconv_and_supply.png", dpi=300, bbox_inches="tight")

    plt.show()


def plot_costs_and_lcoe(
    ind_long: pd.DataFrame,
    ind_det_long: pd.DataFrame,
    cba_long: pd.DataFrame,
    outdir: Path | None = None,
):
    sns.set_theme(style="whitegrid")

    costs = build_cost_long_from_detailed_full(ind_det_long)
    lcoe_comp = build_lcoe_components(ind_det_long, cba_long)

    if costs.empty:
        print("No indicator_accounting_detailed loaded; skipping cost/LCOE plots.")
        return

    cost_comp_order = ["Fuel cost", "Investment", "OMFix", "Slack cost", "Spill penalty"]
    cost_comp_order = [
        c for c in cost_comp_order if c in costs["component"].unique().tolist()
    ]

    lcoe_comp_order = ["Fuel cost", "Investment", "OMFix"]
    lcoe_comp_order = [
        c for c in lcoe_comp_order if not lcoe_comp.empty and c in lcoe_comp["component"].unique().tolist()
    ]

    # revert viridis to full range
    ncols = max(len(cost_comp_order), len(lcoe_comp_order or []))
    colors = evenly_spaced_colors_from_cmap("viridis", n=max(ncols, 1), start=0.0, stop=1.0)
    base_labels = (
        cost_comp_order
        if len(cost_comp_order) >= len(lcoe_comp_order or [])
        else lcoe_comp_order
    )
    cmap = {c: col for c, col in zip(base_labels, colors)}

    fig1, ax1 = plt.subplots(1, 1, figsize=(7.5, 4.2))
    fig1.subplots_adjust(left=0.10, right=0.74, top=0.97, bottom=0.20)

    stacked_grouped_bars(
        ax=ax1,
        df_long=costs,
        value_col="value_bEUR",
        category_col="component",
        y_label="Total annualised system cost (billion EUR)",
        years=YEAR_INTS,
        scen_order=SCEN_ORDER,
        cmap=cmap,
        cat_order=cost_comp_order,
        dx_in=DX_IN,
        year_gap=YEAR_GAP * 0.5,
        year_pad=YEAR_PAD,
        show_scen_labels=True,
    )

    shared_legend(fig1, cost_comp_order, cmap, title="Cost component", x=0.78, y=0.5)

    if EXPORT_FIGURES and outdir is not None:
        outdir.mkdir(parents=True, exist_ok=True)
        fig1.savefig(outdir / "RQ3_total_costs.png", dpi=300, bbox_inches="tight")

    plt.show()

    if lcoe_comp.empty:
        return

    fig2, ax2 = plt.subplots(1, 1, figsize=(7.5, 4.2))
    fig2.subplots_adjust(left=0.10, right=0.74, top=0.97, bottom=0.20)

    stacked_grouped_bars(
        ax=ax2,
        df_long=lcoe_comp,
        value_col="eur_per_mwh",
        category_col="component",
        y_label="LCOE decomposition (EUR/MWh)",
        years=YEAR_INTS,
        scen_order=SCEN_ORDER,
        cmap=cmap,
        cat_order=lcoe_comp_order,
        dx_in=DX_IN,
        year_gap=YEAR_GAP * 0.5,
        year_pad=YEAR_PAD,
        show_scen_labels=True,
    )

    shared_legend(fig2, lcoe_comp_order, cmap, title="Cost component", x=0.78, y=0.5)

    if EXPORT_FIGURES and outdir is not None:
        fig2.savefig(outdir / "RQ3_lcoe_components.png", dpi=300, bbox_inches="tight")

    plt.show()

def plot_co2_grouped(ind_det_long: pd.DataFrame, outdir: Path | None = None):
    sns.set_theme(style="whitegrid")

    # simple labels; group name will be added in legend
    category_map = {
        "TranspDemand_Fossil_LF": "fossil liquid fuel",
        "TranspDemand_REfuel": "e-FTL",
        "TranspDemand_Fossil_CH4": "fossil gas",
        "TranspDemand_e-CH4": "e-methane",
        "HeatDemand_Fossil_LF": "fossil liquid fuel (heat)",
        "HeatDemand_REfuel": "e-FTL (heat)",
        "HeatDemand_Fossil_CH4": "fossil gas (heat)",
        "HeatDemand_e-CH4": "e-methane (heat)",
        "CCGT": "gas",
        "OCGT": "gas",
        "GT": "gas",
        "Thermal_Coal": "coal",
        "Thermal_Diesel": "diesel",
        "DAC": "Direct air capture",
    }

    # colours; transport blues, heat oranges from YlOrBr, power magentas, DAC teal
    ylbr = plt.cm.get_cmap("YlOrBr")
    C = {
        "fossil liquid fuel": "#1851B3",
        "e-FTL": "#3F8CF1",
        "fossil gas": "#6599ED",
        "e-methane": "#6C92BB",
        "fossil liquid fuel (heat)": ylbr(0.80),
        "e-FTL (heat)": ylbr(0.65),
        "fossil gas (heat)": ylbr(0.50),
        "e-methane (heat)": ylbr(0.35),
        "gas": "#C428D9",
        "coal": "#E95CF6",
        "diesel": "#F755CC",
        "Direct air capture": "#B3DDD7",
    }

    alias = {
        "TranspDemandFossilLF": "TranspDemand_Fossil_LF",
        "TranspDemandREfuel": "TranspDemand_REfuel",
        "TranspDemandFossilCH4": "TranspDemand_Fossil_CH4",
        "TranspDemande-CH4": "TranspDemand_e-CH4",
        "HeatDemandFossilLF": "HeatDemand_Fossil_LF",
        "HeatDemandREfuel": "HeatDemand_REfuel",
        "HeatDemandFossilCH4": "HeatDemand_Fossil_CH4",
        "HeatDemande-CH4": "HeatDemand_e-CH4",
        "ThermalCoal": "Thermal_Coal",
        "ThermalDiesel": "Thermal_Diesel",
    }

    cat_order = [
        "fossil liquid fuel",
        "e-FTL",
        "fossil gas",
        "e-methane",
        "fossil liquid fuel (heat)",
        "e-FTL (heat)",
        "fossil gas (heat)",
        "e-methane (heat)",
        "gas",
        "coal",
        "diesel",
        "Direct air capture",
    ]

    df = ind_det_long.copy()
    if "year" not in df.columns:
        for cand in ["accYears", "years"]:
            if cand in df.columns:
                df = df.rename(columns={cand: "year"})
                break
    df["year"] = pd.to_numeric(df["year"], errors="coerce").astype("Int64")
    df = df[df["year"].isin(YEAR_INTS)].copy()

    df = df[df["indicator"].astype(str) == "CO2_emission"].copy()
    df["value_mt"] = df["value"] / 1000.0

    tech_norm = df["techs"].astype(str).map(lambda t: alias.get(t, t))
    df["category"] = tech_norm.map(category_map)
    df = df[df["category"].notna()].copy()

    agg = df.groupby(["year", "scenario", "category"], as_index=False, observed=False)["value_mt"].sum()

    years = list(YEAR_INTS)
    scens = list(SCEN_ORDER)
    x = []
    year_centers = []
    pos = 0.0
    for _y in years:
        xs = [pos + i * DX_IN for i in range(len(scens))]
        x.extend(xs)
        year_centers.append(float(np.mean(xs)))
        pos = xs[-1] + DX_IN + YEAR_GAP
    x = np.array(x)
    year_centers = np.array(year_centers)

    piv = (
        agg.pivot_table(
            index=["year", "scenario"],
            columns="category",
            values="value_mt",
            aggfunc="sum",
            observed=False,
        )
        .reindex(pd.MultiIndex.from_product([years, scens], names=["year", "scenario"]))
        .fillna(0.0)
    )

    cols = [c for c in cat_order if c in piv.columns]
    piv = piv[cols]

    fig, ax = plt.subplots(1, 1, figsize=(7.5, 4.2))
    fig.subplots_adjust(left=0.10, right=0.68, top=0.97, bottom=0.22)

    bottom_pos = np.zeros(len(x))
    bottom_neg = np.zeros(len(x))

    for c in piv.columns:
        vals = piv[c].values
        color = C.get(c, "#9E9E9E")
        pos_vals = np.where(vals > 0, vals, 0.0)
        neg_vals = np.where(vals < 0, vals, 0.0)
        if np.any(pos_vals != 0.0):
            ax.bar(x, pos_vals, bottom=bottom_pos, width=BAR_WIDTH, color=color, edgecolor="none")
            bottom_pos += pos_vals
        if np.any(neg_vals != 0.0):
            ax.bar(x, neg_vals, bottom=bottom_neg, width=BAR_WIDTH, color=color, edgecolor="none")
            bottom_neg += neg_vals

    ax.grid(axis="y", alpha=0.25)
    ax.set_ylabel(r"Total CO$_{2}$ emissions by contributor (Mt CO$_{2}$)")

    scen_labels = [s for _y in years for s in scens]
    ax.set_xticks(x)
    ax.set_xticklabels(
        scen_labels,
        rotation=90,
        va="top",
        fontsize=9,
        fontstretch="condensed",
    )

    if len(x) > 0:
        halfw = BAR_WIDTH / 2.0
        ax.set_xlim(x[0] - halfw, x[-1] + halfw)

    sec = ax.secondary_xaxis("bottom")
    sec.set_xlim(ax.get_xlim())
    sec.set_xticks(year_centers)
    sec.set_xticklabels([str(y) for y in years], rotation=0)
    sec.spines["bottom"].set_visible(False)
    sec.tick_params(axis="x", pad=YEAR_PAD, length=0)

    # Legend with group titles slightly separated vertically, titles not bold, left aligned
    legend_entries = {
        "Transport": [
            ("Fossil liquid fuel", C["fossil liquid fuel"]),
            ("e-FTL", C["e-FTL"]),
            ("Fossil gas", C["fossil gas"]),
            ("e-methane", C["e-methane"]),
        ],
        "Heat": [
            ("Fossil liquid fuel", C["fossil liquid fuel (heat)"]),
            ("e-FTL", C["e-FTL (heat)"]),
            ("Fossil gas", C["fossil gas (heat)"]),
            ("e-methane", C["e-methane (heat)"]),
        ],
        "Electricity generation": [
            ("Gas", C["gas"]),
            ("Coal", C["coal"]),
            ("Diesel", C["diesel"]),
        ],
        "Carbon capture": [
            ("Direct air capture", C["Direct air capture"]),
        ],
    }

    legend_labels = []
    legend_handles = []
    group_titles = set(legend_entries.keys())

    for i, (group_title, items) in enumerate(legend_entries.items()):
        # extra blank line before every group except the first → more space between groups
        if i > 0:
            legend_labels.append("")  # blank separator
            legend_handles.append(plt.Line2D([0], [0], color="none", linewidth=0))

        # group title
        legend_labels.append(group_title)
        legend_handles.append(
            plt.Line2D([0], [0], color="none", linewidth=0)
        )

        # items
        for label, col in items:
            legend_labels.append(label)
            legend_handles.append(
                plt.Line2D([0], [0], color=col, linewidth=8)
            )

    leg = ax.legend(
        legend_handles,
        legend_labels,
        loc="center left",
        bbox_to_anchor=(1.04, 0.50), 
        frameon=True,
        borderpad=0.8,        
        labelspacing=0.35,    
        handlelength=1.8,    
        handletextpad=0.7,   
        framealpha=1.0,
        columnspacing=0.8,  
    )


    if EXPORT_FIGURES and outdir is not None:
        outdir.mkdir(parents=True, exist_ok=True)
        fig.savefig(outdir / "RQ3_co2_grouped.png", dpi=300, bbox_inches="tight")

    plt.show()


def main() -> None:
    outdir = BASE_DIR / "comparison_outputs"

    cba_long, caps_long, ind_long, ind_det_long = load_needed_data()

    plot_capacity_and_generation_grouped(caps_long=caps_long, cba_long=cba_long, outdir=outdir)
    plot_fuelconv_and_supply(caps_long=caps_long, cba_long=cba_long, outdir=outdir)
    plot_costs_and_lcoe(
        ind_long=ind_long, ind_det_long=ind_det_long, cba_long=cba_long, outdir=outdir
    )
    plot_co2_grouped(ind_det_long=ind_det_long, outdir=outdir)


if __name__ == "__main__":
    main()

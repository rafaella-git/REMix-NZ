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

YEAR_INTS = [2020, 2025, 2030, 2035, 2040, 2045, 2050]
SCEN_ORDER = ["GP", "NT", "ELEC+", "BIO+", "H2+"]

DX_IN = 1.0
YEAR_GAP = 1.5
YEAR_PAD = 55
BAR_WIDTH = 0.85

GEN_TECH_MIN_TWH = 0.01
CAP_TECH_MIN_GW = 0.01

FUEL_CONV_TECHS = ["Methanizer", "FTropschSyn", "Electrolyser", "DAC"]
FUEL_CONV_STACK_ORDER = ["Electrolyser", "Methanizer", "FTropschSyn", "DAC"]

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
            f"Could not find a year column. Tried {year_col_candidates}. Columns: {list(df.columns)}"
        )
    df["year"] = pd.to_numeric(df[year_col], errors="coerce")
    df = df[df["year"].isin(YEAR_INTS)].copy()
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
        ind_det_long["scenario"] = pd.Categorical(ind_det_long["scenario"], categories=SCEN_ORDER, ordered=True)

    return cba_long, caps_long, ind_long, ind_det_long


def build_capacity_grouped(caps_long: pd.DataFrame) -> pd.DataFrame:
    df = caps_long.copy()
    if "nodesModel" in df.columns and (df["nodesModel"] == "global").any():
        df = df[df["nodesModel"] == "global"].copy()

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
    if "accNodesModel" in df.columns and (df["accNodesModel"] == "global").any():
        df = df[df["accNodesModel"] == "global"].copy()

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
    if "accNodesModel" in df.columns and (df["accNodesModel"] == "global").any():
        df = df[df["accNodesModel"] == "global"].copy()

    df = df[
        (df["capType"] == "total")
        & (df["year"].isin(YEAR_INTS))
        & (df["value"] > 0)
        & (df["techs"].astype(str).isin(FUEL_CONV_TECHS))
    ].copy()

    df = df[df["value"].abs() > 0.01].copy()
    df["tech"] = df["techs"].astype(str)
    return df[["year", "scenario", "tech", "value"]]


def build_fuel_supply(cba_long: pd.DataFrame) -> pd.DataFrame:
    df = cba_long.copy()
    if "accNodesModel" in df.columns and (df["accNodesModel"] == "global").any():
        df = df[df["accNodesModel"] == "global"].copy()

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


def build_cost_long_from_detailed(ind_det_long: pd.DataFrame) -> pd.DataFrame:
    if ind_det_long.empty:
        return pd.DataFrame(columns=["year", "scenario", "component", "value_bEUR"])

    df = ind_det_long.copy()
    df = df[(df["year"].isin(YEAR_INTS)) & (df["indicator"].isin(COST_COMPONENTS))].copy()
    df = drop_doublecounted_fuel_import_fuelcost(df)

    label_map = {
        "FuelCost": "Fuel cost",
        "Invest": "Investment",
        "OMFix": "OMFix",
        "SlackCost": "Slack cost",
        "SpillPenalty": "Spill penalty",
        "Slack_CO2": r"CO$_{2}$ slack",
        "SystemCost": "System cost",
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
    """
    LCOE = (Invest + OMFix for electricity generators + HV + Battery
            + FuelCost for their fuels) / total positive electricity generation.
    """

    if ind_det_long.empty:
        return pd.DataFrame(columns=["year", "scenario", "component", "eur_per_mwh"])

    df = ind_det_long.copy()
    df = df[df["year"].isin(YEAR_INTS)].copy()
    df = df[df["indicator"].isin(["FuelCost", "Invest", "OMFix"])].copy()

    tech = df["techs"].astype(str)
    ind = df["indicator"].astype(str)

    # generator + network + storage techs in LCOE
    elec_gen_techs = {
        "CCGT", "GT", "H2_CCGT", "Hydro", "OCGT",
        "Thermal_Bio", "Thermal_Coal", "Thermal_Diesel",
        "pv_central_fixed", "pv_decentral",
        "wind_offshore_1", "wind_offshore_2", "wind_offshore_3", "wind_offshore_4",
        "wind_onshore_1", "wind_onshore_2", "wind_onshore_3", "wind_onshore_4",
        # system elements you now want in LCOE:
        "HV", "Battery",
    }

    demand_prefixes = ("HeatDemand_", "TranspDemand_", "Elec_demand")
    slack_mask = (tech == "Elec_slack") | tech.str.startswith("SlackFuel_", na=False)

    df = df[
        ~tech.str.startswith(demand_prefixes, na=False)
        & ~slack_mask
    ].copy()

    elec_fuel_imports = {
        "FuelImport_CH4",
        "FuelImport_Coal",
        "FuelImport_Biofuel",
        "FuelImport_Diesel",
    }

    gen_cost_mask = tech.isin(elec_gen_techs) & ind.isin(["Invest", "OMFix"])
    fuel_cost_mask = tech.isin(elec_fuel_imports) & (ind == "FuelCost")

    df_num = df[gen_cost_mask | fuel_cost_mask].copy()

    # diagnostics (keep, they’re useful)
    print("\n=== LCOE (with HV + Battery): cost summary by year, scenario, indicator ===")
    if not df_num.empty:
        sum_dbg = (
            df_num.groupby(["year", "scenario", "indicator"], observed=False)["value"]
            .sum()
            .reset_index()
            .sort_values(["scenario", "year", "indicator"])
        )
        print(sum_dbg.to_string(index=False))

    label_map = {
        "FuelCost": "Fuel cost",
        "Invest": "Investment",
        "OMFix": "OMFix",
    }
    df_num["component"] = df_num["indicator"].map(lambda x: label_map.get(str(x), str(x)))

    df_num["cost_eur"] = df_num["value"] * 1e6
    cost_agg = df_num.groupby(["year", "scenario", "component"], as_index=False, observed=False)["cost_eur"].sum()

    # denominator: same as before (global, Elec, net, value > 0, drop demands + Battery inflows)
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

    merged = cost_agg.merge(gen_gwh[["year", "scenario", "gen_mwh"]], on=["year", "scenario"], how="left")
    merged["eur_per_mwh"] = (merged["cost_eur"] / merged["gen_mwh"]).replace(
        [np.inf, -np.inf], np.nan
    )

    out = merged.groupby(["year", "scenario", "component"], as_index=False, observed=False)["eur_per_mwh"].sum()
    out = out[out["eur_per_mwh"].abs() > 0.01].copy()

    print("\n=== LCOE (with HV + Battery): components (EUR/MWh) ===")
    if not out.empty:
        print(
            out.sort_values(["scenario", "year", "component"])
            .to_string(index=False)
        )

    return out

def build_cost_long_from_detailed(ind_det_long: pd.DataFrame) -> pd.DataFrame:
    """
    Full-system annualised cost:
    - include: FuelCost, Invest, OMFix, SlackCost, SpillPenalty
    - exclude: CO2_emission (indicator only)
    - drop slack techs (Elec_slack, SlackFuel_*)
    - aggregate across all nodesModel/techs.
    """
    if ind_det_long.empty:
        return pd.DataFrame(columns=["year", "scenario", "component", "value_bEUR"])

    df = ind_det_long.copy()

    # keep all monetary indicators
    df = df[df["indicator"].isin(["FuelCost", "Invest", "OMFix", "SlackCost", "SpillPenalty"])].copy()

    # drop slack techs: Elec_slack + SlackFuel_*
    tech = df["techs"].astype(str)
    slack_mask = (tech == "Elec_slack") | tech.str.startswith("SlackFuel_", na=False)
    df = df[~slack_mask].copy()

    df = df[df["year"].isin(YEAR_INTS)].copy()

    label_map = {
        "FuelCost": "Fuel cost",
        "Invest": "Investment",
        "OMFix": "OMFix",
        "SlackCost": "Slack cost",
        "SpillPenalty": "Spill penalty",
    }
    df["component"] = df["indicator"].map(lambda x: label_map.get(str(x), str(x)))

    agg = df.groupby(["year", "scenario", "component"], as_index=False, observed=False)["value"].sum()
    agg["value_bEUR"] = agg["value"] / 1000.0  # million -> billion
    agg = agg[agg["value_bEUR"].abs() > 0.01].copy()
    return agg[["year", "scenario", "component", "value_bEUR"]]

def plot_capacity_and_generation_grouped(caps_long: pd.DataFrame, cba_long: pd.DataFrame, outdir: Path | None = None):
    sns.set_theme(style="whitegrid")

    cap = build_capacity_grouped(caps_long)
    gen = build_generation_grouped(cba_long)

    if cap.empty or gen.empty:
        print("RQ1.B: capacity or generation data missing.")
        return

    techs = sorted(set(cap["tech"].unique()) | set(gen["tech"].unique()))
    cmap = color_map_for(techs)

    stack_order = [t for t in ["Hydropower", "Geothermal", "Onshore wind", "Solar PV"] if t in techs] + [
        t for t in techs if t not in {"Hydropower", "Geothermal", "Onshore wind", "Solar PV"}
    ]

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


def plot_fuelconv_and_supply(caps_long: pd.DataFrame, cba_long: pd.DataFrame, outdir: Path | None = None):
    sns.set_theme(style="whitegrid")

    cap = build_fuelconv_caps(caps_long)
    sup = build_fuel_supply(cba_long)

    if cap.empty or sup.empty:
        print("RQ1.C: fuel conversion capacity or supply data missing.")
        return

    cap_order = [t for t in FUEL_CONV_STACK_ORDER if t in cap["tech"].unique().tolist()]

    cap_cmap = {
        "Electrolyser": "#FFF7BC",  # pale yellow (YlOrRd-ish)
        "Methanizer":   "#FEE391",  # light yellow-orange
        "FTropschSyn":  "#FEC44F",  # light orange
        "DAC":          "#FE9929",  # light orange-red
    }


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


def plot_costs_and_lcoe(ind_long: pd.DataFrame, ind_det_long: pd.DataFrame,
                        cba_long: pd.DataFrame, outdir: Path | None = None):
    sns.set_theme(style="whitegrid")

    costs = build_cost_long_from_detailed(ind_det_long)
    lcoe_comp = build_lcoe_components(ind_det_long, cba_long)

    if costs.empty:
        print("No indicator_accounting_detailed loaded; skipping cost/LCOE plots.")
        return

    # order for full-system costs (top figure)
    cost_comp_order = ["Fuel cost", "Investment", "OMFix", "Slack cost", "Spill penalty"]
    cost_comp_order = [c for c in cost_comp_order if c in costs["component"].unique().tolist()]

    # order for LCOE (bottom figure) – only components present in lcoe_comp
    lcoe_comp_order = ["Fuel cost", "Investment", "OMFix"]
    lcoe_comp_order = [c for c in lcoe_comp_order if c in lcoe_comp["component"].unique().tolist()]

    # shared colours (so Fuel/Invest/OMFix look the same in both plots)
    ncols = max(len(cost_comp_order), len(lcoe_comp_order))
    colors = evenly_spaced_colors_from_cmap("viridis", n=ncols, start=0.0, stop=1.0)
    base_labels = cost_comp_order if len(cost_comp_order) >= len(lcoe_comp_order) else lcoe_comp_order
    cmap = {c: col for c, col in zip(base_labels, colors)}

    # ---------- Figure 1: total annualised cost ----------
    fig1, ax1 = plt.subplots(1, 1, figsize=(7.5, 4.2))
    fig1.subplots_adjust(left=0.10, right=0.74, top=0.97, bottom=0.20)

    stacked_grouped_bars(
        ax=ax1,
        df_long=costs,
        value_col="value_bEUR",
        category_col="component",
        y_label="Total annualised cost (billion EUR)",
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

    # ---------- Figure 2: LCOE decomposition ----------
    if lcoe_comp.empty:
        print("LCOE components are empty; skipping LCOE plot.")
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

    # map detailed techs to categories
    category_map = {
        "TranspDemand_Fossil_LF": "Transport (fossil liquid fuel)",
        "TranspDemand_REfuel": "Transport (liquid e-fuel)",
        "TranspDemand_Fossil_CH4": "Transport (fossil methane)",
        "TranspDemand_e-CH4": "Transport (e-methane)",
        "HeatDemand_Fossil_LF": "Heat (fossil liquid fuel)",
        "HeatDemand_REfuel": "Heat (e-fuel)",
        "HeatDemand_Fossil_CH4": "Heat (fossil methane)",
        "HeatDemand_e-CH4": "Heat (e-methane)",
        "CCGT": "Electricity generation (gas)",
        "OCGT": "Electricity generation (gas)",
        "GT": "Electricity generation (gas)",
        "Thermal_Coal": "Electricity generation (coal)",
        "Thermal_Diesel": "Electricity generation (diesel)",
        "DAC": "Direct air capture",
    }

    # colours (same as before, you can tweak DAC separately)
    C = {
        "Transport (fossil liquid fuel)": "#1851B3",
        "Transport (liquid e-fuel)": "#3F8CF1",
        "Transport (fossil methane)": "#6599ED",
        "Transport (e-methane)": "#6C92BB",
        "Heat (fossil liquid fuel)": "#FF8C28",
        "Heat (e-fuel)": "#FF962E",
        "Heat (fossil methane)": "#FFBB00",
        "Heat (e-methane)": "#FECF58",
        "Electricity generation (gas)": "#C428D9",
        "Electricity generation (coal)": "#E95CF6",
        "Electricity generation (diesel)": "#F755CC",
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
        "Transport (fossil liquid fuel)",
        "Transport (liquid e-fuel)",
        "Transport (fossil methane)",
        "Transport (e-methane)",
        "Heat (fossil liquid fuel)",
        "Heat (e-fuel)",
        "Heat (fossil methane)",
        "Heat (e-methane)",
        "Electricity generation (gas)",
        "Electricity generation (coal)",
        "Electricity generation (diesel)",
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

    # kt -> Mt, preserve sign
    df["value_mt"] = df["value"] / 1000.0

    tech_norm = df["techs"].astype(str).map(lambda t: alias.get(t, t))
    df["category"] = tech_norm.map(category_map)
    df = df[df["category"].notna()].copy()

    # aggregate by year, scenario, category
    agg = df.groupby(["year", "scenario", "category"], as_index=False, observed=False)["value_mt"].sum()

    # diagnostic: DAC values
    dac = agg[agg["category"] == "Direct air capture"].sort_values(["scenario", "year"])
    if not dac.empty:
        print("DAC CO2_emission used in plot (Mt, signed; negative = net removal):")
        print(dac.to_string(index=False))

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
    fig.subplots_adjust(left=0.10, right=0.78, top=0.97, bottom=0.22)

    # separate positive and negative contributors
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
    ax.set_xticklabels(scen_labels, rotation=90, va="top")

    sec = ax.secondary_xaxis("bottom")
    sec.set_xlim(ax.get_xlim())
    sec.set_xticks(year_centers)
    sec.set_xticklabels([str(y) for y in years], rotation=0)
    sec.spines["bottom"].set_visible(False)
    sec.tick_params(axis="x", pad=YEAR_PAD, length=0)

    handles = [plt.Line2D([0], [0], color=C.get(c, "#9E9E9E"), lw=8) for c in piv.columns]
    ax.legend(
        handles,
        list(piv.columns),
        title="Contributor",
        loc="center left",
        bbox_to_anchor=(0.93, 0.50),
        frameon=True,
        borderpad=0.9,
        labelspacing=0.55,
        handlelength=2.0,
        handletextpad=0.8,
        framealpha=1.0,
    )

    if EXPORT_FIGURES and outdir is not None:
        outdir.mkdir(parents=True, exist_ok=True)
        fig.savefig(outdir / "RQ3_co2_grouped.png", dpi=300, bbox_inches="tight")

    plt.show()

def inspect_indicator_detailed(ind_det_long: pd.DataFrame, max_per_field: int = 200) -> None:
    """
    Print unique values for key fields in indicator_accounting_detailed
    across all scenarios, to understand what to include/exclude in plots.
    """

    if ind_det_long.empty:
        print("indicator_accounting_detailed is empty; nothing to inspect.")
        return

    print("\n=== indicator_accounting_detailed: columns ===")
    print(list(ind_det_long.columns))

    def _print_unique(col: str):
        if col not in ind_det_long.columns:
            print(f"\n[{col}] column not present.")
            return
        vals = ind_det_long[col].dropna().unique()
        vals_sorted = np.sort(vals.astype(str))
        print(f"\n[{col}] unique values (n={len(vals_sorted)}):")
        if len(vals_sorted) <= max_per_field:
            for v in vals_sorted:
                print("  ", v)
        else:
            head = vals_sorted[:max_per_field]
            tail = vals_sorted[-5:]
            for v in head:
                print("  ", v)
            print(f"  ... ({len(vals_sorted) - max_per_field} more)")
            print("  (last 5):")
            for v in tail:
                print("  ", v)

    # indicators
    _print_unique("indicator")

    # nodesModel (if present)
    _print_unique("nodesModel")

    # years / accYears
    if "years" in ind_det_long.columns:
        _print_unique("years")
    if "accYears" in ind_det_long.columns:
        _print_unique("accYears")

    # techs
    _print_unique("techs")

    # quick cross-tab: indicators by tech prefix (first 20 rows)
    print("\n=== Sample indicator x tech summary (first 20 rows) ===")
    tech_str = ind_det_long["techs"].astype(str)
    tech_prefix = tech_str.str.extract(r"^([^_]+)", expand=False)
    ind_det_long["_tech_prefix"] = tech_prefix
    tab = (
        ind_det_long.groupby(["indicator", "_tech_prefix"], observed=False)["value"]
        .size()
        .reset_index(name="n")
        .sort_values(["indicator", "n"], ascending=[True, False])
    )
    print(tab.head(20).to_string(index=False))
    ind_det_long.drop(columns=["_tech_prefix"], inplace=True, errors="ignore")

def load_hourly_commodity_balance(case_dir: str) -> pd.DataFrame | None:
    res = load_results(case_dir)
    if res is None:
        return None
    try:
        df = res["commodity_balance"]  # adjust name if different
    except KeyError:
        print("No hourly commodity_balance in", case_dir)
        return None
    if isinstance(df, pd.Series):
        df = df.to_frame("value")
    if "value" not in df.columns:
        df = df.copy()
        df.columns = ["value"]
    return df.reset_index()


def make_day_hour_carpet(df: pd.DataFrame,
                         time_col: str = "time",
                         value_col: str = "value") -> np.ndarray:
    """
    df: rows for a single scenario/year, one value per hour (0..8759).
    Returns array of shape (365, 24): row = day (0..364), col = hour (0..23).
    """
    arr = df.sort_values(time_col)[value_col].to_numpy()
    if arr.size != 8760:
        raise ValueError(f"Expected 8760 hours, got {arr.size}")
    return arr.reshape(365, 24)

def setup_carpet_axes(ax: plt.Axes, title: str):
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("Day of year")
    ax.set_ylabel("Hour of day")
    ax.set_yticks([0, 6, 12, 18, 23])
    ax.set_yticklabels(["0", "6", "12", "18", "23"])

def build_hourly_profiles_2050() -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Returns:
      elec_op: hourly electrolyser operation (MW or similar) for year=2050
      h2_soc: hourly H2 storage SOC for year=2050
    Columns: scenario, time (0..8759), value.
    """
    elec_list = []
    soc_list = []

    for case_dir in CASE_DIRS:
        scen = scenario_name(case_dir)

        cba = load_hourly_commodity_balance(case_dir)
        soc = load_hourly_storage_soc(case_dir)

        if cba is None or soc is None:
            continue

        # assume columns include 'years' or 'accYears' and 'time'
        # filter to 2050
        year_col = "years" if "years" in cba.columns else "accYears"
        cba_2050 = cba[cba[year_col] == 2050].copy()

        # Electrolyser operation: commodity H2, balanceType "out" or "net" > 0
        # adjust filters to match your exact index set
        elec = cba_2050[
            (cba_2050["techs"].astype(str) == "Electrolyser")
            & (cba_2050["commodity"].astype(str) == "H2")
        ].copy()

        # aggregate across nodesModel etc. to total operation
        elec_grouped = (
            elec.groupby("time", observed=False)["value"]
            .sum()
            .reset_index()
            .sort_values("time")
        )
        elec_grouped["scenario"] = scen
        elec_list.append(elec_grouped[["scenario", "time", "value"]])

        # H2 storage SOC
        year_col_soc = "years" if "years" in soc.columns else "accYears"
        soc_2050 = soc[soc[year_col_soc] == 2050].copy()
        soc_h2 = soc_2050[soc_2050["techs"].astype(str) == "H2_storage"].copy()
        soc_grouped = (
            soc_h2.groupby("time", observed=False)["value"]
            .sum()
            .reset_index()
            .sort_values("time")
        )
        soc_grouped["scenario"] = scen
        soc_list.append(soc_grouped[["scenario", "time", "value"]])

    elec_op = pd.concat(elec_list, ignore_index=True) if elec_list else pd.DataFrame()
    h2_soc = pd.concat(soc_list, ignore_index=True) if soc_list else pd.DataFrame()
    return elec_op, h2_soc

def plot_electrolyser_and_h2soc_2050(elec_op: pd.DataFrame,
                                     h2_soc: pd.DataFrame,
                                     outdir: Path | None = None):
    sns.set_theme(style="white")
    infer = colormaps.get_cmap("inferno")

    fig, axes = plt.subplots(len(SCEN_ORDER), 2, figsize=(10, 10), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.07, right=0.90, top=0.95, bottom=0.08, hspace=0.25, wspace=0.12)

    for i, scen in enumerate(SCEN_ORDER):
        # Electrolyser
        df_e = elec_op[elec_op["scenario"] == scen].copy()
        if not df_e.empty:
            carpet_e = make_day_hour_carpet(df_e, time_col="time", value_col="value")
            im = axes[i, 0].imshow(
                carpet_e.T,  # transpose so x=day, y=hour
                origin="lower",
                aspect="auto",
                cmap=infer,
            )
            setup_carpet_axes(axes[i, 0], f"{scen} – Electrolyser")
        else:
            axes[i, 0].set_visible(False)

        # H2 SOC
        df_s = h2_soc[h2_soc["scenario"] == scen].copy()
        if not df_s.empty:
            carpet_s = make_day_hour_carpet(df_s, time_col="time", value_col="value")
            axes[i, 1].imshow(
                carpet_s.T,
                origin="lower",
                aspect="auto",
                cmap=infer,
            )
            setup_carpet_axes(axes[i, 1], f"{scen} – H$_2$ storage SoC")
        else:
            axes[i, 1].set_visible(False)

    cbar_ax = fig.add_axes([0.92, 0.10, 0.015, 0.80])
    fig.colorbar(im, cax=cbar_ax, label="Value (model units)")

    if EXPORT_FIGURES and outdir is not None:
        outdir.mkdir(parents=True, exist_ok=True)
        fig.savefig(outdir / "RQx_electrolyser_h2soc_2050.png", dpi=300, bbox_inches="tight")

    plt.show()

def build_hourly_dac_and_pv_2050() -> tuple[pd.DataFrame, pd.DataFrame]:
    dac_list = []
    pv_list = []

    for case_dir in CASE_DIRS:
        scen = scenario_name(case_dir)
        cba = load_hourly_commodity_balance(case_dir)
        if cba is None:
            continue

        year_col = "years" if "years" in cba.columns else "accYears"
        cba_2050 = cba[cba[year_col] == 2050].copy()

        # DAC operation: CO2 commodity, DAC tech
        dac = cba_2050[cba_2050["techs"].astype(str) == "DAC"].copy()
        dac_grouped = (
            dac.groupby("time", observed=False)["value"]
            .sum()
            .reset_index()
            .sort_values("time")
        )
        dac_grouped["scenario"] = scen
        dac_list.append(dac_grouped[["scenario", "time", "value"]])

        # Solar PV: aggregate pv_central_fixed + pv_decentral, Elec commodity
        pv = cba_2050[
            cba_2050["techs"].astype(str).isin(["pv_central_fixed", "pv_decentral"])
            & (cba_2050["commodity"].astype(str) == "Elec")
        ].copy()
        pv_grouped = (
            pv.groupby("time", observed=False)["value"]
            .sum()
            .reset_index()
            .sort_values("time")
        )
        pv_grouped["scenario"] = scen
        pv_list.append(pv_grouped[["scenario", "time", "value"]])

    dac_op = pd.concat(dac_list, ignore_index=True) if dac_list else pd.DataFrame()
    pv_op = pd.concat(pv_list, ignore_index=True) if pv_list else pd.DataFrame()
    return dac_op, pv_op


def plot_dac_and_pv_2050(dac_op: pd.DataFrame,
                         pv_op: pd.DataFrame,
                         outdir: Path | None = None):
    sns.set_theme(style="white")
    infer = colormaps.get_cmap("inferno")

    fig, axes = plt.subplots(len(SCEN_ORDER), 2, figsize=(10, 10), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.07, right=0.90, top=0.95, bottom=0.08, hspace=0.25, wspace=0.12)

    for i, scen in enumerate(SCEN_ORDER):
        df_d = dac_op[dac_op["scenario"] == scen].copy()
        if not df_d.empty:
            carpet_d = make_day_hour_carpet(df_d, time_col="time", value_col="value")
            im = axes[i, 0].imshow(
                carpet_d.T,
                origin="lower",
                aspect="auto",
                cmap=infer,
            )
            setup_carpet_axes(axes[i, 0], f"{scen} – DAC")
        else:
            axes[i, 0].set_visible(False)

        df_p = pv_op[pv_op["scenario"] == scen].copy()
        if not df_p.empty:
            carpet_p = make_day_hour_carpet(df_p, time_col="time", value_col="value")
            axes[i, 1].imshow(
                carpet_p.T,
                origin="lower",
                aspect="auto",
                cmap=infer,
            )
            setup_carpet_axes(axes[i, 1], f"{scen} – Solar PV")
        else:
            axes[i, 1].set_visible(False)

    cbar_ax = fig.add_axes([0.92, 0.10, 0.015, 0.80])
    fig.colorbar(im, cax=cbar_ax, label="Value (model units)")

    if EXPORT_FIGURES and outdir is not None:
        outdir.mkdir(parents=True, exist_ok=True)
        fig.savefig(outdir / "RQx_dac_pv_2050.png", dpi=300, bbox_inches="tight")

    plt.show()


def load_hourly_storage_soc(case_dir):
    try:
        return res["storage_soc"]
    except Exception as e:
        print("WARNING: no storage_soc in GDX, skipping SOC.", e)
        return None

def load_hourly_commodity_balance(case_dir: str) -> pd.DataFrame | None:
    """
    Load full-resolution commodity_balance with timeModel (tm1..tm8760) for a case.
    Returns a long df with columns incl. timeModel, years/accYears, techs, commodity, value, etc.
    """
    res = load_results(case_dir)
    if res is None:
        return None
    try:
        df = res["commodity_balance"]  # hourly symbol, NOT the annual one
    except KeyError:
        print("No hourly commodity_balance in", case_dir)
        return None
    if isinstance(df, pd.Series):
        df = df.to_frame("value")
    if "value" not in df.columns:
        df = df.copy()
        df.columns = ["value"]
    return df.reset_index()


FINAL_DEMAND_TECHS = {
    "Elec_demand",
    "HeatDemand_Bio_LF",
    "HeatDemand_REfuel",
    "HeatDemand_e-CH4",
    "TranspDemand_REfuel",
    "TranspDemand_H2",
    "TranspDemand_e-CH4",
}

def rq2_tech_group(tech: str) -> str | None:
    t = str(tech)

    if t.startswith("FuelImport_"):
        return None
    if t in FINAL_DEMAND_TECHS:
        return None

    # generation, electricity
    if t.startswith("Thermal_") or t in {"GT", "OCGT", "CCGT"}:
        return "Thermal and gas turbines"
    if t.startswith("wind_"):
        return "Wind"
    if t.startswith("pv_"):
        return "Solar PV"
    if t == "H2_FC":
        return "H2 fuel cell"

    # storage
    if t == "Battery":
        return "Battery"
    if t == "H2_storage":
        return "H2 storage"

    # fuel/chemical conversion, counted on their Elec input
    if t in {"Electrolyser", "FTropschSyn", "Methanizer", "DAC"}:
        return t  # keep name as-is

    # ignore others (Hydro, Geothermal, etc.) for these plots
    return None

def parse_timeModel_to_hour(time_label: str) -> int:
    # "tm1" -> 0, "tm8760" -> 8759
    s = str(time_label)
    if not s.startswith("tm"):
        raise ValueError(f"Unexpected timeModel label {s}")
    return int(s[2:]) - 1


def make_day_hour_carpet_from_tm(df: pd.DataFrame,
                                 time_col: str = "timeModel",
                                 value_col: str = "value") -> np.ndarray:
    """
    df: single scenario, single year, 8760 hourly values.
    Returns (365, 24) array (day, hour) for carpet plots.
    """
    tmp = df.copy()
    tmp["hour_idx"] = tmp[time_col].map(parse_timeModel_to_hour)
    tmp = tmp.sort_values("hour_idx")
    arr = tmp[value_col].to_numpy()
    if arr.size != 8760:
        raise ValueError(f"Expected 8760 hours for carpet, got {arr.size}")
    return arr.reshape(365, 24)

def build_rq2_hourly_for_year(year: int,
                              left_group: str,
                              right_group: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build hourly series for selected tech groups and year:
    Returns two dfs:
      left_df: scenario, timeModel, value  (for left_group)
      right_df: scenario, timeModel, value (for right_group)
    Includes:
      - row for GP 2020 baseline (scenario label 'GP', year=2020) -> used as "(all)"
      - rows for all scenarios in year 'year'.
    For conversion techs, only Elec commodity is kept; for storage:
      Battery -> Elec, H2_storage -> H2.
    """
    left_list, right_list = [], []

    for case_dir in CASE_DIRS:
        scen = scenario_name(case_dir)
        cba = load_hourly_commodity_balance(case_dir)
        if cba is None:
            continue

        # identify year column
        year_col = "years" if "years" in cba.columns else "accYears"
        time_col = "timeModel"  # your index

        # we need 2020 GP baseline + year for all scenarios
        for target_year, scen_label in [(2020, "GP"), (year, scen)]:
            if target_year not in YEAR_INTS:
                continue
            df_y = cba[cba[year_col] == target_year].copy()
            if df_y.empty:
                continue

            # tech grouping
            df_y["tech_group"] = df_y["techs"].astype(str).map(rq2_tech_group)
            df_y = df_y[df_y["tech_group"].notna()].copy()

            # commodity filters
            tech_g = df_y["tech_group"].astype(str)
            comm = df_y["commodity"].astype(str)

            # for storage: enforce right commodity
            is_batt = tech_g == "Battery"
            is_h2stor = tech_g == "H2 storage"

            df_y = df_y[
                (~is_batt | (comm == "Elec"))
                & (~is_h2stor | (comm == "H2"))
            ].copy()

            # for conversion (Electrolyser, FTropschSyn, Methanizer, DAC), keep only Elec
            is_conv = tech_g.isin(["Electrolyser", "FTropschSyn", "Methanizer", "DAC"])
            df_y = df_y[
                (~is_conv) | (comm == "Elec")
            ].copy()

            # aggregate across nodes and techs for each group
            agg = (
                df_y.groupby([time_col, "tech_group"], observed=False)["value"]
                .sum()
                .reset_index()
            )
            agg["scenario"] = scen_label

            # select left and right groups
            left_g = agg[agg["tech_group"] == left_group][["scenario", time_col, "value"]].copy()
            if not left_g.empty:
                left_list.append(left_g)

            right_g = agg[agg["tech_group"] == right_group][["scenario", time_col, "value"]].copy()
            if not right_g.empty:
                right_list.append(right_g)

    left_df = pd.concat(left_list, ignore_index=True) if left_list else pd.DataFrame()
    right_df = pd.concat(right_list, ignore_index=True) if right_list else pd.DataFrame()

    return left_df, right_df


def plot_rq2_carpet(year: int,
                    left_group: str,
                    right_group: str,
                    outdir: Path | None = None):
    """
    Make a 6×2 carpet figure for a given year:
      row 0: 2020 GP baseline (labelled "(all)")
      rows 1–5: GP, NT, ELEC+, BIO+, H2+ in 'year'
    left column: left_group, right column: right_group.
    """
    left_df, right_df = build_rq2_hourly_for_year(year, left_group, right_group)

    if left_df.empty and right_df.empty:
        print(f"RQ2: no data for groups {left_group} / {right_group} in {year}.")
        return

    sns.set_theme(style="white")
    infer = colormaps.get_cmap("inferno")

    # scenario sequence for rows
    row_labels = ["(all) 2020 GP"] + SCEN_ORDER  # 6 rows total
    nrows = len(row_labels)

    fig, axes = plt.subplots(nrows, 2, figsize=(10, 11), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.07, right=0.90, top=0.95, bottom=0.07, hspace=0.30, wspace=0.12)

    def select_rows(df: pd.DataFrame, scen_label: str, target_year: int, baseline: bool) -> pd.DataFrame:
        if df.empty:
            return df
        if baseline:
            # baseline: scenario == "GP" and year 2020 (we already labelled scenario as "GP" in builder)
            return df[df["scenario"] == "GP"].copy()
        else:
            return df[df["scenario"] == scen_label].copy()

    last_im = None

    for i, label in enumerate(row_labels):
        baseline = (i == 0)
        scen_label = "GP" if baseline else SCEN_ORDER[i - 1]

        # left column
        df_l = select_rows(left_df, scen_label, year, baseline)
        if not df_l.empty:
            carpet_l = make_day_hour_carpet_from_tm(df_l, time_col="timeModel", value_col="value")
            last_im = axes[i, 0].imshow(
                carpet_l.T,
                origin="lower",
                aspect="auto",
                cmap=infer,
            )
            title_l = f"{left_group} – {label}"
            axes[i, 0].set_title(title_l, fontsize=9)
            axes[i, 0].set_ylabel("Hour of day")
            axes[i, 0].set_yticks([0, 6, 12, 18, 23])
            axes[i, 0].set_yticklabels(["0", "6", "12", "18", "23"])
        else:
            axes[i, 0].set_visible(False)

        # right column
        df_r = select_rows(right_df, scen_label, year, baseline)
        if not df_r.empty:
            carpet_r = make_day_hour_carpet_from_tm(df_r, time_col="timeModel", value_col="value")
            last_im = axes[i, 1].imshow(
                carpet_r.T,
                origin="lower",
                aspect="auto",
                cmap=infer,
            )
            title_r = f"{right_group} – {label}"
            axes[i, 1].set_title(title_r, fontsize=9)
        else:
            axes[i, 1].set_visible(False)

        # x-labels only on bottom row
        if i == nrows - 1:
            axes[i, 0].set_xlabel("Day of year")
            axes[i, 1].set_xlabel("Day of year")

    # shared colourbar
    if last_im is not None:
        cbar_ax = fig.add_axes([0.92, 0.10, 0.015, 0.80])
        fig.colorbar(last_im, cax=cbar_ax, label="Value (model units)")

    if EXPORT_FIGURES and outdir is not None:
        outdir.mkdir(parents=True, exist_ok=True)
        fname = f"RQ2_{left_group.replace(' ','_')}_{right_group.replace(' ','_')}_{year}.png"
        fig.savefig(outdir / fname, dpi=300, bbox_inches="tight")

    plt.show()



def main() -> None:
    outdir = BASE_DIR / "comparison_outputs"

    cba_long, caps_long, ind_long, ind_det_long = load_needed_data()

    plot_capacity_and_generation_grouped(caps_long=caps_long, cba_long=cba_long, outdir=outdir)
    plot_fuelconv_and_supply(caps_long=caps_long, cba_long=cba_long, outdir=outdir)
    plot_costs_and_lcoe(ind_long=ind_long, ind_det_long=ind_det_long, cba_long=cba_long, outdir=outdir)
    plot_co2_grouped(ind_det_long=ind_det_long, outdir=outdir)


    outdir = BASE_DIR / "comparison_outputs"

    # Example: Electrolyser vs H2 storage in 2050
    plot_rq2_carpet(
        year=2050,
        left_group="Electrolyser",
        right_group="H2 storage",
        outdir=outdir,
    )

    # Example: DAC vs Solar PV in 2050
    plot_rq2_carpet(
        year=2050,
        left_group="DAC",
        right_group="Solar PV",
        outdir=outdir,
    )


if __name__ == "__main__":
    main()
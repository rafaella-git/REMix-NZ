from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib import colormaps
from matplotlib import cm
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

GEN_TECH_MIN_TWH = 0.01
CAP_TECH_MIN_GW = 0.01
BAR_WIDTH = 0.85

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

CARRIER_LABEL_MAP = {
    "Coal": "Solid (coal)",
    "Diesel": "Liquid fuel (fossil)",
    "Petrol": "Liquid fuel (fossil)",
    "Fossil liquid fuel": "Liquid fuel (fossil)",
    "Natural gas": "Methane (fossil)",
    "Fossil gas": "Methane (fossil)",
    "Bio liquid fuel": "Liquid fuel (bio)",
    "Biogas": "Methane (bio)",
    "Wood": "Solid (wood)",
    "E-fuel liquid": "Liquid fuel (e-FTL)",
    "E-fuel gas": "Methane (e-CH4)",
    "Geothermal": "Thermal energy (geothermal)",
    "Solar": "Thermal energy (solar)",
    "Hydrogen": "Hydrogen",
    "Electricity": "Electricity",
}


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


def build_grouped_x(years: list[int], scen_order: list[str], dx_in: float = 1.0, year_gap: float = 2.0):
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
    s = s.replace("CH_4", r"$\mathrm{CH_4}$")
    s = s.replace("CO_2", r"$\mathrm{CO_2}$")
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

    # If nodesModel exists and has "global", keep only that; otherwise sum all nodes
    if "nodesModel" in df.columns:
        if (df["nodesModel"] == "global").any():
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

    if "accNodesModel" in df.columns:
        if (df["accNodesModel"] == "global").any():
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


def build_fuel_supply(cba_long: pd.DataFrame) -> pd.DataFrame:
    df = cba_long.copy()

    if "accNodesModel" in df.columns:
        if (df["accNodesModel"] == "global").any():
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
            return "E-fuel gas"
        if comm == "REfuel":
            return "E-fuel liquid"
        if comm == "H2":
            return "Hydrogen"

        if t in ["FuelImport_Bio_LF", "FuelImport_Biofuel"]:
            return "Bio liquid fuel"
        if t == "FuelImport_Bio_gas":
            return "Biogas"
        if t in ["FuelImport_CH4", "FuelImport_Fossil_CH4"]:
            return "Natural gas"
        if t == "FuelImport_Coal":
            return "Coal"
        if t in ["FuelImport_Diesel", "FuelImport_Fossil_LF"]:
            return "Diesel"

        return None

    df["carrier_raw"] = df.apply(supply_group, axis=1)
    df = df[df["carrier_raw"].notna()].copy()
    df["carrier"] = df["carrier_raw"].map(lambda x: CARRIER_LABEL_MAP.get(x, x))
    df["value_twh"] = df["value"] / 1000.0
    df = df[df["value_twh"].abs() > 0.01].copy()

    out = df.groupby(["year", "scenario", "carrier"], as_index=False, observed=False)["value_twh"].sum()
    return out


def build_fuelconv_caps(caps_long: pd.DataFrame) -> pd.DataFrame:
    df = caps_long.copy()

    if "accNodesModel" in df.columns:
        if (df["accNodesModel"] == "global").any():
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


def build_lcoe_components(ind_det_long: pd.DataFrame, cba_long: pd.DataFrame) -> pd.DataFrame:
    if ind_det_long.empty:
        return pd.DataFrame(columns=["year", "scenario", "component", "eur_per_mwh"])

    df = ind_det_long.copy()
    df = df[(df["year"].isin(YEAR_INTS)) & (df["indicator"].isin(COST_COMPONENTS))].copy()

    tech = df["techs"].astype(str)
    df = df[~tech.str.startswith(("Transp", "Heat", "FuelImport"), na=False)].copy()

    is_fuelcost = df["indicator"].astype(str) == "FuelCost"
    is_import = df["techs"].astype(str).str.startswith("FuelImport_", na=False)
    df = df[~(is_fuelcost & is_import)].copy()

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

    df["cost_eur"] = df["value"] * 1e6
    cost_agg = df.groupby(["year", "scenario", "component"], as_index=False, observed=False)["cost_eur"].sum()

    gen = cba_long[
        (cba_long["accNodesModel"] == "global")
        & (cba_long["commodity"] == "Elec")
        & (cba_long["balanceType"] == "net")
        & (cba_long["year"].isin(YEAR_INTS))
        & (cba_long["value"] > 0)
    ].copy()
    gen = gen[~gen["techs"].astype(str).str.contains("Battery", case=False, regex=False)].copy()
    gen_gwh = gen.groupby(["year", "scenario"], as_index=False, observed=False)["value"].sum()
    gen_gwh = gen_gwh.rename(columns={"value": "gen_gwh"})
    gen_gwh["gen_mwh"] = gen_gwh["gen_gwh"] * 1000.0

    merged = cost_agg.merge(gen_gwh[["year", "scenario", "gen_mwh"]], on=["year", "scenario"], how="left")
    merged["eur_per_mwh"] = (merged["cost_eur"] / merged["gen_mwh"]).replace(
        [np.inf, -np.inf], np.nan
    )
    out = merged.groupby(["year", "scenario", "component"], as_index=False, observed=False)["eur_per_mwh"].sum()
    out = out[out["eur_per_mwh"].abs() > 0.01].copy()
    return out


def evenly_spaced_colors_from_cmap(cmap_name, n, start=0.0, stop=1.0):
    cmap = cm.get_cmap(cmap_name)
    return [cmap(x) for x in np.linspace(start, stop, n)]


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

    fig, axes = plt.subplots(2, 1, figsize=(6.5, 7.0), sharex=True)
    fig.subplots_adjust(left=0.10, right=0.74, top=0.97, bottom=0.13, hspace=0.05)


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

    shared_legend(fig, stack_order, cmap, title="Technology", x=0.78, y=0.5)

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

    unique_carriers = sorted(sup["carrier"].unique().tolist())
    cmap_inferno = colormaps.get_cmap("inferno")

    def carrier_group(name: str) -> int:
        n = name.lower()
        if any(k in n for k in ["coal", "diesel", "petrol", "fuel", "gas", "lng"]) and "bio" not in n and "wood" not in n:
            return 0
        if "bio" in n or "wood" in n:
            return 1
        if "geothermal" in n:
            return 2
        if "solar" in n:
            return 3
        if any(k in n for k in ["electric", "hydrogen", "e-fuel", "efuel"]):
            return 4
        return 1

    group_to_frac = {
        0: 0.10,
        1: 0.40,
        2: 0.55,
        3: 0.70,
        4: 0.90,
    }

    carrier_groups = {c: carrier_group(c) for c in unique_carriers}
    carrier_colors = {c: cmap_inferno(group_to_frac[carrier_groups[c]]) for c in unique_carriers}

    cap_cmap = {
        "Methanizer": carrier_colors.get("Methane (e-CH4)", carrier_colors.get("Methane (fossil)", DEFAULT_COLOR)),
        "FTropschSyn": carrier_colors.get("Liquid fuel (e-FTL)", carrier_colors.get("Liquid fuel (fossil)", DEFAULT_COLOR)),
        "Electrolyser": carrier_colors.get("Hydrogen", DEFAULT_COLOR),
        "DAC": "#8B5A2B",
    }

    fig, axes = plt.subplots(2, 1, figsize=(6.5, 7.0), sharex=True)
    fig.subplots_adjust(left=0.10, right=0.74, top=0.97, bottom=0.13, hspace=0.05)

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
        category_col="carrier",
        y_label="Supply of fuels and chemicals (TWh)",
        years=YEAR_INTS,
        scen_order=SCEN_ORDER,
        cmap=carrier_colors,
        cat_order=unique_carriers,
        dx_in=DX_IN,
        year_gap=year_gap,
        year_pad=YEAR_PAD,
        show_scen_labels=True,
    )

    LEGEND_X = 0.78
    LEGEND_Y_TOP = 0.78
    LEGEND_Y_BOT = 0.22

    handles0 = [plt.Line2D([0], [0], color=cap_cmap.get(k, DEFAULT_COLOR), lw=10) for k in cap_order]
    leg0 = fig.legend(
        handles0,
        cap_order,
        title="Technology",
        loc="center left",
        bbox_to_anchor=(LEGEND_X, LEGEND_Y_TOP),
        frameon=True,
        borderpad=1.2,
        labelspacing=0.8,
        handlelength=1.8,
        handletextpad=0.8,
    )
    leg0.get_frame().set_linewidth(1.2)

    carrierlabels = [legend_math_subscripts(x).replace(" for ", "\nfor ") for x in unique_carriers]
    handles1 = [
        plt.Line2D([0], [0], color=carrier_colors.get(k, DEFAULT_COLOR), lw=10)
        for k in unique_carriers
    ]
    leg1 = fig.legend(
        handles1,
        carrierlabels,
        title="Carrier",
        loc="center left",
        bbox_to_anchor=(LEGEND_X, LEGEND_Y_BOT),
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


def plot_costs_and_lcoe(ind_long: pd.DataFrame, ind_det_long: pd.DataFrame, cba_long: pd.DataFrame, outdir: Path | None = None):
    sns.set_theme(style="whitegrid")

    costs = build_cost_long_from_detailed(ind_det_long)
    lcoe_comp = build_lcoe_components(ind_det_long, cba_long)

    if costs.empty:
        print("No indicator_accounting loaded; skipping cost/LCOE plots.")
        return

    comp_order = ["Fuel cost", "Investment", "OMFix", "Slack cost", "Spill penalty", r"CO$_{2}$ slack"]
    comp_order = [c for c in comp_order if c in costs["component"].unique().tolist()]

    colors = evenly_spaced_colors_from_cmap("viridis", n=len(comp_order), start=0.0, stop=1.0)
    cmap = {c: col for c, col in zip(comp_order, colors)}

    fig, axes = plt.subplots(2, 1, figsize=(6.5, 7.0), sharex=True)
    fig.subplots_adjust(left=0.10, right=0.74, top=0.97, bottom=0.13, hspace=0.05)


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

    stacked_grouped_bars(
        ax=axes[1],
        df_long=lcoe_comp,
        value_col="eur_per_mwh",
        category_col="component",
        y_label="LCOE decomposition (EUR/MWh)",
        years=YEAR_INTS,
        scen_order=SCEN_ORDER,
        cmap=cmap,
        cat_order=comp_order,
        dx_in=DX_IN,
        year_gap=year_gap,
        year_pad=YEAR_PAD,
        show_scen_labels=True,
    )

    leg = shared_legend(fig, comp_order, cmap, title="Cost component", x=0.78, y=0.5)

    if EXPORT_FIGURES and outdir is not None:
        outdir.mkdir(parents=True, exist_ok=True)
        fig.savefig(outdir / "RQ3_costs_and_lcoe.png", dpi=300, bbox_inches="tight")

    plt.show()


def plot_co2_grouped(ind_det_long: pd.DataFrame, outdir: Path | None = None):
    sns.set_theme(style="whitegrid")

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

    fig, ax = plt.subplots(1, 1, figsize=(9, 4.65))
    fig.subplots_adjust(left=0.08, right=0.78, top=0.98, bottom=0.20)

    bottom = np.zeros(len(x))
    for c in piv.columns:
        vals = piv[c].values
        ax.bar(x, vals, bottom=bottom, width=BAR_WIDTH, color=C.get(c, "#9E9E9E"), edgecolor="none")
        bottom += vals

    ax.grid(axis="y", alpha=0.25)
    ax.set_ylabel("Total CO2 emissions by contributor (Mt CO2)\n(negative = net removal)")

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
        bbox_to_anchor=(0.95, 0.50),
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


def main() -> None:
    outdir = BASE_DIR / "comparison_outputs"

    cba_long, caps_long, ind_long, ind_det_long = load_needed_data()

    plot_capacity_and_generation_grouped(caps_long=caps_long, cba_long=cba_long, outdir=outdir)
    plot_fuelconv_and_supply(caps_long=caps_long, cba_long=cba_long, outdir=outdir)
    plot_costs_and_lcoe(ind_long=ind_long, ind_det_long=ind_det_long, cba_long=cba_long, outdir=outdir)
    plot_co2_grouped(ind_det_long=ind_det_long, outdir=outdir)


if __name__ == "__main__":
    main()

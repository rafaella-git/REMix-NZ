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
CAP_TECH_MIN_GW = 0.01  # 10 MW
BAR_WIDTH = 0.85

FUEL_CONV_TECHS = ["Methanizer", "FTropschSyn", "Electrolyser", "DAC"]
FUEL_CONV_STACK_ORDER = ["Electrolyser", "Methanizer", "FTropschSyn", "DAC"]

COST_COMPONENTS = ["FuelCost", "Invest", "OMFix", "SlackCost",  "SpillPenalty", "Slack_CO2"] #"SystemCost",


LEGEND_CO2_ANCHOR = (0.83, 0.50)
LEGEND_BORDER_PAD = 0.9
LEGEND_LABEL_SPACING = 0.55
LEGEND_FRAME_ALPHA = 1.0

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
        raise KeyError(f"Could not find a year column. Tried {year_col_candidates}. Columns: {list(df.columns)}")

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

def legend_label_mathtext(label: str) -> str:
    # Example: "Methane (e-fuel) for ..." -> "Methane (e-fuel)\nfor ..."
    label = label.replace(" for ", "\nfor ")
    # Enable subscripts where you use "_" in labels
    # e.g., "CH_4" becomes subscript 4 in mathtext
    if "_" in label:
        label = r"$\mathrm{" + label.replace("\n", r"}$" + "\n" + r"$\mathrm{").replace("_", r"_{") + ("}" * label.count("_")) + r"}$"
        # If this gets messy with multiple underscores, tell me what exact strings you want subscripted.
    return label

def legend_math_subscripts(s: str) -> str:
    # Only replace exact tokens; preserve all other "_" and "-" unchanged.
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
        # year_pad=year_pad,
        show_scen=show_scen_labels,
    )

    ymax = float((piv.sum(axis=1)).max())
    ax.set_ylim(0.0, ymax * 1.05 if ymax > 0 else 1.0)


    return piv.columns.tolist()


def shared_legend(fig, labels: list[str], cmap: dict[str, str], x: float, y: float):
    handles = [plt.Line2D([0], [0], color=cmap.get(l, DEFAULT_COLOR), lw=10) for l in labels]
    leg = fig.legend(
        handles,
        labels,
        title="Technology",
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
    df = caps_long[
        (caps_long["nodesModel"] == "global")
        & (caps_long["capType"] == "total")
        & (caps_long["commodity"] == "Elec")
        & (caps_long["year"].isin(YEAR_INTS))
        & (caps_long["value"] > 0)
    ].copy()

    df["tech"] = df["techs"].astype(str)
    df = df[df["tech"] != "Battery"].copy()
    df = df[df["value"] >= CAP_TECH_MIN_GW].copy()

    df["tech"] = df["tech"].map(tech_group_pretty)
    out = df.groupby(["year", "scenario", "tech"], as_index=False, observed=False)["value"].sum()
    return out


def build_generation_grouped(cba_long: pd.DataFrame) -> pd.DataFrame:
    df = cba_long[
        (cba_long["accNodesModel"] == "global")
        & (cba_long["commodity"] == "Elec")
        & (cba_long["balanceType"] == "net")
        & (cba_long["year"].isin(YEAR_INTS))
        & (cba_long["value"] > 0)
    ].copy()

    df["tech"] = df["techs"].astype(str)
    df["tech"] = df["tech"].map(tech_group_pretty)
    df["value_twh"] = df["value"] / 1000.0

    totals = df.groupby("tech", observed=False)["value_twh"].sum()
    keep = totals[totals >= GEN_TECH_MIN_TWH].index.tolist()
    df = df[df["tech"].isin(keep)].copy()

    out = df.groupby(["year", "scenario", "tech"], as_index=False, observed=False)["value_twh"].sum()
    return out


def plot_capacity_and_generation_grouped(caps_long: pd.DataFrame, cba_long: pd.DataFrame):
    sns.set_theme(style="whitegrid")

    cap = build_capacity_grouped(caps_long)
    gen = build_generation_grouped(cba_long)

    techs = sorted(set(cap["tech"].unique()) | set(gen["tech"].unique()))
    cmap = color_map_for(techs)

    stack_order = [t for t in ["Hydropower", "Geothermal", "Onshore wind", "Solar PV"] if t in techs] + [t for t in techs if t not in {"Hydropower", "Geothermal", "Onshore wind", "Solar PV"}]

    fig, axes = plt.subplots(2, 1, figsize=(18, 11), sharex=True)
    fig.subplots_adjust(left=0.08, right=0.75, top=0.98, bottom=0.12, hspace=0.06)

    year_gap = YEAR_GAP * 0.5

    stacked_grouped_bars(
        ax=axes[0],
        df_long=cap,
        value_col="value",
        category_col="tech",
        y_label="Installed electricity \ngeneration capacity (GW)",
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

    shared_legend(fig, stack_order, cmap, x=0.78, y=0.5)
    plt.show()


def make_inferno_palette_carriers(carriers: list[str]):
    infer = colormaps.get_cmap("inferno")
    ordered = [
        "Solid (coal)",
        "Liquid fuel (fossil)",
        "Methane (fossil)",
        "Liquid fuel (bio)",
        "Methane (bio)",
        "Solid (wood)",
        "Liquid fuel (e-fuel)",
        "Methane (e-fuel)",
        "Thermal energy (geothermal)",
        "Hydrogen",
        "Thermal energy (solar)",
        "Electricity",
    ]
    ordered = [c for c in ordered if c in carriers]

    n = max(1, len(ordered))
    base_colors = [infer(i / max(1, (n - 1))) for i in range(n)]
    palette = {c: base_colors[i] for i, c in enumerate(ordered)}

    mid_col = infer(0.5)
    for c in carriers:
        if c not in palette:
            palette[c] = mid_col

    return ordered, palette


def build_fuel_supply(cba_long: pd.DataFrame) -> pd.DataFrame:
    df = cba_long[
        (cba_long["year"].isin(YEAR_INTS))
        & (cba_long["balanceType"] == "net")
        & (cba_long["value"] > 0)
    ].copy()

    df = df[~df["techs"].astype(str).str.startswith("SlackFuel_")].copy()

    def supply_group(row) -> str | None:
        t = str(row["techs"])
        comm = str(row["commodity"])

        if comm == "e-CH4":
            return "Methane (e-fuel)"
        if comm == "REfuel":
            return "Liquid fuel (e-fuel)"
        if comm == "H2":
            return "Hydrogen"

        if t in ["FuelImport_Bio_LF", "FuelImport_Biofuel"]:
            return "Liquid fuel (bio)"
        if t == "FuelImport_Bio_gas":
            return "Methane (bio)"
        if t in ["FuelImport_CH4", "FuelImport_Fossil_CH4"]:
            return "Methane (fossil)"
        if t == "FuelImport_Coal":
            return "Solid (coal)"
        if t in ["FuelImport_Diesel", "FuelImport_Fossil_LF"]:
            return "Liquid fuel (fossil)"

        return None

    df["carrier"] = df.apply(supply_group, axis=1)
    df = df[df["carrier"].notna()].copy()
    df["value_twh"] = df["value"] / 1000.0

    out = df.groupby(["year", "scenario", "carrier"], as_index=False, observed=False)["value_twh"].sum()
    return out


def build_fuelconv_caps(caps_long: pd.DataFrame) -> pd.DataFrame:
    df = caps_long[
        (caps_long["accNodesModel"] == "global")
        & (caps_long["capType"] == "total")
        & (caps_long["year"].isin(YEAR_INTS))
        & (caps_long["value"] > 0)
        & (caps_long["techs"].astype(str).isin(FUEL_CONV_TECHS))
    ].copy()

    df["tech"] = df["techs"].astype(str)
    return df[["year", "scenario", "tech", "value"]]


def plot_fuelconv_and_supply(caps_long: pd.DataFrame, cba_long: pd.DataFrame):
    sns.set_theme(style="whitegrid")

    cap = build_fuelconv_caps(caps_long)
    sup = build_fuel_supply(cba_long)

    cap_order = [t for t in FUEL_CONV_STACK_ORDER if t in cap["tech"].unique().tolist()]


    carrierpresent = sorted(sup["carrier"].unique().tolist())
    carrierorder, carriercmap = make_inferno_palette_carriers(carrierpresent)

    cap_cmap = {
        "Methanizer": carriercmap.get("Methane (e-fuel)", DEFAULT_COLOR),
        "FTropschSyn": carriercmap.get("Liquid fuel (e-fuel)", DEFAULT_COLOR),
        "Electrolyser": carriercmap.get("Hydrogen", DEFAULT_COLOR),
        "DAC": "#8B5A2B",  # brown
    }


    carrierlabels = [legend_math_subscripts(x).replace(" for ", "\nfor ") for x in carrierorder]


    carrier_present = sorted(sup["carrier"].unique().tolist())
    carrier_order, carrier_cmap = make_inferno_palette_carriers(carrier_present)

    LEGEND_X = 0.78
    LEGEND_Y_TOP = 0.78
    LEGEND_Y_BOT = 0.22

    fig, axes = plt.subplots(2, 1, figsize=(18, 11), sharex=True)
    fig.subplots_adjust(left=0.08, right=0.75, top=0.98, bottom=0.12, hspace=0.06)

    year_gap = YEAR_GAP * 0.5

    stacked_grouped_bars(
        ax=axes[0],
        df_long=cap,
        value_col="value",
        category_col="tech",
        y_label = "Total installed capacity\nfor fuel conversion " + r"$\mathrm{GW_{output}}$",
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
        cmap=carrier_cmap,
        cat_order=carrier_order,
        dx_in=DX_IN,
        year_gap=year_gap,
        year_pad=YEAR_PAD,
        show_scen_labels=True,
    )

    handles0 = [plt.Line2D([0], [0], color=cap_cmap.get(k, DEFAULT_COLOR), lw=10) for k in cap_order]
    leg0 = fig.legend(
        handles0, cap_order, title="Technology",
        loc="center left", bbox_to_anchor=(LEGEND_X, LEGEND_Y_TOP),
        frameon=True, borderpad=1.2, labelspacing=0.8, handlelength=1.8, handletextpad=0.8
    )
    leg0.get_frame().set_linewidth(1.2)

    handles1 = [plt.Line2D([0], [0], color=carrier_cmap.get(k, DEFAULT_COLOR), lw=10) for k in carrier_order]

    carrierlabels = [legend_math_subscripts(x).replace(" for ", "\nfor ") for x in carrier_order]
    leg1 = fig.legend(
        handles1, carrierlabels, title="Carrier",
        loc="center left", bbox_to_anchor=(LEGEND_X, LEGEND_Y_BOT),
        frameon=True, borderpad=1.2, labelspacing=0.8, handlelength=1.8, handletextpad=0.8
    )
    leg1.get_frame().set_linewidth(1.2)

    plt.show()


def build_cost_long_from_detailed(ind_det_long: pd.DataFrame) -> pd.DataFrame:
    """
    Build cost components time series from indicator_accounting_detailed,
    removing double-counted FuelCost on FuelImport_* techs.
    Output: year, scenario, component, value_bEUR
    """
    if ind_det_long.empty:
        return pd.DataFrame(columns=["year", "scenario", "component", "value_bEUR"])

    df = ind_det_long.copy()

    # keep only desired components
    df = df[
        (df["year"].isin(YEAR_INTS))
        & (df["indicator"].isin(COST_COMPONENTS))
    ].copy()

    # drop import-side FuelCost (double-count)
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

    # Sum across all remaining dims (regions, techs, etc.)
    agg = df.groupby(["year", "scenario", "component"], as_index=False, observed=False)["value"].sum()

    # million EUR -> billion EUR
    agg["value_bEUR"] = agg["value"] / 1000.0
    return agg[["year", "scenario", "component", "value_bEUR"]]


def build_lcoe_components(ind_det_long: pd.DataFrame, cba_long: pd.DataFrame) -> pd.DataFrame:
    """
    LCOE EUR/MWh by cost component using indicator_accounting_detailed (ind_det_long).
    - Drops Transp*/Heat*/FuelImport* techs (sector exclusion you wanted).
    - Drops FuelImport_* rows for FuelCost (double counting fix).
    - Does NOT assume ind_det_long contains accNodesModel; it aggregates across all node dims present.
    """
    if ind_det_long.empty:
        return pd.DataFrame(columns=["year", "scenario", "component", "eur_per_mwh"])
    
    print("ind_det_long columns:", list(ind_det_long.columns))
    print(ind_det_long.head(2).to_string(index=False))

    # --- COSTS (detailed) ---
    df = ind_det_long.copy()
    df = df[(df["year"].isin(YEAR_INTS)) & (df["indicator"].isin(COST_COMPONENTS))].copy()

    tech = df["techs"].astype(str)
    df = df[~tech.str.startswith(("Transp", "Heat", "FuelImport"), na=False)].copy()

    # remove double-counted FuelCost at imports
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

    # assume value is in million EUR (your previous convention)
    df["cost_eur"] = df["value"] * 1e6

    cost_agg = df.groupby(["year", "scenario", "component"], as_index=False, observed=False)["cost_eur"].sum()

    # --- GENERATION (still safe to use accNodesModel here) ---
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
    merged["eur_per_mwh"] = (merged["cost_eur"] / merged["gen_mwh"]).replace([np.inf, -np.inf], np.nan)

    out = merged.groupby(["year", "scenario", "component"], as_index=False, observed=False)["eur_per_mwh"].sum()
    return out


def evenly_spaced_colors_from_cmap(cmap_name, n, start=0.0, stop=1.0):
    cmap = cm.get_cmap(cmap_name)
    return [cmap(x) for x in np.linspace(start, stop, n)]


def plot_costs_and_lcoe(ind_long: pd.DataFrame, ind_det_long: pd.DataFrame, cba_long: pd.DataFrame):
    costs = build_cost_long_from_detailed(ind_det_long)  # <- changed
    lcoe_comp = build_lcoe_components(ind_det_long, cba_long)  

    if costs.empty:
        print("No indicator_accounting loaded; skipping cost/LCOE plots.")
        return

    # Exclude "System cost" here (since you wanted it gone)
    comp_order = ["Fuel cost", "Investment", "OMFix", "Slack cost", "Spill penalty", r"CO$_{2}$ slack"]
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

    leg = shared_legend(fig, comp_order, cmap, x=0.78, y=0.5)
    leg.set_title("Cost component")

    plt.show()


def make_legend(fig, cmap, labels, title, anchor):
    handles = [plt.Line2D([0], [0], color=cmap[l], lw=8) for l in labels]
    leg = fig.legend(
        handles,
        labels,
        title=title,
        loc="center left",
        bbox_to_anchor=anchor,
        frameon=True,
        borderpad=LEGEND_BORDER_PAD,
        labelspacing=LEGEND_LABEL_SPACING,
        handlelength=2.0,
        handletextpad=0.8,
        framealpha=LEGEND_FRAME_ALPHA,
    )
    leg.get_frame().set_linewidth(1.2)
    return leg


def plot_co2_grouped(ind_det_long: pd.DataFrame, YEAR_INTS, SCEN_ORDER):
    sns.set_theme(style="whitegrid")

    # ----------------- Grouping & aliases -----------------
    category_map = {
        # Transport
        "TranspDemand_Fossil_LF": "Transport (fossil liquid fuel)",
        "TranspDemand_REfuel":    "Transport (liquid e-fuel)",
        "TranspDemand_Fossil_CH4":"Transport (fossil methane)",
        "TranspDemand_e-CH4":     "Transport (e-methane)",

        # Heat
        "HeatDemand_Fossil_LF":   "Heat (fossil liquid fuel)",
        "HeatDemand_REfuel":      "Heat (e-fuel)",
        "HeatDemand_Fossil_CH4":  "Heat (fossil methane)",
        "HeatDemand_e-CH4":       "Heat (e-methane)",

        # Electricity generation – gas
        "CCGT":                   "Electricity generation (gas)",
        "OCGT":                   "Electricity generation (gas)",
        "GT":                     "Electricity generation (gas)",

        # Electricity generation – coal
        "Thermal_Coal":           "Electricity generation (coal)",

        # Electricity generation – diesel
        "Thermal_Diesel":         "Electricity generation (diesel)",

        # DAC
        "DAC":                    "Direct air capture",
    }

    # Your colors
    C = {
        # transport blues
        "Transport (fossil liquid fuel)": "#1851B3",
        "Transport (liquid e-fuel)":      "#3F8CF1",
        "Transport (fossil methane)":     "#6599ED",
        "Transport (e-methane)":          "#6C92BB",

        # heat oranges
        "Heat (fossil liquid fuel)":      "#FF8C28",
        "Heat (e-fuel)":                  "#FF962E",
        "Heat (fossil methane)":          "#FFBB00",
        "Heat (e-methane)":               "#FECF58",

        # electricity purples
        "Electricity generation (gas)":   "#C428D9",
        "Electricity generation (coal)":  "#E95CF6",
        "Electricity generation (diesel)":"#F755CC",

        # DAC pastel teal
        "Direct air capture":             "#B3DDD7",
    }

    # Normalise tech names that sometimes appear without underscores etc.
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

    # Legend / stack order (fixed, readable)
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

    # ----------------- Data prep -----------------
    df = ind_det_long.copy()

    # ensure year exists
    if "year" not in df.columns:
        for cand in ["accYears", "years"]:
            if cand in df.columns:
                df = df.rename(columns={cand: "year"})
                break
    df["year"] = pd.to_numeric(df["year"], errors="coerce").astype("Int64")
    df = df[df["year"].isin(YEAR_INTS)].copy()

    df = df[df["indicator"].astype(str) == "CO2_emission"].copy()

    # IMPORTANT: preserve sign; kt -> Mt keeps sign
    df["value_mt"] = df["value"] / 1000.0

    tech_norm = df["techs"].astype(str).map(lambda t: alias.get(t, t))
    df["category"] = tech_norm.map(category_map)

    df = df[df["category"].notna()].copy()

    # aggregate across regions/nodesModel (signed sum)
    agg = df.groupby(["year", "scenario", "category"], as_index=False, observed=False)["value_mt"].sum()

    # quick sanity print: show DAC (signed) so it’s obvious we are not taking abs()
    dac = agg[agg["category"] == "Direct air capture"].sort_values(["scenario", "year"])
    if not dac.empty:
        print("DAC CO2_emission (Mt, signed; negative = net removal):")
        print(dac.to_string(index=False))

    # ----------------- X geometry (5 scenarios per year) -----------------
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

    # ----------------- Pivot (fill missing with 0) -----------------
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

    # keep only categories in requested order, and only those that exist
    cols = [c for c in cat_order if c in piv.columns]
    piv = piv[cols]

    # ----------------- Plot -----------------
    fig, ax = plt.subplots(1, 1, figsize=(18, 6.2))
    fig.subplots_adjust(left=0.08, right=0.78, top=0.98, bottom=0.20)

    bottom = np.zeros(len(x))
    for c in piv.columns:
        vals = piv[c].values
        ax.bar(x, vals, bottom=bottom, width=BAR_WIDTH, color=C.get(c, "#9E9E9E"), edgecolor="none")
        bottom += vals

    ax.grid(axis="y", alpha=0.25)
    ax.set_ylabel("Total CO$_2$ emissions by contributor (Mt CO$_2$)\n(negative = net removal)")

    # scenario labels
    scen_labels = [s for _y in years for s in scens]
    ax.set_xticks(x)
    ax.set_xticklabels(scen_labels, rotation=90, va="top")

    # year labels on a secondary bottom axis
    sec = ax.secondary_xaxis("bottom")
    sec.set_xlim(ax.get_xlim())
    sec.set_xticks(year_centers)
    sec.set_xticklabels([str(y) for y in years], rotation=0)
    sec.spines["bottom"].set_visible(False)
    sec.tick_params(axis="x", pad=YEAR_PAD, length=0)

    # ----------------- Legend (same order as stacks), moved left -----------------
    handles = [plt.Line2D([0], [0], color=C.get(c, "#9E9E9E"), lw=8) for c in piv.columns]
    ax.legend(
        handles,
        list(piv.columns),
        title="Contributor",
        loc="center left",
        bbox_to_anchor=(0.95, 0.50),  # more left than before (still outside axes)
        frameon=True,
        borderpad=0.9,
        labelspacing=0.55,
        handlelength=2.0,
        handletextpad=0.8,
        framealpha=1.0,
    )

    plt.show()


def _hourly_carpet_tech_group(tech: str) -> str | None:
    """
    Map raw hourly 'techs' from GDX commodity_balance to a small set of groups.
    Return None to drop.
    """
    t = str(tech)

    # Drop imports and all demands (plus explicit demand list)
    if t.startswith("FuelImport_"):
        return None

    drop_demands = {
        "Elec_demand",
        "HeatDemand_Bio_LF",
        "HeatDemand_REfuel",
        "HeatDemand_e-CH4",
        "TranspDemand_REfuel",
        "TranspDemand_H2",
        "TranspDemand_e-CH4",
    }
    if t in drop_demands or "Demand" in t:
        return None

    # Grouped generation
    if t.startswith("Thermal") or t in {"GT", "OCGT", "CCGT"}:
        return "Thermal and gas turbines"
    if t.startswith("wind_"):
        return "Wind"
    if t.startswith("pv_"):
        return "Solar PV"

    # Ignore these for the MVP
    if t in {"Geothermal", "Hydro"}:
        return None

    if t in {"H2_FC"}:
        return "H2_FC"

    # Storage
    if t == "Battery":
        return "Battery"
    if t in {"H2_storage"}:
        return "H2 storage"

    # Fuel/chemical conversion (electricity-consuming operation)
    if t in {"Electrolyser", "FTropschSyn", "Methanizer", "DAC"}:
        return t

    return None


def _parse_tm_to_hour(series: pd.Series) -> pd.Series:
    """
    Parse timeModel labels like 'tm1'..'tm8760' to int hour index 1..8760.
    """
    s = series.astype(str).str.strip()
    return pd.to_numeric(s.str.replace("tm", "", regex=False), errors="coerce")


def _get_hourly_daily_matrix_from_gdx(
    results: GDXEval,
    year: int,
    tech_group: str,
    *,
    debug: bool = False,
) -> np.ndarray | None:
    """
    Returns (365, 24) for one year + one aggregated tech_group, summed nationally.
    Uses symbol: commodity_balance (hourly) and filters by accYears first to avoid explosion.
    """
    cb = read_symbol(results, "commodity_balance")
    if cb is None:
        if debug:
            print(f"[hourly] Missing commodity_balance for year={year}, tech_group={tech_group}")
        return None

    df = cb.reset_index().copy()

    # Year filter FIRST
    if "accYears" not in df.columns:
        raise KeyError(f"[hourly] commodity_balance missing 'accYears'. Columns: {list(df.columns)}")
    df["year"] = pd.to_numeric(df["accYears"], errors="coerce").astype("Int64")
    df = df[df["year"] == int(year)].copy()
    if df.empty:
        if debug:
            print(f"[hourly] No rows after year filter year={year}")
        return None

    # timeModel -> hour
    if "timeModel" not in df.columns:
        raise KeyError(f"[hourly] commodity_balance missing 'timeModel'. Columns: {list(df.columns)}")
    df["hour"] = _parse_tm_to_hour(df["timeModel"]).astype("Int64")
    df = df[df["hour"].between(1, 8760)].copy()
    if df.empty:
        if debug:
            print(f"[hourly] No valid tm1..tm8760 rows for year={year}")
        return None

    # Tech grouping + drop rules
    if "techs" not in df.columns:
        raise KeyError(f"[hourly] commodity_balance missing 'techs'. Columns: {list(df.columns)}")
    df["tech_group"] = df["techs"].astype(str).map(_hourly_carpet_tech_group)
    df = df[df["tech_group"].notna()].copy()
    df = df[df["tech_group"] == tech_group].copy()
    if df.empty:
        if debug:
            print(f"[hourly] No rows for tech_group={tech_group}, year={year}")
        return None

    # Commodity filters depending on selected group
    if "commodity" not in df.columns:
        raise KeyError(f"[hourly] commodity_balance missing 'commodity'. Columns: {list(df.columns)}")

    if tech_group == "Battery":
        df = df[df["commodity"].astype(str) == "Elec"].copy()
    if tech_group == "H2 storage":
        df = df[df["commodity"].astype(str) == "H2"].copy()
    if tech_group in {"Electrolyser", "FTropschSyn", "Methanizer", "DAC"}:
        df = df[df["commodity"].astype(str) == "Elec"].copy()

    if df.empty:
        if debug:
            print(f"[hourly] All rows removed by commodity filter for tech_group={tech_group}, year={year}")
        return None

    # National total by hour
    s = (
        df.groupby("hour", observed=False)["value"]
          .sum()
          .reindex(range(1, 8761))
          .fillna(0.0)
    )
    vals = s.values
    if len(vals) != 8760:
        if debug:
            print(f"[hourly] Got {len(vals)} hours (expected 8760) for tech_group={tech_group}, year={year}")
        return None

    return vals.reshape(365, 24)


def load_results_by_scenario_std() -> dict[str, GDXEval]:
    """
    Loads each scenario GDX once and returns a dict keyed by SCEN_ORDER labels:
    {'GP','NT','ELEC+','BIO+','H2+'}.
    """
    out: dict[str, GDXEval] = {}
    for case_dir in CASE_DIRS:
        scen = scenario_name(case_dir)  # you already have this
        res = load_results(case_dir)    # you already have this
        if res is None:
            continue
        out[scen] = res
    return out


def _prepare_hourly_panel_specs(year: int) -> list[tuple[str, str, int]]:
    """
    6 rows:
      1) 2020 baseline from GP, titled '2020 (all)'
      2-6) year for GP/NT/ELEC+/BIO+/H2+
    """
    return [
        ("2020 (all)", "GP", 2020),
        (f"{year} (GP)", "GP", year),
        (f"{year} (NT)", "NT", year),
        (f"{year} (ELEC+)", "ELEC+", year),
        (f"{year} (BIO+)", "BIO+", year),
        (f"{year} (H2+)", "H2+", year),
    ]


def fig_hourly_carpet_6x2_from_gdx(
    *,
    results_by_scen: dict[str, GDXEval],
    year: int,
    left_group: str,
    right_group: str,
    outpng: Path | None = None,
    debug: bool = False,
) -> None:
    """
    6x2 carpet plot:
      rows: 2020(all GP baseline), then year for 5 scenarios
      cols: left_group vs right_group
    """
    sns.set_theme(style="whitegrid")

    specs = _prepare_hourly_panel_specs(year)

    left_mats, right_mats = [], []
    for _title, scen, y in specs:
        res = results_by_scen.get(scen)
        if res is None:
            left_mats.append(None)
            right_mats.append(None)
            continue
        left_mats.append(_get_hourly_daily_matrix_from_gdx(res, y, left_group, debug=debug))
        right_mats.append(_get_hourly_daily_matrix_from_gdx(res, y, right_group, debug=debug))

    valid = [M for M in (left_mats + right_mats) if M is not None]
    if not valid:
        raise RuntimeError(f"No valid hourly matrices for {left_group=} / {right_group=} / {year=}")

    vmax = max(float(np.nanmax(M)) for M in valid)

    fig, axes = plt.subplots(
        6, 2,
        figsize=(8.5 * 16/9, 8.5),
        sharex=True,
        sharey=True,
    )
    fig.subplots_adjust(left=0.10, right=0.88, top=0.94, bottom=0.20, hspace=0.25, wspace=0.08)

    cmap = colormaps.get_cmap("inferno")
    im = None

    # column headers
    axes[0, 0].set_title(f"{specs[0][0]} — {left_group}", fontsize=10)
    axes[0, 1].set_title(f"{specs[0][0]} — {right_group}", fontsize=10)

    for i, ((title, scen, y), ML, MR) in enumerate(zip(specs, left_mats, right_mats)):
        axL, axR = axes[i, 0], axes[i, 1]

        if i != 0:  # keep row title consistent; already set first row titles above
            axL.set_title(f"{title} — {left_group}", fontsize=10)
            axR.set_title(f"{title} — {right_group}", fontsize=10)

        if ML is not None:
            im = axL.imshow(ML, aspect="auto", origin="lower", cmap=cmap, vmin=0.0, vmax=vmax)
        if MR is not None:
            im = axR.imshow(MR, aspect="auto", origin="lower", cmap=cmap, vmin=0.0, vmax=vmax)

    # x-axis only bottom row
    for ax in axes[-1, :]:
        ax.set_xlabel("Hour of day")
        ax.set_xticks([0, 6, 12, 18, 23])
        ax.set_xticklabels(["0:00", "6:00", "12:00", "18:00", "23:00"])
        ax.xaxis.set_label_coords(0.5, -0.40)

    # y ticks
    mid_day = int(round(365 / 2))
    for r in range(6):
        for c in range(2):
            axes[r, c].set_yticks([0, mid_day, 365])
            axes[r, c].set_yticklabels(["0", str(mid_day), "365"])

    fig.text(0.035, 0.5, "Day of year", rotation="vertical", va="center", ha="center", fontsize=11)

    if im is not None:
        cbar = fig.colorbar(im, ax=axes, orientation="vertical", fraction=0.035, pad=0.02)
        cbar.set_label("Operation (GWh/h ≈ GW)")

    if outpng is not None:
        fig.savefig(outpng, dpi=300, bbox_inches="tight")
        print(f"Saved {outpng}")

    plt.show()


def plot_hourly_carpet_mvp_from_gdx(year: int = 2050, outdir: Path | None = None, debug: bool = False) -> None:
    """
    MVP: make 3 figures:
      1) Wind vs Solar PV
      2) H2 storage vs Battery
      3) FTropschSyn vs DAC
    """
    if outdir is None:
        outdir = BASE_DIR / "comparison_outputs"
    outdir.mkdir(parents=True, exist_ok=True)

    results_by_scen = load_results_by_scenario_std()

    # fig_hourly_carpet_6x2_from_gdx(
    #     results_by_scen=results_by_scen,
    #     year=year,
    #     left_group="Wind",
    #     right_group="Solar PV",
    #     outpng=outdir / f"fig_hourly_carpet_{year}_wind_vs_solar.png",
    #     debug=debug,
    # )

    fig_hourly_carpet_6x2_from_gdx(
        results_by_scen=results_by_scen,
        year=year,
        left_group="H2 storage",
        right_group="Battery",
        outpng=outdir / f"fig_hourly_carpet_{year}_h2storage_vs_battery.png",
        debug=debug,
    )

    fig_hourly_carpet_6x2_from_gdx(
        results_by_scen=results_by_scen,
        year=year,
        left_group="FTropschSyn",
        right_group="DAC",
        outpng=outdir / f"fig_hourly_carpet_{year}_ftropsch_vs_dac.png",
        debug=debug,
    )


def drop_doublecounted_fuel_import_fuelcost(df_cost: pd.DataFrame) -> pd.DataFrame:
    """
    Plotting fix: remove FuelCost entries that come from fuel imports (FuelImport_* techs),
    because the same fuel is already charged again at generator activity under FuelCost.

    Applies only to indicator == 'FuelCost'. Leaves other indicators untouched.
    """
    df = df_cost.copy()

    if "indicator" not in df.columns or "techs" not in df.columns:
        raise KeyError(f"Need columns ['indicator','techs']. Got: {list(df.columns)}")

    is_fuelcost = df["indicator"].astype(str) == "FuelCost"
    is_import = df["techs"].astype(str).str.startswith("FuelImport_", na=False)

    return df[~(is_fuelcost & is_import)].copy()



def main() -> None:
    cba_long, caps_long, ind_long, ind_det_long = load_needed_data()

    # plot_capacity_and_generation_grouped(caps_long=caps_long, cba_long=cba_long)
    # plot_fuelconv_and_supply(caps_long=caps_long, cba_long=cba_long)
    plot_costs_and_lcoe(ind_long=ind_long, ind_det_long=ind_det_long, cba_long=cba_long)

    plot_co2_grouped(ind_det_long, YEAR_INTS, SCEN_ORDER)
    # plot_hourly_carpet_mvp_from_gdx(year=2050, outdir=BASE_DIR / "comparison_outputs", debug=False)



if __name__ == "__main__":
    main()

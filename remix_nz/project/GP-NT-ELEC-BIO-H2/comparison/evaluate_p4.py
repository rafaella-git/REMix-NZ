
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import StrMethodFormatter
from remix.framework.tools.gdx import GDXEval



CASE_DIRS = [
    "nz_case_GP_2020-2025-2030-2035-2040-2045-2050",
    "nz_case_NT_2020-2025-2030-2035-2040-2045-2050",
    "nz_case_ELEC+_2020-2025-2030-2035-2040-2045-2050",
    "nz_case_BIO+_2020-2025-2030-2035-2040-2045-2050",
    "nz_case_H2+_2020-2025-2030-2035-2040-2045-2050",
]

USE_END_YEAR_2045 = False
END_YEAR = 2045 if USE_END_YEAR_2045 else 2050
BASE_YEAR = 2020

PLOT_SHOW = True         
WRITE_CSV = True


PRINT_TABLES = True
PRINT_WARNINGS = True
EXCLUDE_BATTERY_FIG1 = True

BASE_DIR = Path(r"C:\Local\REMix\remix_nz\project\GP-NT-ELEC-BIO-H2")
OUT_DIR = BASE_DIR / "comparison" / "outputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)

YEARS_ALL = [2020, 2025, 2030, 2035, 2040, 2045, 2050]
YEAR_INTS = [y for y in YEARS_ALL if y <= END_YEAR]
SCEN_ORDER = ["GP", "NT", "ELEC+", "BIO+", "H2+"]

STORAGE_TECHS = {"Battery", "H2_storage"}

EXTRA_CAPS = [
    ("DAC", "CO2_feed", "DAC"),
    ("Electrolyser", "H2", "Electrolyser"),
    ("Methanizer", "e-CH4", "Methanizer"),
    ("FTropschSyn", "REfuel", "FT"),
]

def tech_color_map(tech_list):
    # distinct categorical colors (good enough until you assign fixed colors)
    cols = sns.color_palette("tab20", n_colors=len(tech_list))
    return {t: cols[i] for i, t in enumerate(tech_list)}


# ============================================================
# helpers
# ============================================================
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


def tech_group(tech: str) -> str:
    t = str(tech)
    if t.startswith("wind_onshore_"):
        return "wind_onshore"
    if t.startswith("wind_offshore_"):
        return "wind_offshore"
    return t


def print_table(title: str, df: pd.DataFrame, max_rows: int = 30):
    if not PRINT_TABLES:
        return
    print("\n" + "=" * 120)
    print(title)
    print("-" * 120)
    if df is None or df.empty:
        print("(empty)")
        return
    with pd.option_context("display.width", 240, "display.max_rows", max_rows, "display.max_columns", 80):
        print(df)


def read_symbol(results, name: str) -> pd.DataFrame | None:
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


def load_results(case_dir: str):
    gdx_path = BASE_DIR / case_dir / "result" / f"{case_dir}.gdx"
    if not gdx_path.is_file():
        if PRINT_WARNINGS:
            print("Missing:", gdx_path)
        return None
    return GDXEval(str(gdx_path))


def cmap_colors(name: str, n: int, lo: float = 0.15, hi: float = 0.95):
    cmap = plt.get_cmap(name)
    if n <= 1:
        return [cmap((lo + hi) / 2)]
    xs = np.linspace(lo, hi, n)
    return [cmap(x) for x in xs]


def month_from_hour_index(hour_1_to_8760: pd.Series) -> pd.Series:
    month_hours = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]) * 24
    edges = np.cumsum(month_hours)
    h = hour_1_to_8760.astype(int).clip(1, 8760)
    return np.searchsorted(edges, h, side="left") + 1


def savefig(path: Path):
    plt.savefig(path, dpi=300, bbox_inches="tight")
    if PLOT_SHOW:
        plt.show()
    else:
        plt.close()
    print("Wrote:", path)


# ============================================================
# Build long tables
# ============================================================
def build_cba_annual_long(results, scen: str) -> pd.DataFrame:
    cba = read_symbol(results, "commodity_balance_annual")
    if cba is None:
        return pd.DataFrame()
    df = cba.reset_index()
    df["scenario"] = scen
    df["year"] = pd.to_numeric(df["accYears"], errors="coerce")
    df = df[df["year"].isin(YEAR_INTS)].copy()
    df["year"] = df["year"].astype(int)
    df["tech_group"] = df["techs"].map(tech_group)
    return df


def build_caps_long(results, scen: str) -> pd.DataFrame:
    caps = read_symbol(results, "converter_caps")
    if caps is None:
        return pd.DataFrame()
    caps = caps[caps["value"] > 0.01].dropna()
    df = caps.reset_index()
    df["scenario"] = scen
    df["year"] = pd.to_numeric(df["accYears"], errors="coerce")
    df = df[df["year"].isin(YEAR_INTS)].copy()
    df["year"] = df["year"].astype(int)
    df["tech_group"] = df["techs"].map(tech_group)
    return df


def build_indicator_long(results, scen: str) -> pd.DataFrame:
    ind = read_symbol(results, "indicator_accounting")
    if ind is None:
        return pd.DataFrame()
    df = ind.reset_index()
    df["scenario"] = scen
    df["year"] = pd.to_numeric(df["accYears"], errors="coerce")
    df = df[df["year"].isin(YEAR_INTS)].copy()
    df["year"] = df["year"].astype(int)
    return df


# ============================================================
# Derived annual datasets
# ============================================================
def global_elec_balance_net(cba_annual_long: pd.DataFrame) -> pd.DataFrame:
    df = cba_annual_long.copy()
    df = df[
        (df["accNodesModel"] == "global")
        & (df["commodity"] == "Elec")
        & (df["balanceType"] == "net")
        & (df["year"].isin(YEAR_INTS))
    ].copy()
    df["tech"] = df["tech_group"]
    return df.groupby(["scenario", "year", "tech"], as_index=False, observed=False)["value"].sum()


def global_caps_elec_total(caps_long: pd.DataFrame) -> pd.DataFrame:
    df = caps_long.copy()
    df = df[
        (df["accNodesModel"] == "global")
        & (df["capType"] == "total")
        & (df["commodity"] == "Elec")
        & (df["year"].isin(YEAR_INTS))
    ].copy()
    df["tech"] = df["tech_group"]
    return df.groupby(["scenario", "year", "tech"], as_index=False, observed=False)["value"].sum()


def global_caps_multioutput_total(caps_long: pd.DataFrame) -> pd.DataFrame:
    df = caps_long.copy()
    df = df[
        (df["accNodesModel"] == "global")
        & (df["capType"] == "total")
        & (df["year"].isin(YEAR_INTS))
    ].copy()

    rename_map = {(t, c): lab for (t, c, lab) in EXTRA_CAPS}

    def label_row(row):
        key = (row["techs"], row["commodity"])
        if key in rename_map:
            return rename_map[key]
        return row["tech_group"]

    df["tech"] = df.apply(label_row, axis=1)

    is_elec = df["commodity"] == "Elec"
    is_extra = df.apply(lambda r: (r["techs"], r["commodity"]) in rename_map, axis=1)
    df = df[is_elec | is_extra].copy()

    return df.groupby(["scenario", "year", "tech"], as_index=False, observed=False)["value"].sum()


# ============================================================
# Better, short prints (insightful only)
# ============================================================
def print_key_checks(bal: pd.DataFrame, ind_long: pd.DataFrame):
    # 1) Net elec balance totals: should be ~0
    tot_bal = (bal.groupby(["scenario", "year"], as_index=False, observed=False)["value"].sum()
               .pivot(index="year", columns="scenario", values="value")
               .reindex(YEAR_INTS))
    print_table("Net Elec balance (sum over techs) [GWh]", tot_bal.round(6), max_rows=20)

    # 2) Indicator quick view (global): show only a small set, avoid massive pivot
    if ind_long is None or ind_long.empty:
        return
    ind = ind_long[ind_long["accNodesModel"] == "global"].copy()
    keep = ["CO2_emission", "FuelCost", "Invest", "OMFix", "SlackCost", "Slack_CO2", "SpillPenalty", "SystemCost"]
    ind = ind[ind["indicator"].isin(keep)].copy()
    piv = (ind.pivot_table(index=["year", "indicator"], columns="scenario", values="value", aggfunc="sum", observed=False)
              .reindex(pd.MultiIndex.from_product([YEAR_INTS, keep], names=["year", "indicator"])))
    print_table("Indicators (global) [native units]", piv.round(4), max_rows=80)


# ============================================================
# Figure 1: Demand vs Supply, Battery removed completely
# ============================================================
def plot_balance_demand_supply_grid(bal: pd.DataFrame, out_png: Path,
                                   label_threshold_gwh: float = 10.0):
    """
    Demand/Supply side-by-side, using abs() for consumption.
    Battery removed entirely (no split).
    Values shown in TWh.
    Legend excludes techs whose total annual magnitude < label_threshold_gwh (GWh).
    """
    sns.set_theme(style="whitegrid")
    layout = [["GP", "NT", None], ["ELEC+", "BIO+", "H2+"]]
    fig, axes = plt.subplots(2, 3, figsize=(14, 7), sharey=True)
    fig.subplots_adjust(left=0.07, right=0.82, top=0.98, bottom=0.18, wspace=0.22, hspace=0.35)

    df = bal.copy()
    df["tech_plot"] = df["tech"].astype(str)

    # Drop Battery altogether
    df = df[df["tech_plot"] != "Battery"].copy()

    # Side + magnitude in GWh first
    df["side"] = np.where(df["value"] >= 0, "Supply", "Consumption")
    df["mag_gwh"] = np.where(df["value"] >= 0, df["value"], -df["value"])

    # Filter techs by total magnitude (GWh) to reduce legend clutter
    tech_tot = (df.groupby("tech_plot", observed=False)["mag_gwh"].sum()
                  .sort_values(ascending=False))
    keep_tech = tech_tot[tech_tot >= label_threshold_gwh].index.tolist()
    df = df[df["tech_plot"].isin(keep_tech)].copy()

    # Convert to TWh for plotting
    df["mag"] = df["mag_gwh"] / 1000.0  # TWh

    def ordered_techs(side_name):
        return (df[df["side"] == side_name]
                .groupby("tech_plot", observed=False)["mag"].sum()
                .sort_values(ascending=False)
                .index.tolist())

    cons_order = ordered_techs("Consumption")
    supp_order = ordered_techs("Supply")

    cons_cmap = tech_color_map(cons_order)
    supp_cmap = tech_color_map(supp_order)

    years = YEAR_INTS
    x_base = np.arange(len(years))
    width = 0.26
    dist = 0.23
    x_cons = x_base - dist
    x_supp = x_base + dist

    ymax = 0.0
    for sc in SCEN_ORDER:
        sub = df[df["scenario"] == sc]
        sup = (sub[sub["side"] == "Supply"]
               .pivot_table(index="year", columns="tech_plot", values="mag", aggfunc="sum", observed=False)
               .reindex(years).fillna(0.0))
        con = (sub[sub["side"] == "Consumption"]
               .pivot_table(index="year", columns="tech_plot", values="mag", aggfunc="sum", observed=False)
               .reindex(years).fillna(0.0))
        ymax = max(ymax, sup.sum(axis=1).max(), con.sum(axis=1).max())
    ymax = ymax * 1.10 if ymax > 0 else 1.0

    for r in range(2):
        for c in range(3):
            sc = layout[r][c]
            ax = axes[r, c]
            if sc is None:
                ax.axis("off")
                continue

            sub = df[df["scenario"] == sc].copy()

            sup = (sub[sub["side"] == "Supply"]
                   .pivot_table(index="year", columns="tech_plot", values="mag", aggfunc="sum", observed=False)
                   .reindex(years).fillna(0.0))
            con = (sub[sub["side"] == "Consumption"]
                   .pivot_table(index="year", columns="tech_plot", values="mag", aggfunc="sum", observed=False)
                   .reindex(years).fillna(0.0))

            sup = sup[[t for t in supp_order if t in sup.columns]]
            con = con[[t for t in cons_order if t in con.columns]]

            bottom = np.zeros(len(years))
            for t in con.columns:
                vals = con[t].values
                ax.bar(x_cons, vals, bottom=bottom, width=width, color=cons_cmap[t], edgecolor="none")
                bottom += vals

            bottom = np.zeros(len(years))
            for t in sup.columns:
                vals = sup[t].values
                ax.bar(x_supp, vals, bottom=bottom, width=width, color=supp_cmap[t], edgecolor="none")
                bottom += vals

            ax.set_title(sc, fontsize=11)
            ax.set_ylim(0, ymax)
            ax.grid(axis="y", alpha=0.25)

            ticks = np.r_[x_cons, x_supp]
            labels = [f"Consumption\n{y}" for y in years] + [f"Supply\n{y}" for y in years]
            order = np.argsort(ticks)
            ax.set_xticks(ticks[order])
            ax.set_xticklabels([labels[i] for i in order], rotation=90)

            ax.yaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))

            if c == 0:
                ax.set_ylabel("Electricity balance, absolute net (TWh)")

    cons_handles = [plt.Line2D([0], [0], color=cons_cmap[t], lw=8) for t in cons_order]
    supp_handles = [plt.Line2D([0], [0], color=supp_cmap[t], lw=8) for t in supp_order]
    fig.legend(cons_handles, cons_order, loc="upper right", bbox_to_anchor=(0.98, 0.92), frameon=True, title="Consumption")
    fig.legend(supp_handles, supp_order, loc="lower right", bbox_to_anchor=(0.98, 0.08), frameon=True, title="Supply")

    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.show()
    print("Wrote:", out_png)

# ============================================================
# Figure 2: End-year stacked capacities by tech with 2-level x-axis (Scenario + Year)
# ============================================================
def plot_elec_caps_all_years_grouped_stacked(cap_elec_gen: pd.DataFrame, out_png: Path):
    """
    Single panel.
    For each year: group of scenario bars (GP, NT, ELEC+, BIO+, H2+).
    Each bar is stacked by tech.
    Hydro + Geothermal forced to the bottom of stacks.
    """
    sns.set_theme(style="whitegrid")
    df = cap_elec_gen.copy()
    df = df[df["year"].isin(YEAR_INTS)].copy()

    # Pivot: (year, scenario) x tech
    piv = (df.pivot_table(index=["year", "scenario"], columns="tech", values="value",
                          aggfunc="sum", observed=False)
             .reindex(pd.MultiIndex.from_product([YEAR_INTS, SCEN_ORDER], names=["year", "scenario"]))
             .fillna(0.0))

    # Tech order: force base techs first, then by total descending
    base_techs = [t for t in ["Hydro", "Geothermal"] if t in piv.columns]
    rest = [t for t in piv.columns if t not in base_techs]
    rest_sorted = piv[rest].sum(axis=0).sort_values(ascending=False).index.tolist()
    tech_order = base_techs + rest_sorted
    piv = piv[tech_order]

    # X positions: groups by year, within year by scenario
    years = YEAR_INTS
    scen = SCEN_ORDER
    x_year = np.arange(len(years))

    width = 0.13
    gap = 0.02
    group_width = len(scen) * width + (len(scen) - 1) * gap
    offsets = np.linspace(-group_width/2 + width/2, group_width/2 - width/2, len(scen))

    # Build color map (still inferno for now, but categorical order)
    colors = cmap_colors("inferno", len(tech_order))
    cmap = {t: colors[i] for i, t in enumerate(tech_order)}

    fig, ax = plt.subplots(1, 1, figsize=(16, 6))
    ax.grid(axis="y", alpha=0.25)

    for i_s, s in enumerate(scen):
        xs = x_year + offsets[i_s]
        bottoms = np.zeros(len(years))
        for t in tech_order:
            vals = []
            for y in years:
                vals.append(float(piv.loc[(y, s), t]) if (y, s) in piv.index else 0.0)
            vals = np.array(vals)
            ax.bar(xs, vals, bottom=bottoms, width=width, color=cmap[t], edgecolor="none")
            bottoms += vals

    # Tick labels: scenario rotated 90 with year below (two-line)
    ticks = []
    labels = []
    for yi, y in enumerate(years):
        for si, s in enumerate(scen):
            ticks.append(x_year[yi] + offsets[si])
            labels.append(f"{s}\n{y}")
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, rotation=90, va="top")

    ax.set_ylabel("Capacity installed of electricity generation technologies (GW)")
    ax.yaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))

    handles = [plt.Line2D([0], [0], color=cmap[t], lw=8) for t in tech_order]
    fig.legend(handles, tech_order, loc="center right", bbox_to_anchor=(0.98, 0.5),
               frameon=True, title="Technology")

    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.show()
    print("Wrote:", out_png)
# ============================================================
# Figure 3
# ============================================================

def plot_caps_multioutput_total_grouped(cap_multi: pd.DataFrame, out_png: Path):
    """
    Single panel: total (sum over tech) capacity per scenario-year.
    Grouped by year; within each year, bars for scenarios.
    """
    sns.set_theme(style="whitegrid")

    df = cap_multi.copy()
    df = df[df["year"].isin(YEAR_INTS)].copy()

    tot = (df.groupby(["year", "scenario"], as_index=False, observed=False)["value"]
             .sum())

    years = YEAR_INTS
    scen = SCEN_ORDER
    x = np.arange(len(years))

    width = 0.13
    gap = 0.02
    group_width = len(scen) * width + (len(scen) - 1) * gap
    offsets = np.linspace(-group_width/2 + width/2, group_width/2 - width/2, len(scen))

    colors = cmap_colors("inferno", len(scen))
    sc_cmap = {s: colors[i] for i, s in enumerate(scen)}

    fig, ax = plt.subplots(1, 1, figsize=(12.5, 4.8))
    ax.grid(axis="y", alpha=0.25)

    for i, s in enumerate(scen):
        vals = []
        for y in years:
            v = tot[(tot["scenario"] == s) & (tot["year"] == y)]["value"]
            vals.append(float(v.iloc[0]) if len(v) else 0.0)
        ax.bar(x + offsets[i], vals, width=width, color=sc_cmap[s], edgecolor="none", label=s)

    ax.set_xticks(x)
    ax.set_xticklabels([str(y) for y in years], rotation=0)
    ax.set_ylabel("Cap")
    ax.yaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))

    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=True)

    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.show()
    print("Wrote:", out_png)

# ============================================================
# Other plots: shorten y-labels
# ============================================================
def plot_grid_stacked_positive(df: pd.DataFrame, out_png: Path, ylab: str, cmap_name: str):
    sns.set_theme(style="whitegrid")
    layout = [["GP", "NT", None], ["ELEC+", "BIO+", "H2+"]]
    fig, axes = plt.subplots(2, 3, figsize=(14, 7), sharey=True)
    fig.subplots_adjust(left=0.07, right=0.82, top=0.98, bottom=0.18, wspace=0.22, hspace=0.35)

    tech_order = (
        df.groupby("tech", observed=False)["value"].sum()
        .sort_values(ascending=False)
        .index.tolist()
    )
    colors = cmap_colors(cmap_name, max(1, len(tech_order)))
    cmap = {t: colors[i] for i, t in enumerate(tech_order)}

    ymax = 0.0
    for sc in SCEN_ORDER:
        s = df[df["scenario"] == sc]
        piv = s.pivot_table(index="year", columns="tech", values="value", aggfunc="sum", observed=False).reindex(YEAR_INTS).fillna(0.0)
        ymax = max(ymax, piv.sum(axis=1).max())
    ymax = ymax * 1.10 if ymax > 0 else 1.0

    for r in range(2):
        for c in range(3):
            sc = layout[r][c]
            ax = axes[r, c]
            if sc is None:
                ax.axis("off")
                continue

            s = df[df["scenario"] == sc].copy()
            piv = s.pivot_table(index="year", columns="tech", values="value", aggfunc="sum", observed=False).reindex(YEAR_INTS).fillna(0.0)
            piv = piv[[t for t in tech_order if t in piv.columns]]

            bottom = np.zeros(len(piv.index))
            for tech in piv.columns:
                vals = piv[tech].values
                ax.bar(piv.index.values, vals, bottom=bottom, width=3.6, color=cmap[tech], edgecolor="none")
                bottom += vals

            ax.set_title(sc, fontsize=11)
            ax.set_ylim(0, ymax)
            ax.set_xticks(YEAR_INTS)
            ax.set_xticklabels([str(y) for y in YEAR_INTS], rotation=90)
            ax.grid(axis="y", alpha=0.25)
            if c == 0:
                ax.set_ylabel(ylab)

    handles = [plt.Line2D([0], [0], color=cmap[t], lw=8) for t in tech_order]
    fig.legend(handles, tech_order, loc="center right", bbox_to_anchor=(0.98, 0.5), frameon=True)

    savefig(out_png)


def plot_elec_caps_vs_generation_endyear(cap_elec: pd.DataFrame, cba_long: pd.DataFrame, out_png: Path, end_year: int):
    sns.set_theme(style="whitegrid")

    cap = cap_elec[cap_elec["year"] == end_year].copy()
    cap_piv = (cap.pivot_table(index="scenario", columns="tech", values="value", aggfunc="sum", observed=False)
               .reindex(SCEN_ORDER).fillna(0.0))

    gen = cba_long.copy()
    gen = gen[
        (gen["accNodesModel"] == "global")
        & (gen["commodity"] == "Elec")
        & (gen["balanceType"] == "net")
        & (gen["year"] == end_year)
        & (gen["value"] > 0)
    ].copy()
    gen["tech"] = gen["tech_group"]
    gen["value_twh"] = gen["value"] / 1000.0  # GWh -> TWh

    gen_piv = (gen.groupby(["scenario", "tech"], as_index=False, observed=False)["value_twh"].sum()
               .pivot(index="scenario", columns="tech", values="value_twh")
               .reindex(SCEN_ORDER).fillna(0.0))

    techs = sorted(set(cap_piv.columns) | set(gen_piv.columns))
    score = cap_piv.reindex(columns=techs, fill_value=0).sum(axis=0) + gen_piv.reindex(columns=techs, fill_value=0).sum(axis=0)
    tech_order = score.sort_values(ascending=False).index.tolist()

    colors = cmap_colors("inferno", len(tech_order))
    cmap = {t: colors[i] for i, t in enumerate(tech_order)}

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    fig.subplots_adjust(left=0.07, right=0.83, top=0.98, bottom=0.26, wspace=0.25)
    x = np.arange(len(SCEN_ORDER))

    ax = axes[0]
    bottom = np.zeros(len(x))
    for t in tech_order:
        vals = cap_piv.reindex(columns=tech_order, fill_value=0)[t].values
        ax.bar(x, vals, bottom=bottom, width=0.75, color=cmap[t], edgecolor="none")
        bottom += vals
    ax.set_ylabel("Cap (GW)")
    ax.set_xticks(x)
    ax.set_xticklabels(SCEN_ORDER)
    ax.grid(axis="y", alpha=0.25)

    ax = axes[1]
    bottom = np.zeros(len(x))
    for t in tech_order:
        vals = gen_piv.reindex(columns=tech_order, fill_value=0)[t].values
        ax.bar(x, vals, bottom=bottom, width=0.75, color=cmap[t], edgecolor="none")
        bottom += vals
    ax.set_ylabel("Gen (TWh)")
    ax.set_xticks(x)
    ax.set_xticklabels(SCEN_ORDER)
    ax.grid(axis="y", alpha=0.25)

    handles = [plt.Line2D([0], [0], color=cmap[t], lw=8) for t in tech_order]
    fig.legend(handles, tech_order, loc="center right", bbox_to_anchor=(0.98, 0.5), frameon=True)

    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.show()
    print("Wrote:", out_png)

def plot_fuelconv_caps_vs_supply_endyear(cap_multi: pd.DataFrame, cba_long: pd.DataFrame, out_png: Path, end_year: int):
    sns.set_theme(style="whitegrid")

    # LEFT: conversion capacities (GW_output)
    cap = cap_multi[cap_multi["year"] == end_year].copy()

    left_keep = ["Water electrolysis", "Methanizer", "FT", "DAC"]
    rename_left = {
        "Water electrolysis": "Water electrolysis",
        "Methanizer": "Methanation",
        "FT": "Fischer Tropsch",
        "DAC": "DAC",
    }



    cap = cap[cap["tech"].isin(left_keep)].copy()
    cap["tech"] = cap["tech"].map(rename_left)
    cap_piv = (cap.pivot_table(index="scenario", columns="tech", values="value", aggfunc="sum", observed=False)
               .reindex(SCEN_ORDER).fillna(0.0))

    # RIGHT: fuels & chemicals supply (TWh) from commodity_balance_annual
    df = cba_long.copy()
    df = df[(df["year"] == end_year) & (df["balanceType"] == "net") & (df["value"] > 0)].copy()

    # warn + exclude SlackFuel_*
    slack = df[df["techs"].astype(str).str.startswith("SlackFuel_")].copy()
    if not slack.empty:
        slack_sum = slack.groupby("scenario", as_index=False, observed=False)["value"].sum()
        for _, row in slack_sum.iterrows():
            if row["value"] > 1.0:
                print(f"WARNING: {row['scenario']} SlackFuel_* supply is {row['value']:.3f} GWh (> 1 GWh).")
    df = df[~df["techs"].astype(str).str.startswith("SlackFuel_")].copy()

    def supply_group(row):
        t = str(row["techs"])
        comm = str(row["commodity"])

        # direct commodities
        if comm == "e-CH4":
            return "e-methane"
        if comm == "REfuel":
            return "e-FTL fuels"
        if comm == "H2":
            return "e-hydrogen (H2)"

        # imports mapped via tech name
        if t in ["FuelImport_Bio_LF", "FuelImport_Biofuel"]:
            return "Biofuels liquid"
        if t == "FuelImport_Bio_gas":
            return "Biogas"
        if t in ["FuelImport_CH4", "FuelImport_Fossil_CH4"]:
            return "Fossil methane"
        if t == "FuelImport_Coal":
            return "Fossil coal"
        if t in ["FuelImport_Diesel", "FuelImport_Fossil_LF"]:
            return "Fossil oil"
        return None

    df["supply_group"] = df.apply(supply_group, axis=1)
    df = df[df["supply_group"].notna()].copy()
    df["value_twh"] = df["value"] / 1000.0

    sup_piv = (df.groupby(["scenario", "supply_group"], as_index=False, observed=False)["value_twh"].sum()
               .pivot(index="scenario", columns="supply_group", values="value_twh")
               .reindex(SCEN_ORDER).fillna(0.0))

    # shared color map across both panels (union of categories)
    left_items = list(cap_piv.columns)
    right_items = list(sup_piv.columns)
    items = left_items + [x for x in right_items if x not in left_items]
    colors = cmap_colors("inferno", len(items))
    cmap = {k: colors[i] for i, k in enumerate(items)}

    fig, axes = plt.subplots(1, 2, figsize=(14, 4.8))
    fig.subplots_adjust(left=0.07, right=0.83, top=0.98, bottom=0.26, wspace=0.25)
    x = np.arange(len(SCEN_ORDER))

    ax = axes[0]
    bottom = np.zeros(len(x))
    for k in left_items:
        vals = cap_piv[k].values
        ax.bar(x, vals, bottom=bottom, width=0.75, color=cmap[k], edgecolor="none")
        bottom += vals
    ax.set_ylabel(r"Total installed capacity for fuel conversion [GW$_{output}$]")
    ax.set_xticks(x)
    ax.set_xticklabels(SCEN_ORDER)
    ax.grid(axis="y", alpha=0.25)

    ax = axes[1]
    bottom = np.zeros(len(x))
    for k in right_items:
        vals = sup_piv[k].values
        ax.bar(x, vals, bottom=bottom, width=0.75, color=cmap[k], edgecolor="none")
        bottom += vals
    ax.set_ylabel("Supply of fuels and chemicals [TWh]")
    ax.set_xticks(x)
    ax.set_xticklabels(SCEN_ORDER)
    ax.grid(axis="y", alpha=0.25)

    handles = [plt.Line2D([0], [0], color=cmap[k], lw=8) for k in items]
    fig.legend(handles, items, loc="center right", bbox_to_anchor=(0.98, 0.5), frameon=True)

    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.show()
    print("Wrote:", out_png)

def plot_elec_generation_all_years_grouped_stacked(cba_annual_long: pd.DataFrame, out_png: Path):
    """
    Elec generation by tech (TWh), from commodity_balance_annual:
      commodity==Elec, balanceType==net, value>0, accNodesModel==global
    Grouped by year, scenario bars within year, stacked by tech.
    Hydro + Geothermal forced to base.
    """
    sns.set_theme(style="whitegrid")
    df = cba_annual_long.copy()
    df = df[
        (df["accNodesModel"] == "global")
        & (df["commodity"] == "Elec")
        & (df["balanceType"] == "net")
        & (df["value"] > 0)
        & (df["year"].isin(YEAR_INTS))
    ].copy()
    df["tech"] = df["tech_group"]
    df["value_twh"] = df["value"] / 1000.0

    piv = (df.pivot_table(index=["year", "scenario"], columns="tech", values="value_twh",
                          aggfunc="sum", observed=False)
             .reindex(pd.MultiIndex.from_product([YEAR_INTS, SCEN_ORDER], names=["year", "scenario"]))
             .fillna(0.0))

    base_techs = [t for t in ["Hydro", "Geothermal"] if t in piv.columns]
    rest = [t for t in piv.columns if t not in base_techs]
    rest_sorted = piv[rest].sum(axis=0).sort_values(ascending=False).index.tolist()
    tech_order = base_techs + rest_sorted
    piv = piv[tech_order]

    years = YEAR_INTS
    scen = SCEN_ORDER
    x_year = np.arange(len(years))
    width = 0.13
    gap = 0.02
    group_width = len(scen) * width + (len(scen) - 1) * gap
    offsets = np.linspace(-group_width/2 + width/2, group_width/2 - width/2, len(scen))

    colors = cmap_colors("inferno", len(tech_order))
    cmap = {t: colors[i] for i, t in enumerate(tech_order)}

    fig, ax = plt.subplots(1, 1, figsize=(16, 6))
    ax.grid(axis="y", alpha=0.25)

    for i_s, s in enumerate(scen):
        xs = x_year + offsets[i_s]
        bottoms = np.zeros(len(years))
        for t in tech_order:
            vals = []
            for y in years:
                vals.append(float(piv.loc[(y, s), t]) if (y, s) in piv.index else 0.0)
            vals = np.array(vals)
            ax.bar(xs, vals, bottom=bottoms, width=width, color=cmap[t], edgecolor="none")
            bottoms += vals

    ticks, labels = [], []
    for yi, y in enumerate(years):
        for si, s in enumerate(scen):
            ticks.append(x_year[yi] + offsets[si])
            labels.append(f"{s}\n{y}")
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, rotation=90, va="top")

    ax.set_ylabel(f"Electricity generation technologies ({END_YEAR} horizon shown), TWh")
    ax.yaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))

    handles = [plt.Line2D([0], [0], color=cmap[t], lw=8) for t in tech_order]
    fig.legend(handles, tech_order, loc="center right", bbox_to_anchor=(0.98, 0.5),
               frameon=True, title="Technology")

    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.show()
    print("Wrote:", out_png)

def plot_flh_vre_side_by_side(cba_annual_long: pd.DataFrame, out_png: Path, scenario: str = "GP"):
    """
    FLH for VRE from commodity_balance_annual (balanceType='flh', commodity='Elec').
    Bars are side-by-side (not stacked), grouped by year.
    """
    sns.set_theme(style="whitegrid")

    df = cba_annual_long.copy()
    df = df[
        (df["accNodesModel"] == "global")
        & (df["commodity"] == "Elec")
        & (df["balanceType"] == "flh")
        & (df["scenario"] == scenario)
        & (df["year"].isin(YEAR_INTS))
    ].copy()

    # Keep each onshore/offshore wind variant separately by using techs,
    # and PV separately too
    keep = df["techs"].astype(str).str.startswith("wind_onshore_") | \
           df["techs"].astype(str).str.startswith("wind_offshore_") | \
           df["techs"].astype(str).isin(["pv_central_fixed", "pv_decentral"])
    df = df[keep].copy()

    # Average FLH over nodesModel/regions if present
    s = (df.groupby(["year", "techs"], as_index=False, observed=False)["value"]
           .mean())

    tech_order = (s.groupby("techs", observed=False)["value"].mean()
                    .sort_values(ascending=False).index.tolist())

    years = YEAR_INTS
    techs = tech_order

    x = np.arange(len(years))
    width = min(0.8 / max(1, len(techs)), 0.08)  # keep reasonable
    offsets = np.linspace(-0.4 + width/2, 0.4 - width/2, len(techs))

    colors = sns.color_palette("tab20", n_colors=len(techs))
    cmap = {t: colors[i] for i, t in enumerate(techs)}

    fig, ax = plt.subplots(1, 1, figsize=(14, 4.8))
    ax.grid(axis="y", alpha=0.25)

    for i, t in enumerate(techs):
        vals = []
        for y in years:
            v = s[(s["year"] == y) & (s["techs"] == t)]["value"]
            vals.append(float(v.iloc[0]) if len(v) else 0.0)
        ax.bar(x + offsets[i], vals, width=width, color=cmap[t], edgecolor="none", label=t)

    ax.set_xticks(x)
    ax.set_xticklabels([str(y) for y in years])
    ax.set_ylabel(f"FLH (h), {scenario}")
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=True, title="Technology")

    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.show()
    print("Wrote:", out_png)


# ============================================================
# Main
# ============================================================
def main():
    results_by_scen = {}
    cba_list, caps_list, ind_list = [], [], []

    for case_dir in CASE_DIRS:
        scen = scenario_name(case_dir)
        res = load_results(case_dir)
        if res is None:
            continue
        results_by_scen[scen] = res

        cba = build_cba_annual_long(res, scen)
        caps = build_caps_long(res, scen)
        ind = build_indicator_long(res, scen)

        if not cba.empty:
            cba_list.append(cba)
        if not caps.empty:
            caps_list.append(caps)
        if not ind.empty:
            ind_list.append(ind)

    if not cba_list or not caps_list:
        print("No data loaded.")
        return

    cba_annual_long = pd.concat(cba_list, ignore_index=True)
    caps_long = pd.concat(caps_list, ignore_index=True)
    ind_long = pd.concat(ind_list, ignore_index=True) if ind_list else pd.DataFrame()

    # consistent ordering
    cba_annual_long["scenario"] = pd.Categorical(cba_annual_long["scenario"], categories=SCEN_ORDER, ordered=True)
    caps_long["scenario"] = pd.Categorical(caps_long["scenario"], categories=SCEN_ORDER, ordered=True)
    if not ind_long.empty:
        ind_long["scenario"] = pd.Categorical(ind_long["scenario"], categories=SCEN_ORDER, ordered=True)

    # derived tables
    bal = global_elec_balance_net(cba_annual_long)
    cap_elec = global_caps_elec_total(caps_long)
    cap_elec_gen = cap_elec[~cap_elec["tech"].isin(STORAGE_TECHS)].copy()
    cap_multi = global_caps_multioutput_total(caps_long)

    # prints (minimal but useful)
    print(f"\nRUN CONFIG: BASE_YEAR={BASE_YEAR}, END_YEAR={END_YEAR}, YEARS={YEAR_INTS}")
    print_key_checks(bal, ind_long)

    # # CSVs
    # if WRITE_CSV:
    #     cba_annual_long.to_csv(OUT_DIR / f"commodity_balance_annual__long__to{END_YEAR}.csv", index=False)
    #     caps_long.to_csv(OUT_DIR / f"converter_caps__long__to{END_YEAR}.csv", index=False)
    #     if not ind_long.empty:
    #         ind_long.to_csv(OUT_DIR / f"indicator_accounting__long__to{END_YEAR}.csv", index=False)
    #     bal.to_csv(OUT_DIR / f"global_elec_balance_net_by_tech__to{END_YEAR}.csv", index=False)
    #     cap_elec.to_csv(OUT_DIR / f"global_caps_elec_total_by_tech__to{END_YEAR}.csv", index=False)
    #     cap_multi.to_csv(OUT_DIR / f"global_caps_multioutput_total__to{END_YEAR}.csv", index=False)
    #     print("Wrote CSVs to:", OUT_DIR)

    # === FIGURES ===
    # plot_balance_demand_supply_grid(
    #     bal=bal,
    #     out_png=OUT_DIR / f"fig1_elec_balance_demand_supply__to{END_YEAR}.png",
    # )


    plot_elec_caps_all_years_grouped_stacked(
        cap_elec_gen,
        out_png=OUT_DIR / f"fig2_elec_caps_all_years_grouped_stacked__to{END_YEAR}.png",
    )




    # converter caps grid (short y-label)
    plot_grid_stacked_positive(
        df=cap_multi,
        out_png=OUT_DIR / f"fig3_caps_multioutput__to{END_YEAR}.png",
        ylab="Cap",
        cmap_name="inferno",
    )

    plot_fuelconv_caps_vs_supply_endyear(
        cap_multi=cap_multi,
        cba_long=cba_annual_long,
        out_png=OUT_DIR / f"fig5_fuelconv_caps_vs_supply__{END_YEAR}__inferno.png",
        end_year=END_YEAR,
    )

    plot_flh_vre_side_by_side(
        cba_annual_long,
        out_png=OUT_DIR / f"fig_flh_vre_side_by_side__{END_YEAR}__GP.png",
        scenario="GP",
    )


if __name__ == "__main__":
    main()

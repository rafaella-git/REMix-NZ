# compare_runs_win.py
# -----------------------------------------------------------------------------
# Windows-friendly REMix comparison (Base vs Scenario) + optional weekly plots.
# - Fixes gdxpds import warning (import order).
# - Normalizes value strings so 'Elec'/'H2' work.
# - Adds readable stats + optional REMix tutorial-style weekly plots.
# -----------------------------------------------------------------------------

# 1) Import gdxpds BEFORE pandas to avoid the warning
import gdxpds  # must be first

# Now the rest
import os
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Optional: use REMix's high-level reader if available (for weekly plots)
try:
    from remix.framework import GDXEval
    HAS_GDXEVAL = True
except Exception:
    HAS_GDXEVAL = False

# --------------------------- User settings -----------------------------------
GROUP       = "hadi"
BASE_CASE   = "pypsa-cascade"
SCENARIO    = "dry-year"
YEAR        = "2030"
LABEL_BASE  = "Base"
LABEL_SCEN  = "Dry-Year"

# Pin to your scenario file (dispatch/dry-year/...)
SCENARIO_GDX = Path(
    f"../project/{GROUP}/{BASE_CASE}/dispatch/{SCENARIO}/{BASE_CASE}_{SCENARIO}_dispatch.gdx"
)

SAVE_OUTPUT = True
IMG_DPI     = 180
MAKE_WEEKLY_PLOTS = True  # set False to skip the extra tutorial-style charts

# --------------------------- Paths -------------------------------------------
os.chdir(Path(__file__).parent.resolve())

base_gdx_path = Path(f"../project/{GROUP}/{BASE_CASE}/result/{BASE_CASE}_opt.gdx")
if not base_gdx_path.exists():
    raise FileNotFoundError(f"Base GDX not found: {base_gdx_path}")

scenario_gdx_path = SCENARIO_GDX
if not scenario_gdx_path.exists():
    raise FileNotFoundError(f"Scenario GDX not found: {scenario_gdx_path}")

OUTDIR = Path("_compare_out")
if SAVE_OUTPUT:
    OUTDIR.mkdir(exist_ok=True)

# --------------------------- Helpers -----------------------------------------
def load_gdx(path: Path) -> dict:
    data = gdxpds.to_dataframes(str(path))
    print(f"Loaded {len(data)} tables from {path.name}")
    return data

def _lower_cols(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = df.columns.str.lower()
    return df

def _norm_strings(df: pd.DataFrame, cols) -> pd.DataFrame:
    df = df.copy()
    for c in cols:
        if c in df.columns:
            df[c] = df[c].astype(str).str.strip().str.lower()
    return df

def _fmt(x, d=2): return f"{x:.{d}f}"

def _print_section(t):
    print("\n" + t)
    print("-" * len(t))

def _save_csv(df: pd.DataFrame, name: str):
    if SAVE_OUTPUT:
        p = OUTDIR / f"{name}.csv"
        df.to_csv(p, index=False)
        print(f"  • Saved: {p}")

def _annotate_barh(ax: plt.Axes):
    for c in ax.containers:
        ax.bar_label(c, fmt="%.2f", padding=4)
    ax.set_axisbelow(True)
    ax.grid(True, axis="x", linestyle="--", alpha=0.5)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

def _report_table(name: str, df: pd.DataFrame):
    d = _lower_cols(df.reset_index())
    d = _norm_strings(d, ["accyears", "commodity", "balancetype", "accnodesmodel", "captype"])
    yrs   = sorted(d["accyears"].unique())        if "accyears"      in d else []
    comms = sorted(d["commodity"].unique())       if "commodity"     in d else []
    bals  = sorted(d["balancetype"].unique())     if "balancetype"   in d else []
    nodes = sorted(d["accnodesmodel"].unique())   if "accnodesmodel" in d else []
    _print_section(f"TABLE REPORT: {name}")
    print(f"Rows: {len(d)}")
    if yrs:   print(f"Years: {', '.join(yrs[:12])}" + (" ..." if len(yrs) > 12 else ""))
    if comms: print(f"Commodities: {', '.join(comms[:12])}" + (" ..." if len(comms) > 12 else ""))
    if bals:  print(f"Balance types: {', '.join(bals)}")
    if nodes: print(f"Nodes (sample): {', '.join(nodes[:12])}" + (" ..." if len(nodes) > 12 else ""))

# ------------------------ Comparisons -----------------------------------------
def compare_converter_capacity(cap_base, cap_scen, year, label1, label2):
    d1 = _norm_strings(_lower_cols(cap_base.reset_index()), ["accyears", "captype", "commodity", "accnodesmodel"])
    d2 = _norm_strings(_lower_cols(cap_scen.reset_index()), ["accyears", "captype", "commodity", "accnodesmodel"])

    filt = lambda d: d[
        (d["accyears"] == str(year)) &
        (d["captype"] == "total") &
        (d["commodity"] == "elec") &
        (~d["accnodesmodel"].eq("global"))
    ]
    d1 = filt(d1).groupby("techs", as_index=False)["value"].sum()
    d2 = filt(d2).groupby("techs", as_index=False)["value"].sum()

    comp = pd.merge(d1, d2, on="techs", how="outer",
                    suffixes=(f"_{label1}", f"_{label2}")).fillna(0)
    comp["Difference"] = comp[f"value_{label2}"] - comp[f"value_{label1}"]
    comp = comp.sort_values("Difference", ascending=False).reset_index(drop=True)

    _print_section("Installed Capacity Comparison (GW)")
    totals = {label1: comp[f"value_{label1}"].sum(), label2: comp[f"value_{label2}"].sum()}
    change_abs = totals[label2] - totals[label1]
    change_pct = 0 if totals[label1] == 0 else 100 * change_abs / totals[label1]

    view = comp.rename(columns={
        "techs": "Technology",
        f"value_{label1}": f"{label1} (GW)",
        f"value_{label2}": f"{label2} (GW)"
    })
    print(view.to_string(index=False, formatters={
        f"{label1} (GW)": lambda x: _fmt(x, 2),
        f"{label2} (GW)": lambda x: _fmt(x, 2),
        "Difference":    lambda x: _fmt(x, 2)
    }))

    print(f"\nTotal {label1}: {_fmt(totals[label1])} GW")
    print(f"Total {label2}: {_fmt(totals[label2])} GW")
    print(f"Change: {change_abs:+.2f} GW ({change_pct:+.2f}%)")

    print("\nTop increases (GW):")
    print(view.nlargest(5, "Difference")[["Technology", "Difference"]].to_string(index=False))
    print("\nTop decreases (GW):")
    print(view.nsmallest(5, "Difference")[["Technology", "Difference"]].to_string(index=False))

    _save_csv(view, "capacity_comparison_gw")
    return comp

def compare_commodity_balance(bal_base, bal_scen, year, label1, label2):
    d1 = _norm_strings(_lower_cols(bal_base.reset_index()), ["accyears", "accnodesmodel", "balancetype", "commodity"])
    d2 = _norm_strings(_lower_cols(bal_scen.reset_index()), ["accyears", "accnodesmodel", "balancetype", "commodity"])

    filt = lambda d: d[
        (d["accyears"] == str(year)) &
        (~d["accnodesmodel"].eq("global")) &
        (d["balancetype"] == "net")
    ]
    d1 = filt(d1)
    d2 = filt(d2)

    cascade_mode = bool(d1["commodity"].str.contains("_in").any())
    print("\nInflows mode:", "Detailed cascade inflows detected" if cascade_mode
          else "Aggregated inflow (water_in)")

    def agg(df):
        if cascade_mode:
            df = df.copy()
            df["commodity_group"] = df["commodity"].where(~df["commodity"].str.contains("_in"),
                                                          "Total_Hydro_Inflow")
        else:
            df = df.assign(commodity_group=df["commodity"])
        return df.groupby("commodity_group", as_index=False)["value"].sum()

    d1s = agg(d1)
    d2s = agg(d2)

    comp = pd.merge(d1s, d2s, on="commodity_group", how="outer",
                    suffixes=(f"_{label1}", f"_{label2}")).fillna(0)
    comp["Difference"] = comp[f"value_{label2}"] - comp[f"value_{label1}"]
    comp = comp.sort_values("Difference", ascending=False).reset_index(drop=True)

    _print_section("Commodity Balance Comparison (TWh)")
    view = comp.rename(columns={
        "commodity_group": "Commodity",
        f"value_{label1}": f"{label1} (TWh)",
        f"value_{label2}": f"{label2} (TWh)"
    })
    print(view.to_string(index=False, formatters={
        f"{label1} (TWh)": lambda x: _fmt(x, 2),
        f"{label2} (TWh)": lambda x: _fmt(x, 2),
        "Difference":     lambda x: _fmt(x, 2)
    }))

    tot1 = comp[f"value_{label1}"].sum()
    tot2 = comp[f"value_{label2}"].sum()
    print(f"\nTotal flow {label1}: {_fmt(tot1)} TWh")
    print(f"Total flow {label2}: {_fmt(tot2)} TWh")
    print(f"Change: {tot2 - tot1:+.2f} TWh")

    print("\nTop increases (TWh):")
    print(view.nlargest(5, "Difference")[["Commodity", "Difference"]].to_string(index=False))
    print("\nTop decreases (TWh):")
    print(view.nsmallest(5, "Difference")[["Commodity", "Difference"]].to_string(index=False))

    _save_csv(view, "commodity_balance_comparison_twh")
    return comp

def compare_global_dispatch(bal_base, bal_scen, label1, label2, year):
    def extract(df):
        d = _norm_strings(_lower_cols(df.reset_index()),
                          ["accnodesmodel", "commodity", "balancetype", "accyears"])
        d = d[
            (d["accnodesmodel"] == "global") &
            (d["commodity"].isin(["elec", "h2"])) &
            (d["balancetype"] == "net") &
            (d["accyears"] == str(year))
        ]
        return d.groupby(["techs", "commodity"], as_index=False)["value"].sum()

    d1 = extract(bal_base)
    d2 = extract(bal_scen)
    comp = pd.merge(d1, d2, on=["techs", "commodity"], how="outer",
                    suffixes=(f"_{label1}", f"_{label2}")).fillna(0)
    comp["Difference"] = comp[f"value_{label2}"] - comp[f"value_{label1}"]
    comp = comp.sort_values("Difference", ascending=False).reset_index(drop=True)

    _print_section("Global Dispatch Comparison (Elec + H2, TWh)")
    view = comp.rename(columns={
        "techs": "Technology",
        "commodity": "Commodity",
        f"value_{label1}": f"{label1} (TWh)",
        f"value_{label2}": f"{label2} (TWh)"
    })
    print(view.to_string(index=False, formatters={
        f"{label1} (TWh)": lambda x: _fmt(x, 2),
        f"{label2} (TWh)": lambda x: _fmt(x, 2),
        "Difference":     lambda x: _fmt(x, 2)
    }))

    _save_csv(view, "global_dispatch_comparison_twh")
    return comp

# ------------------------ Plotting --------------------------------------------
def plot_results(cap_df, bal_df, disp_df, label1, label2):
    sns.set_theme(style="whitegrid", context="talk")
    palette = sns.color_palette("colorblind", 2)

    def _maybe_plot(df, id_col, val_cols, title, xlab, fname):
        if df is None or df.empty:
            print(f"  • Skipping plot '{title}' (no data after filters).")
            return
        fig, ax = plt.subplots(figsize=(10, max(5, 0.4 * len(df))))
        plot = df.copy().sort_values("Difference").rename(columns={id_col: "Label"})
        plot_m = plot.melt(id_vars=["Label"], value_vars=val_cols, var_name="Scenario", value_name=xlab)
        plot_m["Scenario"] = plot_m["Scenario"].map({val_cols[0]: label1, val_cols[1]: label2})
        sns.barplot(data=plot_m, y="Label", x=xlab, hue="Scenario", palette=palette, ax=ax, errorbar=None)
        ax.set_title(title)
        _annotate_barh(ax)
        ax.legend(title="Scenario", frameon=False, loc="lower right")
        fig.tight_layout()
        if SAVE_OUTPUT:
            fig.savefig(OUTDIR / f"{fname}.png", dpi=IMG_DPI)

    _maybe_plot(cap_df,  "techs",           [f"value_{label1}", f"value_{label2}"], "Installed Capacity (GW)", "GW",  "plot_capacity_gw")
    _maybe_plot(bal_df,  "commodity_group", [f"value_{label1}", f"value_{label2}"], "Commodity Balance (TWh)", "TWh", "plot_commodity_balance_twh")
    _maybe_plot(disp_df, "techs",           [f"value_{label1}", f"value_{label2}"], "Global Dispatch (TWh)",   "TWh", "plot_global_dispatch_twh")

    plt.show()

# ----------------------- Weekly plots (tutorial-style) ------------------------
def weekly_extremes_plot(gdx_path: Path, year: str):
    if not HAS_GDXEVAL:
        print("Skipping weekly plots: remix.framework.GDXEval not available.")
        return
    # Load and follow the tutorial logic
    results = GDXEval(str(gdx_path))
    idx = pd.IndexSlice

    # Capacities just to prove the index works (optional print)
    caps = results["converter_caps"]
    print("\nSample capacities (Elec/total):")
    try:
        print(caps.loc[idx[:, :, :, "Elec", "total"], :].head(8))
    except Exception:
        pass  # indexing can differ across cases; it's just a demo print

    # Commodity balances
    commodities = results["commodity_balance"]

    # Generation (positive) and demand (negative) for the chosen year and node model
    # Adjust 'R1_model' if your model uses a different accNodesModel.
    node_model = "R1_model"

    generation = (
        commodities[commodities > 0]
        .loc[idx[:, node_model, year, :, "Elec"], :]
        .dropna()
        .groupby(["timeModel", "techs"])
        .sum()
        .unstack("techs")
        .fillna(0)
    )
    generation.columns = generation.columns.get_level_values(1)

    demand = (
        commodities[commodities < 0]
        .loc[idx[:, node_model, year, :, "Elec"], :]
        .dropna()
        .groupby(["timeModel", "techs"])
        .sum()
        .unstack("techs")
        .fillna(0)
    )
    demand.columns = demand.columns.get_level_values(1)

    # Rolling weekly extremes
    hours_per_interval = 168
    # If you only want renewables, subset columns here (e.g., ["PV","WindOnshore",...])
    gen_total = generation.sum(axis=1)
    rolling_mean = gen_total.rolling(hours_per_interval).mean()
    mean_max = int(rolling_mean.argmax())
    mean_min = int(rolling_mean.argmin())

    technology_colors = None  # let matplotlib choose; set a dict to force colors
    background_color = "#fffaeb"

    def _plot_window(end_hour: int, title: str, fname: str):
        start = max(0, end_hour - hours_per_interval)
        timeslice = range(start, end_hour)
        fig, ax1 = plt.subplots(figsize=(10, 6))
        fig.patch.set_facecolor(background_color)
        ax1.set_title(title)
        ax1.set_facecolor(background_color)
        generation.iloc[timeslice].plot.area(stacked=True, ax=ax1, color=technology_colors)
        plt.legend(loc=(0.0, 1.05), ncols=3)
        plt.ylabel("Generation in GWh_el")

        ax2 = ax1.twinx()
        demand.iloc[timeslice].mul(-1).plot(kind="line", ax=ax2, color="black")
        plt.legend(loc=(0.8, 1.05))
        plt.ylabel("Demand in GWh_el")
        ax2.set_ylim(ax1.get_ylim())
        fig.tight_layout()
        if SAVE_OUTPUT:
            fig.savefig(OUTDIR / f"{fname}.png", dpi=IMG_DPI)
        plt.show()

    print(f"\nWeek with highest renewable feed-in ends in hour {mean_max}.")
    print(f"Week with lowest renewable feed-in ends in hour {mean_min}.")
    _plot_window(mean_max,  "Week with the highest renewable feed-in", "week_highest_feed_in")
    _plot_window(mean_min,  "Week with the lowest renewable feed-in",  "week_lowest_feed_in")

# ------------------------ Main ------------------------------------------------
def main():
    print("\n--- Comparing REMix runs (Windows) ---")
    print(f"Base     : {base_gdx_path}")
    print(f"Scenario : {scenario_gdx_path}")
    print(f"Year     : {YEAR}")

    data_base = load_gdx(base_gdx_path)
    data_scen = load_gdx(scenario_gdx_path)

    # Reports (helps catch naming/case issues)
    for name, tbl in [("converter_caps", data_base.get("converter_caps")),
                      ("commodity_balance_annual (base)", data_base.get("commodity_balance_annual")),
                      ("commodity_balance_annual (scen)", data_scen.get("commodity_balance_annual"))]:
        if tbl is None:
            print(f"WARNING: missing table '{name}'")
        else:
            _report_table(name, tbl)

    cap_base = data_base.get("converter_caps")
    cap_scen = data_scen.get("converter_caps")
    bal_base = data_base.get("commodity_balance_annual")
    bal_scen = data_scen.get("commodity_balance_annual")
    missing = [n for n, t in {
        "converter_caps (base)"          : cap_base,
        "converter_caps (scenario)"      : cap_scen,
        "commodity_balance_annual (base)": bal_base,
        "commodity_balance_annual (scen)": bal_scen
    }.items() if t is None]
    if missing:
        raise KeyError("Missing required tables: " + ", ".join(missing))

    cap_df  = compare_converter_capacity(cap_base, cap_scen, YEAR, LABEL_BASE, LABEL_SCEN)
    bal_df  = compare_commodity_balance(bal_base, bal_scen, YEAR, LABEL_BASE, LABEL_SCEN)
    disp_df = compare_global_dispatch(bal_base, bal_scen, LABEL_BASE, LABEL_SCEN, YEAR)

    plot_results(cap_df, bal_df, disp_df, LABEL_BASE, LABEL_SCEN)

    if MAKE_WEEKLY_PLOTS:
        # Make weekly extremes chart for the scenario file (you can switch to base_gdx_path)
        weekly_extremes_plot(scenario_gdx_path, YEAR)

    print("\nDone.\n")

if __name__ == "__main__":
    main()

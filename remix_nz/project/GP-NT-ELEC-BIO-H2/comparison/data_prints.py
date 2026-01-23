from __future__ import annotations

from pathlib import Path
import re
import pandas as pd
import warnings


from remix.framework.tools.gdx import GDXEval  # type: ignore

# ============================================================
# Config (edit only this section)
# ============================================================

CASE_RESULT_DIRS = [
    r"C:\Local\REMix\remix_nz\project\GP-NT-ELEC-BIO-H2\nz_case_GP_2020-2025-2030-2035-2040-2045-2050\result",
    r"C:\Local\REMix\remix_nz\project\GP-NT-ELEC-BIO-H2\nz_case_NT_2020-2025-2030-2035-2040-2045-2050\result",
    r"C:\Local\REMix\remix_nz\project\GP-NT-ELEC-BIO-H2\nz_case_ELEC+_2020-2025-2030-2035-2040-2045-2050\result",
    r"C:\Local\REMix\remix_nz\project\GP-NT-ELEC-BIO-H2\nz_case_BIO+_2020-2025-2030-2035-2040-2045-2050\result",
    r"C:\Local\REMix\remix_nz\project\GP-NT-ELEC-BIO-H2\nz_case_H2+_2020-2025-2030-2035-2040-2045-2050\result",
]

SCEN_ORDER = ["GP", "NT", "ELEC", "BIO", "H2"]

YEAR_INTS = [2020, 2025, 2030, 2035, 2040, 2045, 2050]

# Exclude tiny stuff everywhere
VALUE_EPS = 0.01

# Slack fuel warning threshold: 1 GWh
SLACK_WARN_GWH = 1.0

# Output folder (single combined CSV + warnings)
OUTDIR = Path("outputs")
OUT_CSV = OUTDIR / "rq_tables_long.csv"
WARN_TXT = OUTDIR / "warnings.txt"

# ============================================================
# Utilities
# ============================================================

warnings.filterwarnings(
    "ignore",
    message=r"The GAMS version .* differs from the API version .*",
    category=UserWarning,
    module=r"gams\.transfer\..*",
)

warnings.filterwarnings(
    "ignore",
    message=r"The default of observed=False is deprecated",
    category=FutureWarning,
)

def print_header(title: str) -> None:
    print("\n" + "=" * 120)
    print(title)
    print("=" * 120)

def scenario_from_dir(result_dir: str) -> str:
    s = str(result_dir).replace("\\", "/")
    m = re.search(r"/nz_case_([A-Za-z0-9]+)(?:\+)?_", s)
    return m.group(1) if m else Path(result_dir).name


def find_gdx_in_result_dir(result_dir: str) -> Path | None:
    d = Path(result_dir)
    if not d.exists():
        return None
    # prefer newest to avoid grabbing stale -old.gdx if a newer exists
    cand = list(d.glob("*.gdx")) + list(d.glob("*.gxd"))
    if not cand:
        return None
    cand = sorted(cand, key=lambda p: p.stat().st_mtime, reverse=True)
    return cand[0]

def read_symbol(results: GDXEval, name: str) -> pd.DataFrame:
    obj = results[name]  # may raise KeyError
    if isinstance(obj, pd.Series):
        df = obj.to_frame("value")
    else:
        df = obj.copy()
    if "value" not in df.columns:
        if len(df.columns) == 1:
            df = df.rename(columns={df.columns[0]: "value"})
    return df

def coerce_year_col(dfr: pd.DataFrame, col: str, newcol: str = "year") -> pd.DataFrame:
    dfr[newcol] = pd.to_numeric(dfr[col], errors="coerce")
    dfr = dfr[dfr[newcol].notna()].copy()
    dfr[newcol] = dfr[newcol].astype(int)
    return dfr

def pretty_table(df: pd.DataFrame, decimals: int = 2) -> pd.DataFrame:
    # Only formatting for console display (CSV stays numeric)
    out = df.copy()
    with pd.option_context("display.width", 240, "display.max_rows", 40, "display.max_columns", 80):
        pass
    return out.round(decimals)

def table_to_long_cells(table_id: str, title: str, unit: str, df_wide: pd.DataFrame) -> pd.DataFrame:
    d = df_wide.copy()

    # stack into one "value" column; reset_index will create N index-level cols + M col-level cols + value
    out = d.stack(dropna=False).reset_index()

    # last column is the stacked values
    out = out.rename(columns={out.columns[-1]: "value"})

    # name the index level columns
    idx_names = list(d.index.names)
    idx_names = [n if n is not None else f"idx{i}" for i, n in enumerate(idx_names)]
    for i, name in enumerate(idx_names):
        out = out.rename(columns={out.columns[i]: name})

    # name the column level columns (e.g., 'year')
    n_idx = d.index.nlevels
    col_names = list(d.columns.names)
    col_names = [n if n is not None else f"col{i}" for i, n in enumerate(col_names)]
    for j, name in enumerate(col_names):
        out = out.rename(columns={out.columns[n_idx + j]: name})

    out.insert(0, "table_id", table_id)
    out.insert(1, "title", title)
    out.insert(2, "unit", unit)
    return out


def warn(msg: str, warnings: list[str]) -> None:
    warnings.append(msg)
    print(f"WARNING: {msg}")

# ============================================================
# Load + normalize long tables
# ============================================================

def load_case(result_dir: str) -> dict[str, pd.DataFrame]:
    scen = scenario_from_dir(result_dir)
    gdx = find_gdx_in_result_dir(result_dir)
    if gdx is None:
        raise FileNotFoundError(f"No .gdx/.gxd found in: {result_dir}")

    results = GDXEval(str(gdx))

    caps = read_symbol(results, "converter_caps").reset_index()
    cba = read_symbol(results, "commodity_balance_annual").reset_index()
    ind = read_symbol(results, "indicator_accounting").reset_index()
    ind_det = read_symbol(results, "indicator_accounting_detailed").reset_index()

    caps["scenario"] = scen
    cba["scenario"] = scen
    ind["scenario"] = scen
    ind_det["scenario"] = scen

    # years
    caps = coerce_year_col(caps, "accYears", "year")
    cba = coerce_year_col(cba, "accYears", "year")

    # indicator_accounting has accYears incl. 'horizon' -> drop non-numeric
    ind = coerce_year_col(ind, "accYears", "year")

    # detailed uses 'years'
    ind_det = coerce_year_col(ind_det, "years", "year")

    # year filter
    caps = caps[caps["year"].isin(YEAR_INTS)].copy()
    cba = cba[cba["year"].isin(YEAR_INTS)].copy()
    ind = ind[ind["year"].isin(YEAR_INTS)].copy()
    ind_det = ind_det[ind_det["year"].isin(YEAR_INTS)].copy()

    return {"caps": caps, "cba": cba, "ind": ind, "ind_det": ind_det, "gdx_path": str(gdx), "scenario": scen}

# ============================================================
# RQ1 tables
# ============================================================

def rq1_table1_installed_elec_cap_gw(caps: pd.DataFrame) -> pd.DataFrame:
    # Installed electricity generation capacity (GW)
    df = caps.copy()
    df = df[(df["accNodesModel"] == "global") &
            (df["commodity"] == "Elec") &
            (df["capType"] == "total")].copy()
    df = df[df["value"].abs() > VALUE_EPS].copy()

    piv = (df.pivot_table(index=["scenario", "techs"], columns="year", values="value",
                          aggfunc="sum", observed=False)
             .fillna(0.0))

    # add Total across techs per scenario
    total = piv.groupby(level=0).sum()
    total.index = pd.MultiIndex.from_product([total.index, ["Total"]], names=["scenario", "techs"])
    out = pd.concat([piv, total]).sort_index()
    return out

def rq1_table2_elec_generation_twh(cba: pd.DataFrame) -> pd.DataFrame:
    # Electricity generation (TWh) from commodity_balance_annual net Elec, only Value>0
    df = cba.copy()
    df = df[(df["accNodesModel"] == "global") &
            (df["commodity"] == "Elec") &
            (df["balanceType"] == "net") &
            (df["value"] > 0)].copy()

    # Convert GWh -> TWh
    df["value_twh"] = df["value"] / 1000.0
    df = df[df["value_twh"] > VALUE_EPS].copy()

    piv = (df.pivot_table(index=["scenario", "techs"], columns="year", values="value_twh",
                          aggfunc="sum", observed=False)
             .fillna(0.0))

    total = piv.groupby(level=0).sum()
    total.index = pd.MultiIndex.from_product([total.index, ["Total"]], names=["scenario", "techs"])
    out = pd.concat([piv, total]).sort_index()
    return out

def rq1_table3_fuelconv_caps_gw_output(caps: pd.DataFrame) -> pd.DataFrame:
    # Total installed capacity for fuel conversion [GW_output]
    keep = {"Methanizer", "FTropschSyn", "Electrolyser", "DAC"}

    df = caps.copy()
    # make sure column exists and is string
    df["techs"] = df["techs"].astype(str)

    df = df[(df["accNodesModel"] == "global") &
            (df["capType"] == "total") &
            (df["techs"].isin(keep))].copy()

    # sanity check
    print("RQ1-3 techs present:", sorted(df["techs"].unique()))

    # NOTE: unit already "GW_output" by construction in your convention
    df = df[df["value"].abs() > VALUE_EPS].copy()

    piv = (df.pivot_table(index=["scenario", "techs", "commodity"], columns="year", values="value",
                          aggfunc="sum", observed=False)
             .fillna(0.0))

    # total across the four techs per scenario (still by commodity)
    total = piv.groupby(level=[0, 2]).sum()
    total.index = pd.MultiIndex.from_product(
        [total.index.get_level_values(0).unique(), ["Total"], total.index.get_level_values(1).unique()],
        names=["scenario", "techs", "commodity"]
    )
    # safer construction (avoid misalignment)
    total = piv.groupby(level=[0, 2]).sum()
    total.index = pd.MultiIndex.from_arrays(
        [total.index.get_level_values(0), ["Total"] * len(total), total.index.get_level_values(1)],
        names=["scenario", "techs", "commodity"]
    )

    out = pd.concat([piv, total]).sort_index()
    return out

def _fuel_supply_group(row) -> str | None:
    t = str(row["techs"])
    comm = str(row["commodity"])

    # exclude SlackFuel_* (handled elsewhere)
    if t.startswith("SlackFuel_"):
        return None

    # direct commodities
    if comm == "e-CH4":
        return "e-CH4"
    if comm == "REfuel":
        return "e-FTL fuels"
    if comm == "H2":
        return "e-hydrogen"

    # group by tech name (imports)
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

def rq1_table4_fuels_supply_twh(cba: pd.DataFrame, warnings: list[str]) -> pd.DataFrame:
    # Supply of fuels and chemicals [TWh]
    df = cba.copy()
    df = df[(df["accNodesModel"] == "global") &
            (df["balanceType"] == "net")].copy()

    # Slack warning: SlackFuel_* above 1 GWh (aggregate)
    slack = df[df["techs"].astype(str).str.startswith("SlackFuel_")].copy()
    if not slack.empty:
        slack_sum = (slack.groupby(["scenario", "year"], observed=False)["value"].sum()
                          .reset_index())
        for _, r in slack_sum.iterrows():
            if float(r["value"]) > SLACK_WARN_GWH:
                warn(f"{r['scenario']} {int(r['year'])}: SlackFuel_* total = {float(r['value']):,.2f} GWh (> {SLACK_WARN_GWH} GWh)", warnings)

    # map to groups
    df["fuel_group"] = df.apply(_fuel_supply_group, axis=1)
    df = df[df["fuel_group"].notna()].copy()

    # Convert GWh -> TWh
    df["value_twh"] = df["value"] / 1000.0
    df = df[df["value_twh"].abs() > VALUE_EPS].copy()

    piv = (df.pivot_table(index=["scenario", "fuel_group"], columns="year", values="value_twh",
                          aggfunc="sum", observed=False)
             .fillna(0.0))
    return piv.sort_index()

# ============================================================
# RQ3 tables
# ============================================================

def rq3_table1_total_annualised_cost_bil_eur(ind: pd.DataFrame) -> pd.DataFrame:
    # Total annualised cost (billion EUR) from Million EUR components FuelCost, Invest, OMFix
    keep = ["FuelCost", "Invest", "OMFix"]
    df = ind.copy()
    df = df[(df["accNodesModel"] == "global") & (df["indicator"].isin(keep))].copy()

    # Convert Million EUR -> Billion EUR
    df["value_bil_eur"] = df["value"] / 1000.0
    df = df[df["value_bil_eur"].abs() > VALUE_EPS].copy()

    piv = (df.pivot_table(index=["scenario", "indicator"], columns="year", values="value_bil_eur",
                          aggfunc="sum", observed=False)
             .fillna(0.0))

    # add Total across the three cost components
    total = piv.groupby(level=0).sum()
    total.index = pd.MultiIndex.from_product([total.index, ["TotalCost"]], names=["scenario", "indicator"])
    out = pd.concat([piv, total]).sort_index()
    return out

def rq3_table2_lcoe_eur_per_mwh(ind: pd.DataFrame, cba: pd.DataFrame) -> pd.DataFrame:
    # LCOE (EUR/MWh) = total annual system cost / total electricity produced.
    # Here: use (FuelCost+Invest+OMFix) from indicator_accounting [Million EUR]
    # divided by Elec produced from commodity_balance_annual (net Elec, only >0) [GWh].
    # => (Million EUR / GWh) == EUR/MWh. (because 1 million / 1,000 == 1,000; units cancel neatly)
    costs = ind[(ind["accNodesModel"] == "global") &
                (ind["indicator"].isin(["FuelCost", "Invest", "OMFix"]))].copy()
    costs = (costs.groupby(["scenario", "year"], observed=False)["value"].sum()
                  .reset_index()
                  .rename(columns={"value": "cost_million_eur"}))

    gen = cba[(cba["accNodesModel"] == "global") &
              (cba["commodity"] == "Elec") &
              (cba["balanceType"] == "net") &
              (cba["value"] > 0)].copy()
    gen = (gen.groupby(["scenario", "year"], observed=False)["value"].sum()
              .reset_index()
              .rename(columns={"value": "elec_gwh"}))

    df = costs.merge(gen, on=["scenario", "year"], how="outer")
    df["lcoe_eur_per_mwh"] = df["cost_million_eur"] / df["elec_gwh"]
    df = df[df["lcoe_eur_per_mwh"].notna()].copy()
    df = df[df["lcoe_eur_per_mwh"] > 0].copy()

    piv = (df.pivot_table(index="scenario", columns="year", values="lcoe_eur_per_mwh",
                          aggfunc="mean", observed=False)
             .reindex(SCEN_ORDER)
             .fillna(0.0))
    return piv

def rq3_table3_total_co2_mt(ind: pd.DataFrame) -> pd.DataFrame:
    # Total CO2 emissions (Mt CO2), input is kt CO2
    df = ind[(ind["accNodesModel"] == "global") & (ind["indicator"] == "CO2_emission")].copy()
    df["value_mt"] = df["value"] / 1000.0  # kt -> Mt
    df = df[df["value_mt"].abs() > VALUE_EPS].copy()

    # pick the CO2 indicator name that exists
    cand = [c for c in ind["indicator"].astype(str).unique() if "CO2" in c]
    print("CO2 indicators:", cand)


    piv = (df.pivot_table(index="scenario", columns="year", values="value_mt",
                          aggfunc="sum", observed=False)
             .reindex(SCEN_ORDER)
             .fillna(0.0))
    return piv

def rq3_table4_co2_by_contributor_2050_mt(ind_det: pd.DataFrame) -> pd.DataFrame:
    # Total CO2 emissions by contributor in 2050 (Mt CO2)
    # indicator_accounting_detailed: indicator, nodesModel, years, techs, value (kt CO2)
    df = ind_det[(ind_det["indicator"] == "CO2_emission") & (ind_det["year"] == 2050)].copy()

    # aggregate nodes -> "global-like"
    df = (df.groupby(["scenario", "techs"], observed=False)["value"].sum()
            .reset_index())

    # kt -> Mt
    df["value_mt"] = df["value"] / 1000.0
    df = df[df["value_mt"].abs() > VALUE_EPS].copy()

    # Category grouping (as requested; duplicates removed)
    cat_map = {
        "HeatDemand_Fossil_LF": "Heat (fossil LF)",
        "HeatDemand_Fossil_CH4": "Heat (fossil CH4)",
        "HeatDemand_REfuel": "Heat (e-fuels)",
        "HeatDemand_e-CH4": "Heat (e-CH4)",
        "TranspDemand_Fossil_LF": "Transport (fossil LF)",
        "TranspDemand_REfuel": "Transport (e-fuels)",
        "CCGT": "Power thermal",
        "GT": "Power thermal",
        "OCGT": "Power thermal",
        "Thermal_Coal": "Power thermal",
        "Thermal_Diesel": "Power thermal",
        "DAC": "DAC",
    }

    df["category"] = df["techs"].map(cat_map).fillna(df["techs"].astype(str))
    out = (df.groupby(["scenario", "category"], observed=False)["value_mt"].sum()
             .reset_index())

    piv = (out.pivot_table(index="category", columns="scenario", values="value_mt",
                           aggfunc="sum", observed=False)
              .reindex(columns=SCEN_ORDER)
              .fillna(0.0))

    # add Total row
    piv.loc["Total"] = piv.sum(axis=0)
    return piv.sort_index()

# ============================================================
# Main: build all requested tables, print + write one CSV
# ============================================================

def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    warnings: list[str] = []
    all_tables_long: list[pd.DataFrame] = []

    cases = []
    for rd in CASE_RESULT_DIRS:
        scen = scenario_from_dir(rd)
        try:
            case = load_case(rd)
            cases.append(case)
            print_header(f"Loaded {scen} | {case['gdx_path']}")
        except Exception as e:
            warn(f"{scen}: failed to load from {rd} ({e})", warnings)

    if not cases:
        raise RuntimeError("No scenarios loaded successfully.")

    # concatenate long data across scenarios
    caps = pd.concat([c["caps"] for c in cases], ignore_index=True)
    cba = pd.concat([c["cba"] for c in cases], ignore_index=True)
    ind = pd.concat([c["ind"] for c in cases], ignore_index=True)
    ind_det = pd.concat([c["ind_det"] for c in cases], ignore_index=True)

    # enforce scenario order
    for df in [caps, cba, ind, ind_det]:
        df["scenario"] = pd.Categorical(df["scenario"], categories=SCEN_ORDER, ordered=True)

    # ---------------- RQ1 ----------------
    t = rq1_table1_installed_elec_cap_gw(caps)
    print_header("RQ1-1 Installed electricity generation capacity (GW) [converter_caps, global, Elec, capType=total]")
    with pd.option_context("display.width", 240, "display.max_rows", 60, "display.max_columns", 80):
        print(t.round(2))
    all_tables_long.append(table_to_long_cells("RQ1-1", "Installed electricity generation capacity", "GW", t))

    t = rq1_table2_elec_generation_twh(cba)
    print_header("RQ1-2 Electricity generation (TWh) [commodity_balance_annual, global, Elec, balanceType=net, value>0]")
    with pd.option_context("display.width", 240, "display.max_rows", 60, "display.max_columns", 80):
        print(t.round(2))
    all_tables_long.append(table_to_long_cells("RQ1-2", "Electricity generation", "TWh", t))

    t = rq1_table3_fuelconv_caps_gw_output(caps)
    print_header("RQ1-3 Fuel conversion installed capacity (GW_output) [converter_caps, global, capType=total, tech in {Methanizer, FTropschSyn, Electrolyser, DAC}]")
    with pd.option_context("display.width", 240, "display.max_rows", 60, "display.max_columns", 80):
        print(t.round(2))
    all_tables_long.append(table_to_long_cells("RQ1-3", "Fuel conversion installed capacity", "GW_output", t))

    t = rq1_table4_fuels_supply_twh(cba, warnings)
    print_header("RQ1-4 Supply of fuels and chemicals (TWh) [commodity_balance_annual, global, balanceType=net; SlackFuel_* excluded + warnings]")
    with pd.option_context("display.width", 240, "display.max_rows", 60, "display.max_columns", 80):
        print(t.round(2))
    all_tables_long.append(table_to_long_cells("RQ1-4", "Supply of fuels and chemicals", "TWh", t))

    # ---------------- RQ3 ----------------
    t = rq3_table1_total_annualised_cost_bil_eur(ind)
    print_header("RQ3-1 Total annualised cost (billion EUR) [indicator_accounting, global; Million EUR -> billion EUR]")
    with pd.option_context("display.width", 240, "display.max_rows", 60, "display.max_columns", 80):
        print(t.round(3))
    all_tables_long.append(table_to_long_cells("RQ3-1", "Total annualised cost", "billion EUR", t))

    t = rq3_table2_lcoe_eur_per_mwh(ind, cba)
    print_header("RQ3-2 LCOE (EUR/MWh) = (FuelCost+Invest+OMFix)[Million EUR] / Elec generation[GWh]")
    with pd.option_context("display.width", 240, "display.max_rows", 60, "display.max_columns", 80):
        print(t.round(2))
    all_tables_long.append(table_to_long_cells("RQ3-2", "LCOE", "EUR/MWh", t))

    t = rq3_table3_total_co2_mt(ind)
    print_header("RQ3-3 Total CO2 emissions (Mt CO2) [indicator_accounting, global; kt -> Mt]")
    with pd.option_context("display.width", 240, "display.max_rows", 60, "display.max_columns", 80):
        print(t.round(3))
    all_tables_long.append(table_to_long_cells("RQ3-3", "Total CO2 emissions", "Mt CO2", t))

    t = rq3_table4_co2_by_contributor_2050_mt(ind_det)
    print_header("RQ3-4 CO2 emissions by contributor in 2050 (Mt CO2) [indicator_accounting_detailed; nodes aggregated]")
    with pd.option_context("display.width", 240, "display.max_rows", 80, "display.max_columns", 80):
        print(t.round(3))
    all_tables_long.append(table_to_long_cells("RQ3-4", "CO2 emissions by contributor in 2050", "Mt CO2", t))

    # ---------------- Write one combined CSV ----------------
    out_long = pd.concat(all_tables_long, ignore_index=True)
    out_long.to_csv(OUT_CSV, index=False)

    # warnings
    WARN_TXT.write_text("\n".join(warnings) + ("\n" if warnings else ""), encoding="utf-8")

    print_header(f"WROTE outputs: {OUT_CSV}")
    print(f"WROTE warnings: {WARN_TXT}")

if __name__ == "__main__":
    main()

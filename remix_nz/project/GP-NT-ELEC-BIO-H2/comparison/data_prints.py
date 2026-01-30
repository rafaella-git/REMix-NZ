import gdxpds
import pandas as pd
from pathlib import Path


CASE_RESULT_DIRS = [
    r"C:\Local\REMix\remix_nz\project\GP-NT-ELEC-BIO-H2\nz_case_GP_2020-2025-2030-2035-2040-2045-2050\result",
    r"C:\Local\REMix\remix_nz\project\GP-NT-ELEC-BIO-H2\nz_case_NT_2020-2025-2030-2035-2040-2045-2050\result",
    r"C:\Local\REMix\remix_nz\project\GP-NT-ELEC-BIO-H2\nz_case_ELEC+_2020-2025-2030-2035-2040-2045-2050\result",
    r"C:\Local\REMix\remix_nz\project\GP-NT-ELEC-BIO-H2\nz_case_BIO+_2020-2025-2030-2035-2040-2045-2050\result",
    r"C:\Local\REMix\remix_nz\project\GP-NT-ELEC-BIO-H2\nz_case_H2+_2020-2025-2030-2035-2040-2045-2050\result",
]

SCEN_ORDER = ["GP", "NT", "ELEC+", "BIO+", "H2+"]

YEAR_INTS = [2020, 2025, 2030, 2035, 2040, 2045, 2050]

DEFAULT_SYMBOLS_FAST = [
    "converter_caps",
    "commodity_balance_annual",
    "indicator_accounting",
    "indicator_accounting_detailed",
    "storage_caps",
    "transfer_caps",
]


def _banner(msg: str) -> None:
    print("=" * 120)
    print(msg)
    print("=" * 120)


def _scenario_sort_key(scen_name: str) -> int:
    for i, key in enumerate(SCEN_ORDER):
        if key in scen_name:
            return i
    return 999


def _detect_col(df: pd.DataFrame, candidates) -> str:
    cols = set(df.columns)
    for c in candidates:
        if c in cols:
            return c
    raise KeyError(f"Could not find any of {candidates} in columns {sorted(df.columns)}")


def _as_int_years(df: pd.DataFrame, ycol: str) -> pd.DataFrame:
    out = df.copy()
    out[ycol] = pd.to_numeric(out[ycol], errors="coerce").astype("Int64")
    out = out.dropna(subset=[ycol])
    out[ycol] = out[ycol].astype(int)
    return out


def _drop_all_zero_rows_cols(pv: pd.DataFrame, tol: float = 1e-12) -> pd.DataFrame:
    if pv.empty:
        return pv
    pv = pv.copy()
    pv = pv.loc[(pv.abs() > tol).any(axis=1), :]
    pv = pv.loc[:, (pv.abs() > tol).any(axis=0)]
    return pv


def _maybe_to_gw(series: pd.Series) -> pd.Series:
    vmax = float(pd.to_numeric(series, errors="coerce").max())
    if vmax > 500:
        return series / 1000.0
    return series


def _find_gdx_in_result_dir(result_dir: str | Path) -> Path:
    result_dir = Path(result_dir)
    gdxs = sorted(result_dir.glob("*.gdx"))
    if len(gdxs) == 1:
        return gdxs[0]
    if len(gdxs) > 1:
        return gdxs[-1]
    raise FileNotFoundError(f"No .gdx found in {result_dir}")


def load_scenarios(case_result_dirs, symbols=DEFAULT_SYMBOLS_FAST) -> dict:
    gdx_paths = [_find_gdx_in_result_dir(d) for d in case_result_dirs]
    scenarios = [(g.stem, g) for g in gdx_paths]
    scenarios.sort(key=lambda x: (_scenario_sort_key(x[0]), x[0]))

    data = {}
    for scen, gdx in scenarios:
        _banner(f"Loading scenario: {scen}  |  file={gdx.name}")
        raw = gdxpds.to_dataframes(str(gdx))
        sd = {}
        for sym in symbols:
            df = raw.get(sym)
            if df is None:
                df = pd.DataFrame()
            sd[sym] = df
            print(f"  {sym}: {len(df):,} rows")
        data[scen] = sd
        print()
    return data


def rq1_1_gen_capacity_gw(data: dict, year_list) -> pd.DataFrame:
    rows = []
    for scen, sd in data.items():
        df = sd.get("converter_caps", pd.DataFrame())
        if df.empty:
            continue
        ycol = _detect_col(df, ["accYears", "year", "years"])
        tcol = _detect_col(df, ["techs", "tech"])
        vcol = _detect_col(df, ["Value", "value"])
        df = _as_int_years(df, ycol)
        df = df[df[ycol].isin(year_list)]
        if df.empty:
            continue
        tmp = df[[tcol, ycol, vcol]].copy()
        tmp[vcol] = pd.to_numeric(tmp[vcol], errors="coerce").fillna(0.0)
        tmp[vcol] = _maybe_to_gw(tmp[vcol])
        tmp = tmp.groupby([tcol, ycol], as_index=False)[vcol].sum()
        tmp = tmp.rename(columns={tcol: "techs", ycol: "accYears", vcol: "Value"})
        tmp["scenario"] = scen
        rows.append(tmp)

    if not rows:
        return pd.DataFrame()

    all_df = pd.concat(rows, ignore_index=True)
    pv = all_df.pivot_table(index=["scenario", "techs"], columns="accYears", values="Value", aggfunc="sum", fill_value=0.0)
    totals = pv.groupby(level=0).sum()
    totals.index = pd.MultiIndex.from_product([totals.index, ["Total"]], names=pv.index.names)
    pv = pd.concat([pv, totals]).sort_index()
    pv = pv.reindex(sorted(pv.columns), axis=1)
    pv = _drop_all_zero_rows_cols(pv)
    return pv


def rq1_2_elec_supply_demand_twh(data: dict, year_list) -> tuple[pd.DataFrame, pd.DataFrame]:
    supply_rows = []
    demand_rows = []
    for scen, sd in data.items():
        df = sd.get("commodity_balance_annual", pd.DataFrame())
        if df.empty:
            continue
        ycol = _detect_col(df, ["accYears", "year", "years"])
        tcol = _detect_col(df, ["techs", "tech"])
        vcol = _detect_col(df, ["Value", "value"])
        ccol = None
        for cand in ["commodity", "commodities", "accCommodities", "accCommodity"]:
            if cand in df.columns:
                ccol = cand
                break

        df = _as_int_years(df, ycol)
        df = df[df[ycol].isin(year_list)]
        if ccol is not None:
            df = df[df[ccol].astype(str).str.contains("Elec", case=False, na=False)]
        if df.empty:
            continue

        tmp = df[[tcol, ycol, vcol]].copy()
        tmp[vcol] = pd.to_numeric(tmp[vcol], errors="coerce").fillna(0.0)

        pos = tmp[tmp[vcol] > 1e-12].copy()
        neg = tmp[tmp[vcol] < -1e-12].copy()

        if not pos.empty:
            pos = pos.groupby([tcol, ycol], as_index=False)[vcol].sum()
            pos = pos.rename(columns={tcol: "techs", ycol: "accYears", vcol: "Value"})
            pos["scenario"] = scen
            supply_rows.append(pos)

        if not neg.empty:
            neg = neg.groupby([tcol, ycol], as_index=False)[vcol].sum()
            neg = neg.rename(columns={tcol: "techs", ycol: "accYears", vcol: "Value"})
            neg["scenario"] = scen
            demand_rows.append(neg)

    supply = pd.DataFrame()
    demand = pd.DataFrame()
    if supply_rows:
        supply = pd.concat(supply_rows, ignore_index=True)
        supply = supply.pivot_table(index=["scenario", "techs"], columns="accYears", values="Value", aggfunc="sum", fill_value=0.0)
        supply = supply.reindex(sorted(supply.columns), axis=1)
        supply = _drop_all_zero_rows_cols(supply)
    if demand_rows:
        demand = pd.concat(demand_rows, ignore_index=True)
        demand = demand.pivot_table(index=["scenario", "techs"], columns="accYears", values="Value", aggfunc="sum", fill_value=0.0)
        demand = demand.reindex(sorted(demand.columns), axis=1)
        demand = _drop_all_zero_rows_cols(demand)
    return supply, demand


def rq1_3_fuel_conversion_capacity_gw_out(data: dict, year_list) -> pd.DataFrame:
    keep = {"DAC", "Electrolyser", "FTropschSyn", "Methanizer"}
    rows = []
    for scen, sd in data.items():
        df = sd.get("converter_caps", pd.DataFrame())
        if df.empty:
            continue
        ycol = _detect_col(df, ["accYears", "year", "years"])
        tcol = _detect_col(df, ["techs", "tech"])
        vcol = _detect_col(df, ["Value", "value"])
        df = _as_int_years(df, ycol)
        df = df[df[ycol].isin(year_list)]
        df = df[df[tcol].isin(keep)]
        if df.empty:
            continue
        tmp = df[[tcol, ycol, vcol]].copy()
        tmp[vcol] = pd.to_numeric(tmp[vcol], errors="coerce").fillna(0.0)
        tmp[vcol] = _maybe_to_gw(tmp[vcol])
        tmp = tmp.groupby([tcol, ycol], as_index=False)[vcol].sum()
        tmp = tmp.rename(columns={tcol: "techs", ycol: "accYears", vcol: "Value"})
        tmp["scenario"] = scen
        rows.append(tmp)

    if not rows:
        return pd.DataFrame()

    all_df = pd.concat(rows, ignore_index=True)
    pv = all_df.pivot_table(index=["scenario", "techs"], columns="accYears", values="Value", aggfunc="sum", fill_value=0.0)
    pv = pv.reindex(sorted(pv.columns), axis=1)
    pv = _drop_all_zero_rows_cols(pv)
    return pv


def rq3_1_system_cost_components_billion_eur(data: dict, year_list) -> tuple[pd.DataFrame, list[str]]:
    warnings = []
    rows = []
    for scen, sd in data.items():
        df = sd.get("indicator_accounting", pd.DataFrame())
        if df.empty:
            continue
        ycol = _detect_col(df, ["accYears", "year", "years"])
        icol = _detect_col(df, ["indicator", "indicators"])
        vcol = _detect_col(df, ["Value", "value"])
        ncol = None
        for cand in ["nodesModel", "accNodesModel", "nodeModel", "node", "region"]:
            if cand in df.columns:
                ncol = cand
                break

        df = _as_int_years(df, ycol)
        df = df[df[ycol].isin(year_list)]
        if ncol is not None:
            df = df[df[ncol].astype(str).str.lower().eq("global")]
        if df.empty:
            continue

        keep = ["FuelCost", "Invest", "OMFix", "SlackCost", "Slack_CO2", "SpillPenalty", "SystemCost"]
        tmp = df[df[icol].isin(keep)].copy()
        tmp[vcol] = pd.to_numeric(tmp[vcol], errors="coerce").fillna(0.0)
        tmp = tmp.groupby([icol, ycol], as_index=False)[vcol].sum()
        tmp = tmp.rename(columns={icol: "indicator", ycol: "accYears", vcol: "Value"})
        tmp["scenario"] = scen
        rows.append(tmp)

        check = tmp.pivot(index="accYears", columns="indicator", values="Value").fillna(0.0)
        for yr in year_list:
            if yr not in check.index:
                continue
            comp = float(check.loc[yr, ["FuelCost", "Invest", "OMFix", "SlackCost", "Slack_CO2", "SpillPenalty"]].sum())
            sysc = float(check.loc[yr, "SystemCost"])
            diff = sysc - comp
            if abs(diff) > 1e-3:
                warnings.append(
                    f"[SystemCost mismatch] scenario={scen} year={yr} SystemCost={sysc:,.2f} sum_components={comp:,.2f} diff={diff:,.2f} (Million EUR)"
                )

    if not rows:
        return pd.DataFrame(), warnings

    all_df = pd.concat(rows, ignore_index=True)
    pv = all_df.pivot_table(index=["scenario", "indicator"], columns="accYears", values="Value", aggfunc="sum", fill_value=0.0)
    pv = pv / 1000.0
    pv = pv.reindex(sorted(pv.columns), axis=1)
    pv = _drop_all_zero_rows_cols(pv)
    return pv, warnings


def rq3_2_generation_only_lcoe_eur_per_mwh(data: dict, year_list) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows_cost = []
    rows_gen = []

    for scen, sd in data.items():
        ia = sd.get("indicator_accounting_detailed", pd.DataFrame())
        cb = sd.get("commodity_balance_annual", pd.DataFrame())
        if ia.empty or cb.empty:
            continue

        ycol_a = _detect_col(ia, ["accYears", "year", "years"])
        icol_a = _detect_col(ia, ["indicator", "indicators"])
        tcol_a = _detect_col(ia, ["techs", "tech"])
        vcol_a = _detect_col(ia, ["Value", "value"])

        ia = _as_int_years(ia, ycol_a)
        ia = ia[ia[ycol_a].isin(year_list)]
        ia = ia[ia[icol_a].isin(["FuelCost", "Invest", "OMFix", "SlackCost"])].copy()
        if ia.empty:
            continue
        ia[vcol_a] = pd.to_numeric(ia[vcol_a], errors="coerce").fillna(0.0)

        end_use_prefixes = ("HeatDemand_", "TranspDemand_", "Elec_demand")
        ia = ia[~ia[tcol_a].astype(str).str.startswith(end_use_prefixes)]

        fuel_keep_prefixes = ("FuelImport_",)
        ia_fuel = ia[(ia[icol_a] == "FuelCost") & (ia[tcol_a].astype(str).str.startswith(fuel_keep_prefixes))]
        ia_nonfuel = ia[ia[icol_a] != "FuelCost"]
        ia = pd.concat([ia_fuel, ia_nonfuel], ignore_index=True)

        ia = ia.groupby([icol_a, ycol_a], as_index=False)[vcol_a].sum()
        ia = ia.rename(columns={icol_a: "indicator", ycol_a: "accYears", vcol_a: "Value"})
        ia["scenario"] = scen
        rows_cost.append(ia)

        ycol_c = _detect_col(cb, ["accYears", "year", "years"])
        tcol_c = _detect_col(cb, ["techs", "tech"])
        vcol_c = _detect_col(cb, ["Value", "value"])
        ccol = None
        for cand in ["commodity", "commodities", "accCommodities", "accCommodity"]:
            if cand in cb.columns:
                ccol = cand
                break

        cb = _as_int_years(cb, ycol_c)
        cb = cb[cb[ycol_c].isin(year_list)]
        if ccol is not None:
            cb = cb[cb[ccol].astype(str).str.contains("Elec", case=False, na=False)]
        cb[vcol_c] = pd.to_numeric(cb[vcol_c], errors="coerce").fillna(0.0)

        gen = cb[cb[vcol_c] > 1e-12].copy()
        gen = gen[~gen[tcol_c].astype(str).str.startswith(end_use_prefixes)]
        gen = gen.groupby([ycol_c], as_index=False)[vcol_c].sum()
        gen = gen.rename(columns={ycol_c: "accYears", vcol_c: "TWh"})
        gen["scenario"] = scen
        rows_gen.append(gen)

    if not rows_cost or not rows_gen:
        return pd.DataFrame(), pd.DataFrame()

    cost = pd.concat(rows_cost, ignore_index=True)
    gen = pd.concat(rows_gen, ignore_index=True)

    cost_pv = cost.pivot_table(index=["scenario", "indicator"], columns="accYears", values="Value", aggfunc="sum", fill_value=0.0)
    gen_pv = gen.pivot_table(index=["scenario"], columns="accYears", values="TWh", aggfunc="sum", fill_value=0.0)

    lcoe = pd.DataFrame(index=gen_pv.index, columns=gen_pv.columns, dtype=float)
    comp_rows = []
    for scen in gen_pv.index:
        for yr in gen_pv.columns:
            e_twh = float(gen_pv.loc[scen, yr])
            if abs(e_twh) < 1e-12:
                continue
            fuel = float(cost_pv.loc[(scen, "FuelCost"), yr]) if (scen, "FuelCost") in cost_pv.index else 0.0
            inv = float(cost_pv.loc[(scen, "Invest"), yr]) if (scen, "Invest") in cost_pv.index else 0.0
            om = float(cost_pv.loc[(scen, "OMFix"), yr]) if (scen, "OMFix") in cost_pv.index else 0.0
            sl = float(cost_pv.loc[(scen, "SlackCost"), yr]) if (scen, "SlackCost") in cost_pv.index else 0.0
            total_meur = fuel + inv + om + sl

            lcoe.loc[scen, yr] = total_meur / (e_twh * 1_000_000.0) * 1_000_000_000.0
            comp_rows.append(
                {
                    "scenario": scen,
                    "year": yr,
                    "Fuel": fuel / (e_twh * 1_000_000.0) * 1_000_000_000.0,
                    "Invest": inv / (e_twh * 1_000_000.0) * 1_000_000_000.0,
                    "OMFix": om / (e_twh * 1_000_000.0) * 1_000_000_000.0,
                    "Slack": sl / (e_twh * 1_000_000.0) * 1_000_000_000.0,
                }
            )

    lcoe = lcoe.reindex(sorted(lcoe.columns), axis=1)
    lcoe.index.name = "scenario"
    lcoe.columns.name = "year"
    comp = pd.DataFrame(comp_rows).set_index(["scenario", "year"]).sort_index()
    return lcoe, comp


def rq3_3_total_co2_mt(data: dict, year_list) -> pd.DataFrame:
    rows = []
    for scen, sd in data.items():
        df = sd.get("indicator_accounting", pd.DataFrame())
        if df.empty:
            continue
        ycol = _detect_col(df, ["accYears", "year", "years"])
        icol = _detect_col(df, ["indicator", "indicators"])
        vcol = _detect_col(df, ["Value", "value"])
        ncol = None
        for cand in ["nodesModel", "accNodesModel", "nodeModel", "node", "region"]:
            if cand in df.columns:
                ncol = cand
                break

        df = _as_int_years(df, ycol)
        df = df[df[ycol].isin(year_list)]
        if ncol is not None:
            df = df[df[ncol].astype(str).str.lower().eq("global")]
        df = df[df[icol].eq("CO2_emission")].copy()
        if df.empty:
            continue
        df[vcol] = pd.to_numeric(df[vcol], errors="coerce").fillna(0.0)
        tmp = df.groupby([ycol], as_index=False)[vcol].sum()
        tmp = tmp.rename(columns={ycol: "accYears", vcol: "ktCO2"})
        tmp["scenario"] = scen
        rows.append(tmp)

    if not rows:
        return pd.DataFrame()

    all_df = pd.concat(rows, ignore_index=True)
    pv = all_df.pivot_table(index="scenario", columns="accYears", values="ktCO2", aggfunc="sum", fill_value=0.0)
    pv = pv / 1000.0
    pv = pv.reindex(sorted(pv.columns), axis=1)
    pv.index.name = "scenario"
    pv.columns.name = "accYears"
    pv = _drop_all_zero_rows_cols(pv)
    return pv


def rq3_4_co2_by_contributor_2050(data: dict, year: int = 2050) -> pd.DataFrame:
    rows = []
    for scen, sd in data.items():
        df = sd.get("indicator_accounting_detailed", pd.DataFrame())
        if df.empty:
            continue
        ycol = _detect_col(df, ["accYears", "year", "years"])
        icol = _detect_col(df, ["indicator", "indicators"])
        vcol = _detect_col(df, ["Value", "value"])
        tcol = _detect_col(df, ["techs", "tech"])

        df = _as_int_years(df, ycol)
        df = df[(df[ycol] == year) & (df[icol].eq("CO2_emission"))].copy()
        if df.empty:
            continue
        df[vcol] = pd.to_numeric(df[vcol], errors="coerce").fillna(0.0)

        grp = df.groupby([tcol], as_index=False)[vcol].sum()
        grp = grp.rename(columns={tcol: "techs", vcol: "ktCO2"})
        grp["scenario"] = scen
        rows.append(grp)

    if not rows:
        return pd.DataFrame()

    all_df = pd.concat(rows, ignore_index=True)
    all_df["MtCO2"] = all_df["ktCO2"] / 1000.0
    all_df["sign"] = all_df["MtCO2"].apply(lambda x: "negative" if x < 0 else "positive")

    pv = all_df.pivot_table(index=["scenario", "sign", "techs"], values="MtCO2", aggfunc="sum", fill_value=0.0).sort_index()
    pv = _drop_all_zero_rows_cols(pv)
    return pv


def run_all(case_result_dirs):
    data = load_scenarios(case_result_dirs)

    _banner("RQ1-1 Installed electricity generation capacity (GW)")
    t11 = rq1_1_gen_capacity_gw(data, YEAR_INTS)
    print(t11 if not t11.empty else "Empty")

    _banner("RQ1-2A Electricity supply by tech (TWh)")
    t2s, t2d = rq1_2_elec_supply_demand_twh(data, YEAR_INTS)
    print(t2s if not t2s.empty else "Empty")

    _banner("RQ1-2B Electricity demand by tech (TWh)")
    print(t2d if not t2d.empty else "Empty")

    _banner("RQ1-3 Total installed capacity for fuel conversion (GW_[output])")
    t13 = rq1_3_fuel_conversion_capacity_gw_out(data, [y for y in YEAR_INTS if y != 2020])
    print(t13 if not t13.empty else "Empty")

    _banner("RQ3-1 Total annualised costs (billion EUR) [components + SystemCost separate]")
    t31, w31 = rq3_1_system_cost_components_billion_eur(data, YEAR_INTS)
    print(t31 if not t31.empty else "Empty")

    if w31:
        _banner("Warnings: SystemCost consistency")
        for w in w31:
            print(w)

    _banner("RQ3-2 Generation-only LCOE (EUR/MWh) + composition")
    lcoe, comp = rq3_2_generation_only_lcoe_eur_per_mwh(data, YEAR_INTS)
    print(lcoe if not lcoe.empty else "Empty")
    print()
    print(comp if not comp.empty else "Empty")

    _banner("RQ3-3 Total CO2 emissions (Mt CO2)")
    t33 = rq3_3_total_co2_mt(data, YEAR_INTS)
    print(t33 if not t33.empty else "Empty")

    _banner("RQ3-4 CO2 emissions by contributor in 2050 (Mt CO2) [aggregated across nodes]")
    t34 = rq3_4_co2_by_contributor_2050(data, 2050)
    print(t34 if not t34.empty else "Empty")


if __name__ == "__main__":
    run_all(CASE_RESULT_DIRS)

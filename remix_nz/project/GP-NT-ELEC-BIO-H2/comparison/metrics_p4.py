from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
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

SCEN_ORDER = ["GP", "NT", "ELEC+", "BIO+", "H2+"]
YEAR = 2050

NODE_GLOBAL = "global"
ELEC = "Elec"
H2 = "H2"

WIND_TECHS = [
    "wind_onshore_1",
    "wind_onshore_2",
    "wind_onshore_3",
    "wind_onshore_4",
    "wind_offshore_1",
    "wind_offshore_2",
    "wind_offshore_3",
    "wind_offshore_4",
]
SOLAR_TECHS = ["pv_central_fixed", "pv_decentral"]

BATTERY_TECH = "Battery"
ELEC_SLACK_TECH = "Elec_slack"
HV_TECH = "HV"
H2_STORAGE_TECH = "H2_storage"

ELECTROLYSER_TECH = "Electrolyser"
FT_TECH = "FTropschSyn"


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


def gdx_path_for_scenario(scenario: str) -> Path:
    case_dir = next((d for d in CASE_DIRS if scenario_name(d) == scenario), None)
    if case_dir is None:
        raise FileNotFoundError(f"No case_dir found for scenario: {scenario}")
    gdx_path = BASE_DIR / case_dir / "result" / f"{case_dir}.gdx"
    if not gdx_path.is_file():
        raise FileNotFoundError(f"Missing: {gdx_path}")
    return gdx_path


def _ensure_cols(df: pd.DataFrame, cols: list[str], name: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise KeyError(f"{name} missing columns {missing}. Got: {df.columns.tolist()}")


def _filter_global_year(df: pd.DataFrame, year_col: str, year: int) -> pd.DataFrame:
    out = df.copy()
    out[year_col] = out[year_col].astype(int)
    out = out[(out[year_col] == int(year)) & (out["accNodesModel"].astype(str) == NODE_GLOBAL)].copy()
    return out


def _read_cb(res: GDXEval) -> pd.DataFrame:
    cb = res["commodity_balance"].reset_index()
    _ensure_cols(cb, ["timeModel", "accNodesModel", "accYears", "techs", "commodity", "value"], "commodity_balance")
    cb = _filter_global_year(cb, "accYears", YEAR)
    cb["techs"] = cb["techs"].astype(str)
    cb["commodity"] = cb["commodity"].astype(str)
    cb["timeModel"] = cb["timeModel"].astype(str)
    cb["value"] = pd.to_numeric(cb["value"], errors="coerce").fillna(0.0).astype(float)
    return cb


def _read_caps(res: GDXEval) -> pd.DataFrame:
    caps = res["converter_caps"].reset_index()
    _ensure_cols(caps, ["accNodesModel", "accYears", "techs", "commodity", "capType", "value"], "converter_caps")
    caps = _filter_global_year(caps, "accYears", YEAR)
    caps["techs"] = caps["techs"].astype(str)
    caps["commodity"] = caps["commodity"].astype(str)
    caps["capType"] = caps["capType"].astype(str)
    caps["value"] = pd.to_numeric(caps["value"], errors="coerce").fillna(0.0).astype(float)
    return caps


def _read_tf(res: GDXEval) -> pd.DataFrame:
    tf = res["transfer_flows"].reset_index()
    _ensure_cols(tf, ["timeModel", "nodesModel_start", "nodesModel_end", "linksModel", "years", "techs", "commodity", "value"], "transfer_flows")
    tf["years"] = tf["years"].astype(int)
    tf = tf[tf["years"] == int(YEAR)].copy()
    tf["timeModel"] = tf["timeModel"].astype(str)
    tf["nodesModel_start"] = tf["nodesModel_start"].astype(str)
    tf["nodesModel_end"] = tf["nodesModel_end"].astype(str)
    tf["linksModel"] = tf["linksModel"].astype(str)
    tf["techs"] = tf["techs"].astype(str)
    tf["commodity"] = tf["commodity"].astype(str)
    tf["value"] = pd.to_numeric(tf["value"], errors="coerce").fillna(0.0).astype(float)
    return tf


def _read_ind_det(res: GDXEval) -> pd.DataFrame:
    det = res["indicator_accounting_detailed"].reset_index()
    _ensure_cols(det, ["years", "indicator", "techs", "value"], "indicator_accounting_detailed")
    det["years"] = det["years"].astype(int)
    det = det[det["years"] == int(YEAR)].copy()
    det["indicator"] = det["indicator"].astype(str)
    det["techs"] = det["techs"].astype(str)
    det["value"] = pd.to_numeric(det["value"], errors="coerce").fillna(0.0).astype(float)
    return det


def _cap_gw(caps: pd.DataFrame, tech: str, commodity: str, captype: str = "total") -> float:
    sub = caps[(caps["techs"] == tech) & (caps["commodity"] == commodity) & (caps["capType"] == captype)]
    return float(sub["value"].sum())


def _cap_group_gw(caps: pd.DataFrame, techs: list[str], commodity: str, captype: str = "total") -> float:
    sub = caps[(caps["techs"].isin(techs)) & (caps["commodity"] == commodity) & (caps["capType"] == captype)]
    return float(sub["value"].sum())


def _sum_gwh(cb: pd.DataFrame, techs: list[str] | str, commodity: str, sign: str | None = None) -> float:
    if isinstance(techs, str):
        techs = [techs]
    sub = cb[(cb["techs"].isin(techs)) & (cb["commodity"] == commodity)]
    if sign == "pos":
        sub = sub[sub["value"] > 0]
    elif sign == "neg":
        sub = sub[sub["value"] < 0]
    return float(sub["value"].sum())


def _tx_metrics(tf: pd.DataFrame, commodity: str = ELEC) -> tuple[float, float]:
    sub = tf[(tf["commodity"] == commodity)]
    abs_vals = sub["value"].abs()
    abs_sum_gwh = float(abs_vals.sum())
    peak_abs_gw = float(abs_vals.max()) if not abs_vals.empty else 0.0
    return abs_sum_gwh, peak_abs_gw


def _cycling_metrics_from_series(vals: np.ndarray) -> tuple[float, float, float]:
    # vals are hourly (GW). charge is negative in commodity_balance for Elec storage load, discharge positive.
    pos = vals[vals > 0].sum()
    neg = -vals[vals < 0].sum()
    cycles = float(min(pos, neg))
    return float(neg), float(pos), cycles


def _plant_factor(gen_gwh: float, cap_gw: float) -> float:
    if cap_gw <= 0:
        return 0.0
    return float(gen_gwh / (cap_gw * 8760.0))


def _storage_series(cb: pd.DataFrame, tech: str, commodity: str) -> np.ndarray:
    sub = cb[(cb["techs"] == tech) & (cb["commodity"] == commodity)][["timeModel", "value"]].copy()
    if sub.empty:
        return np.zeros(8760, dtype=float)
    sub["tm_i"] = sub["timeModel"].str.replace("tm", "", regex=False).astype(int) - 1
    out = np.zeros(8760, dtype=float)
    for i, v in zip(sub["tm_i"].values, sub["value"].values):
        if 0 <= i < 8760:
            out[i] += float(v)
    return out


def compute_metrics_for_scenario(scenario: str) -> dict[str, float | str]:
    gdx_path = gdx_path_for_scenario(scenario)
    res = GDXEval(str(gdx_path))

    cb = _read_cb(res)
    caps = _read_caps(res)

    try:
        tf = _read_tf(res)
    except Exception:
        tf = pd.DataFrame(columns=["commodity", "value"])

    try:
        ind_det = _read_ind_det(res)
    except Exception:
        ind_det = pd.DataFrame(columns=["indicator", "value"])

    # --- VRE ---
    wind_cap = _cap_group_gw(caps, WIND_TECHS, ELEC)
    solar_cap = _cap_group_gw(caps, SOLAR_TECHS, ELEC)

    wind_gen = _sum_gwh(cb, WIND_TECHS, ELEC, sign="pos")
    solar_gen = _sum_gwh(cb, SOLAR_TECHS, ELEC, sign="pos")

    wind_pf = _plant_factor(wind_gen, wind_cap)
    solar_pf = _plant_factor(solar_gen, solar_cap)

    # --- Unserved energy ---
    unserved = -_sum_gwh(cb, ELEC_SLACK_TECH, ELEC, sign="neg")

    # --- Battery cycling (Elec) ---
    bat_cap = _cap_gw(caps, BATTERY_TECH, ELEC)
    bat_series = _storage_series(cb, BATTERY_TECH, ELEC)
    bat_charge, bat_discharge, bat_cycles = _cycling_metrics_from_series(bat_series)
    bat_flh = (bat_discharge / bat_cap) if bat_cap > 0 else 0.0

    # --- H2 storage cycling (H2) ---
    h2_series = _storage_series(cb, H2_STORAGE_TECH, H2)
    h2_charge, h2_discharge, h2_cycles = _cycling_metrics_from_series(h2_series)

    # --- Electrolyser output (H2 commodity) ---
    ely_cap = _cap_gw(caps, ELECTROLYSER_TECH, H2)
    ely_out = _sum_gwh(cb, ELECTROLYSER_TECH, H2, sign="pos")
    ely_pf = _plant_factor(ely_out, ely_cap)

    # --- FT output (REfuel commodity) ---
    ft_cap = _cap_gw(caps, FT_TECH, "REfuel")
    ft_out = _sum_gwh(cb, FT_TECH, "REfuel", sign="pos")
    ft_pf = _plant_factor(ft_out, ft_cap)

    # --- Transmission (transfer_flows) ---
    tx_abs_sum, tx_peak_abs = _tx_metrics(tf, commodity=ELEC)

    # --- Hydro spill penalty (cost proxy) ---
    spill_penalty_meur = 0.0
    if not ind_det.empty and "indicator" in ind_det.columns:
        spill_penalty_meur = float(ind_det[ind_det["indicator"] == "SpillPenalty"]["value"].sum() * 1000.0)

    # --- HV activity from commodity_balance (system-wide proxy) ---
    hv_abs = float(cb[(cb["techs"] == HV_TECH) & (cb["commodity"] == ELEC)]["value"].abs().sum())

    return {
        "Scenario": scenario,
        "Wind cap (GW)": wind_cap,
        "Wind gen (GWh)": wind_gen,
        "Wind PF": wind_pf,
        "Solar cap (GW)": solar_cap,
        "Solar gen (GWh)": solar_gen,
        "Solar PF": solar_pf,
        "Unserved elec (GWh)": unserved,
        "Battery cap (GW)": bat_cap,
        "Battery charge (GWh)": bat_charge,
        "Battery discharge (GWh)": bat_discharge,
        "Battery cycles (GWh_throughput)": bat_cycles,
        "Battery full-load hours": bat_flh,
        "H2 storage charge (GWh)": h2_charge,
        "H2 storage discharge (GWh)": h2_discharge,
        "H2 storage cycles (GWh_throughput)": h2_cycles,
        "Electrolyser cap (GW_H2)": ely_cap,
        "Electrolyser H2 out (GWh)": ely_out,
        "Electrolyser PF(out)": ely_pf,
        "FT cap (GW_REfuel)": ft_cap,
        "FT out (GWh)": ft_out,
        "FT PF(out)": ft_pf,
        "TX abs flow sum (GWh)": tx_abs_sum,
        "TX peak abs (GW)": tx_peak_abs,
        "HV abs in system balance (GWh)": hv_abs,
        "Hydro spill penalty (MEUR, proxy)": spill_penalty_meur,
    }


def main() -> None:
    rows = []
    for scen in SCEN_ORDER:
        rows.append(compute_metrics_for_scenario(scen))

    df = pd.DataFrame(rows)

    # formatting
    float_cols = [c for c in df.columns if c != "Scenario"]
    for c in float_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", None)

    print("\n=== System-wide + tech-specific operation metrics (YEAR=2050, global where applicable) ===\n")
    print(df.to_string(index=False, justify="right", float_format=lambda x: f"{x:,.6f}"))


if __name__ == "__main__":
    main()

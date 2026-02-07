from __future__ import annotations

from pathlib import Path
from typing import Callable

import numpy as np
import matplotlib.pyplot as plt

from gams.core import gdx as gdxcc


# =========================
# Paths / scenarios
# =========================

BASE_DIR = Path(r"C:\Local\REMix\remix_nz\project\GP-NT-ELEC-BIO-H2")

CASE_DIRS = [
    "nz_case_GP_2020-2025-2030-2035-2040-2045-2050",
    "nz_case_NT_2020-2025-2030-2035-2040-2045-2050",
    "nz_case_ELEC+_2020-2025-2030-2035-2040-2045-2050",
    "nz_case_BIO+_2020-2025-2030-2035-2040-2045-2050",
    "nz_case_H2+_2020-2025-2030-2035-2040-2045-2050",
]

SCEN_ORDER = ["GP", "NT", "ELEC+", "BIO+", "H2+"]


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


# =========================
# GDX helpers (robust)
# =========================

def _as_int(x) -> int:
    if isinstance(x, int):
        return x
    if isinstance(x, (list, tuple)):
        for v in reversed(x):
            if isinstance(v, int):
                return v
    return int(x)


def tm_to_idx(tm: str) -> int:
    # tm1..tm8760 -> 0..8759
    tm = str(tm)
    if tm.startswith("tm"):
        return int(tm[2:]) - 1
    return int(tm) - 1


# =========================
# Tech grouping
# =========================

def _is_wind(t: str) -> bool:
    t = str(t)
    return "wind_" in t.lower()


def _is_solar(t: str) -> bool:
    t = str(t)
    return t in {"pv_central_fixed", "pv_decentral"} or "pv" in t.lower()


def _is_hydro(t: str) -> bool:
    t = str(t)
    return t in {"Hydro", "Hydro_reservoir"}


def _category_for_tech(tech: str) -> str:
    t = str(tech)

    if _is_wind(t):
        return "Wind"
    if _is_solar(t):
        return "Solar PV"
    if t == "Geothermal":
        return "Geothermal"
    if _is_hydro(t):
        return "Hydro"
    if t == "Battery":
        return "Battery"
    if t == "Electrolyser":
        return "Electrolysis"
    if t == "FTropschSyn":
        return "Fischer-Tropsch"
    if t == "HV":
        return "HV"
    if t == "Elec_demand":
        return "Demand"
    if t == "Elec_slack":
        return "Shedding (slack)"
    if t in {"CCGT", "GT", "OCGT", "Thermal_Bio", "Thermal_Coal", "Thermal_Diesel", "H2_CCGT", "H2_FC"}:
        return "Dispatchables"
    return "Other"


# =========================
# Streaming read (fast)
# =========================

def read_cb_categories_8760(
    gdx_path: Path,
    year: int,
    node: str = "global",
    commodity: str = "Elec",
) -> dict[str, np.ndarray]:
    out: dict[str, np.ndarray] = {}

    def get_arr(k: str) -> np.ndarray:
        if k not in out:
            out[k] = np.zeros(8760, dtype=float)
        return out[k]

    gdx_handle = gdxcc.new_gdxHandle_tp()
    ok = gdxcc.gdxCreate(gdx_handle, 0)
    if not ok:
        raise RuntimeError("Failed to create GDX handle")

    try:
        ok = gdxcc.gdxOpenRead(gdx_handle, str(gdx_path))
        if not ok:
            raise RuntimeError(f"Failed to open GDX: {gdx_path}")

        sym = _as_int(gdxcc.gdxFindSymbol(gdx_handle, "commodity_balance"))
        if sym <= 0:
            raise KeyError("Symbol 'commodity_balance' not found in GDX")

        dim = _as_int(gdxcc.gdxSymbolDim(gdx_handle, sym))
        if dim != 5:
            raise RuntimeError(f"Expected commodity_balance dim=5, got dim={dim}")

        ok = gdxcc.gdxDataReadStrStart(gdx_handle, sym)
        if not ok:
            raise RuntimeError("Failed to start reading commodity_balance")

        while True:
            rec = gdxcc.gdxDataReadStr(gdx_handle)
            rc = rec[0]
            if rc == 0:
                break
            uels = rec[1]
            vals = rec[2]

            tm, acc_node, acc_year, tech, comm = uels[0], uels[1], uels[2], uels[3], uels[4]

            if str(comm) != commodity:
                continue
            if str(acc_node) != node:
                continue
            if int(acc_year) != int(year):
                continue

            idx = tm_to_idx(tm)
            if idx < 0 or idx >= 8760:
                continue

            cat = _category_for_tech(str(tech))
            get_arr(cat)[idx] += float(vals[0])

        gdxcc.gdxDataReadDone(gdx_handle)

    finally:
        try:
            gdxcc.gdxClose(gdx_handle)
        except Exception:
            pass
        try:
            gdxcc.gdxFree(gdx_handle)
        except Exception:
            pass

    return out


# =========================
# Plot: sorted balance (like your example)
# =========================

def _pos(a: np.ndarray) -> np.ndarray:
    return np.where(a > 0, a, 0.0)


def _neg(a: np.ndarray) -> np.ndarray:
    return np.where(a < 0, a, 0.0)


def _stack_fill(ax, x, series_list, colors, labels, where_positive: bool) -> None:
    bottom = np.zeros_like(x, dtype=float)
    for s, c, lab in zip(series_list, colors, labels):
        if where_positive:
            v = _pos(s)
        else:
            v = _neg(s)
        if np.allclose(v, 0.0):
            continue
        ax.fill_between(x, bottom, bottom + v, step="pre", color=c, linewidth=0.0, label=lab)
        bottom = bottom + v


def plot_sorted_balance_like_example(
    scenario: str = "ELEC+",
    year: int = 2050,
    node: str = "global",
    commodity: str = "Elec",
) -> None:
    gdx_path = gdx_path_for_scenario(scenario)
    cats = read_cb_categories_8760(gdx_path=gdx_path, year=year, node=node, commodity=commodity)

    Z = np.zeros(8760, dtype=float)

    wind = cats.get("Wind", Z)
    solar = cats.get("Solar PV", Z)
    hydro = cats.get("Hydro", Z)
    geo = cats.get("Geothermal", Z)
    disp = cats.get("Dispatchables", Z)
    other = cats.get("Other", Z)

    battery = cats.get("Battery", Z)
    electro = cats.get("Electrolysis", Z)
    ft = cats.get("Fischer-Tropsch", Z)
    hv = cats.get("HV", Z)

    demand = cats.get("Demand", Z)
    shedding = cats.get("Shedding (slack)", Z)

    # split sign-sensitive series into nicer legend items
    batt_dis = _pos(battery)
    batt_chg = _neg(battery)

    hv_imp = _pos(hv)
    hv_exp = _neg(hv)

    # define "residual demand" (net load): demand + sinks - non-dispatchable supply
    demand_pos = -_neg(demand)  # demand is negative in commodity_balance -> positive magnitude
    sinks = (-batt_chg) + (-_neg(electro)) + (-_neg(ft)) + (-hv_exp)  # magnitudes
    nondisp_supply = _pos(wind) + _pos(solar) + _pos(hydro) + _pos(geo)

    residual = (demand_pos + sinks) - nondisp_supply  # <0 means surplus, >0 means deficit
    sort_idx = np.argsort(residual)  # most surplus (most negative) -> most deficit

    def re(a: np.ndarray) -> np.ndarray:
        return np.asarray(a, dtype=float)[sort_idx]

    x = np.arange(8760, dtype=int)

    # series for plotting (keep original signs)
    wind_s = re(wind)
    solar_s = re(solar)
    hydro_s = re(hydro)
    geo_s = re(geo)
    disp_s = re(disp)
    other_s = re(other)

    batt_dis_s = re(batt_dis)
    batt_chg_s = re(batt_chg)

    hv_imp_s = re(hv_imp)
    hv_exp_s = re(hv_exp)

    electro_s = re(electro)
    ft_s = re(ft)

    demand_s = re(demand)
    shedding_s = re(shedding)

    residual_s = re(residual)

    # plot
    fig, ax = plt.subplots(figsize=(14, 6))
    fig.subplots_adjust(left=0.08, right=0.78, top=0.94, bottom=0.12)

    # positives (generation / supply)
    pos_series = [
        wind_s,
        solar_s,
        hydro_s,
        geo_s,
        batt_dis_s,
        hv_imp_s,
        disp_s,
        other_s,
        shedding_s,  # if slack ever positive
    ]
    pos_labels = [
        "Wind",
        "Solar PV",
        "Hydro",
        "Geothermal",
        "Battery (discharge)",
        "Import (HV)",
        "Dispatchables",
        "Other",
        "Shedding (slack)",
    ]
    pos_colors = [
        "#1f77b4",  # wind
        "#ffbb00",  # solar
        "#2ca02c",  # hydro
        "#9467bd",  # geo
        "#17becf",  # batt dis
        "#ff7f0e",  # import
        "#7f7f7f",  # dispatchables
        "#c7c7c7",  # other
        "#d62728",  # slack
    ]

    # negatives (loads / sinks)
    neg_series = [
        demand_s,
        electro_s,
        ft_s,
        batt_chg_s,
        hv_exp_s,
    ]
    neg_labels = [
        "Demand",
        "Electrolysis",
        "Fischer-Tropsch",
        "Battery (charge)",
        "Export (HV)",
    ]
    neg_colors = [
        "#000000",  # demand
        "#2ca02c",  # electro
        "#8c564b",  # FT
        "#17becf",  # batt chg
        "#ff7f0e",  # export
    ]

    _stack_fill(ax, x, pos_series, pos_colors, pos_labels, where_positive=True)
    _stack_fill(ax, x, neg_series, neg_colors, neg_labels, where_positive=False)

    ax.plot(x, residual_s, color="black", linewidth=2.0, label="Residual demand")

    ax.set_xlabel("Hour of the year (sorted: most surplus → most deficit)")
    ax.set_ylabel("GW")
    ax.grid(True, axis="y", alpha=0.25)

    # legend
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False)

    plt.show()


if __name__ == "__main__":
    plot_sorted_balance_like_example(scenario="ELEC+", year=2050, node="global", commodity="Elec")

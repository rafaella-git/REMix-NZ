from __future__ import annotations

from pathlib import Path
from dataclasses import dataclass
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import FuncFormatter

from gams.core import gdx as gdxcc


BASE_DIR = Path(r"C:\Local\REMix\remix_nz\project\GP-NT-ELEC-BIO-H2")

CASE_DIRS = [
    "nz_case_GP_2020-2025-2030-2035-2040-2045-2050",
    "nz_case_NT_2020-2025-2030-2035-2040-2045-2050",
    "nz_case_ELEC+_2020-2025-2030-2035-2040-2045-2050",
    "nz_case_BIO+_2020-2025-2030-2035-2040-2045-2050",
    "nz_case_H2+_2020-2025-2030-2035-2040-2045-2050",
]

SCEN_ORDER = ["GP", "NT", "ELEC+", "BIO+", "H2+"]
ACC_NODE = "global"


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


def _as_int(x):
    if isinstance(x, int):
        return x
    if isinstance(x, (list, tuple)):
        for v in reversed(x):
            if isinstance(v, int):
                return v
    return int(x)


def _gdx_create(handle):
    try:
        return gdxcc.gdxCreate(handle, None, 0)
    except TypeError:
        try:
            return gdxcc.gdxCreate(handle, 0)
        except TypeError:
            return gdxcc.gdxCreate(handle)


def _gdx_open_read(handle, path: str):
    try:
        return gdxcc.gdxOpenRead(handle, path)
    except TypeError:
        return gdxcc.gdxOpenRead(handle, path, "")


def tm_to_idx(tm: str) -> int:
    if tm.startswith("tm"):
        return int(tm[2:]) - 1
    return int(tm) - 1


def to_carpet(arr_8760: np.ndarray) -> np.ndarray:
    return np.asarray(arr_8760, dtype=float).reshape((365, 24)).T


def _truncate_cmap(cmap, a: float, b: float, n: int = 256):
    xs = np.linspace(a, b, n)
    return LinearSegmentedColormap.from_list(f"{cmap.name}_{a:.2f}_{b:.2f}", cmap(xs))


class WhiteBandNorm(plt.Normalize):
    def __init__(self, vmin, vmax, eps=0.01, clip=False):
        super().__init__(vmin=vmin, vmax=vmax, clip=clip)
        self.eps = float(eps)

    def __call__(self, value, clip=None):
        v = np.ma.array(value, copy=False)
        vmin, vmax, eps = float(self.vmin), float(self.vmax), self.eps
        if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax == vmin:
            return np.ma.zeros_like(v, dtype=float)

        out = np.ma.empty_like(v, dtype=float)

        neg = v < -eps
        mid = (v >= -eps) & (v <= eps)
        pos = v > eps

        out[mid] = 0.5

        if np.any(neg):
            denom = (-eps - vmin) if (-eps - vmin) != 0 else 1.0
            out[neg] = 0.5 * (v[neg] - vmin) / denom

        if np.any(pos):
            denom = (vmax - eps) if (vmax - eps) != 0 else 1.0
            out[pos] = 0.5 + 0.5 * (v[pos] - eps) / denom

        return out


def diverging_no_yellow(cmap_name="RdYlBu", center_width=0.14):
    base = colormaps.get_cmap(cmap_name)
    a = 0.5 - center_width / 2
    b = 0.5 + center_width / 2
    left = base(np.linspace(0.0, a, 128))
    right = base(np.linspace(b, 1.0, 128))
    w = np.array([1.0, 1.0, 1.0, 1.0])
    left[-1] = w
    right[0] = w
    return LinearSegmentedColormap.from_list(f"{cmap_name}_noyellow_white", np.vstack([left, right]))


def pick_cmap_and_norm(values, eps_white=0.01, cmap_name="RdYlBu", center_width=0.14):
    vmin = float(np.nanmin(values))
    vmax = float(np.nanmax(values))
    cmap_full = diverging_no_yellow(cmap_name=cmap_name, center_width=center_width)

    if not np.isfinite(vmin) or not np.isfinite(vmax) or (vmin == 0 and vmax == 0):
        return cmap_full, WhiteBandNorm(-1.0, 1.0, eps=eps_white)

    if vmin >= 0:
        cmap = _truncate_cmap(cmap_full, 0.5, 1.0)
        norm = plt.Normalize(vmin=max(eps_white, 0.0), vmax=vmax if vmax > eps_white else eps_white * 2)
        return cmap, norm

    if vmax <= 0:
        cmap = _truncate_cmap(cmap_full, 0.0, 0.5)
        norm = plt.Normalize(vmin=vmin if vmin < -eps_white else -eps_white * 2, vmax=-eps_white)
        return cmap, norm

    vmax_abs = max(abs(vmin), abs(vmax))
    if not np.isfinite(vmax_abs) or vmax_abs <= 0:
        vmax_abs = 1.0
    return cmap_full, WhiteBandNorm(-vmax_abs, vmax_abs, eps=eps_white)


def _safe_set_cursor_precision(ax, decimals: int = 4):
    def _fmt(x):
        if not np.isfinite(x):
            return "nan"
        return f"{x:.{decimals}f}"
    ax.format_coord = lambda x, y: f"x={_fmt(x)}, y={_fmt(y)}"


def _format_axes_grid(fig, axes):
    for ax in axes.ravel():
        ax.set_yticks([0, 6, 12, 18, 23])
        ax.set_yticklabels(["0", "6", "12", "18", "23"])
        ax.set_xticks([0, 60, 120, 180, 240, 300, 364])
        ax.set_xticklabels(["1", "61", "121", "181", "241", "301", "365"])
    fig.text(0.04, 0.5, "Hour of day", rotation="vertical", va="center", ha="center")
    fig.text(0.5, 0.04, "Day of year", ha="center")


def _month_ticks_hours():
    month_days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    letters = list("JFMAMJJASOND")
    starts = [0]
    s = 0
    for d in month_days[:-1]:
        s += d
        starts.append(s)
    xs = [d * 24 for d in starts]
    return xs, letters


def read_symbol_8760_filter(
    gdx_path: Path,
    symbol: str,
    year: int,
    node: str,
    tech: str,
    commodity: str,
    expect_dim: int = 5,
) -> np.ndarray:
    out = np.zeros(8760, dtype=float)

    gdx_handle = gdxcc.new_gdxHandle_tp()
    rc = _gdx_create(gdx_handle)
    rc = rc[0] if isinstance(rc, (list, tuple)) else rc
    if not rc:
        raise RuntimeError("Failed to create GDX handle")

    try:
        rc = _gdx_open_read(gdx_handle, str(gdx_path))
        rc = rc[0] if isinstance(rc, (list, tuple)) else rc
        if not rc:
            raise RuntimeError(f"Failed to open GDX: {gdx_path}")

        sym_idx = _as_int(gdxcc.gdxFindSymbol(gdx_handle, symbol))
        if sym_idx <= 0:
            raise KeyError(f"Symbol '{symbol}' not found in GDX")

        dim = _as_int(gdxcc.gdxSymbolDim(gdx_handle, sym_idx))
        if dim != expect_dim:
            raise RuntimeError(f"Expected {symbol} dim={expect_dim}, got dim={dim}")

        rc = gdxcc.gdxDataReadStrStart(gdx_handle, sym_idx)
        rc = rc[0] if isinstance(rc, (list, tuple)) else rc
        if not rc:
            raise RuntimeError(f"Failed to start reading {symbol}")

        def _data_read_str(handle):
            res = gdxcc.gdxDataReadStr(handle)
            if not isinstance(res, (list, tuple)):
                return res, None, None
            rc0 = res[0]
            if isinstance(rc0, (list, tuple)):
                rc0 = rc0[0]
            uels = res[1] if len(res) > 1 else None
            vals = res[2] if len(res) > 2 else None
            return rc0, uels, vals

        while True:
            rc, uels, vals = _data_read_str(gdx_handle)
            if rc == 0:
                break

            tm, acc_node, acc_year, t, comm = uels[0], uels[1], uels[2], uels[3], uels[4]

            if str(acc_node) != node:
                continue
            if int(acc_year) != year:
                continue
            if str(t) != tech:
                continue
            if str(comm) != commodity:
                continue

            idx = tm_to_idx(str(tm))
            if not (0 <= idx < 8760):
                continue

            v = float(vals[0]) if isinstance(vals, (list, tuple)) else float(vals)
            out[idx] += v

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


def read_storage_cap_scalar(
    gdx_path: Path,
    year: int,
    node: str,
    tech: str,
    commodity: str,
    prefer: tuple[str, ...] = ("total", "upperLimit"),
) -> float | None:
    gdx_handle = gdxcc.new_gdxHandle_tp()
    rc = _gdx_create(gdx_handle)
    rc = rc[0] if isinstance(rc, (list, tuple)) else rc
    if not rc:
        raise RuntimeError("Failed to create GDX handle")

    best = None
    best_rank = None

    try:
        rc = _gdx_open_read(gdx_handle, str(gdx_path))
        rc = rc[0] if isinstance(rc, (list, tuple)) else rc
        if not rc:
            raise RuntimeError(f"Failed to open GDX: {gdx_path}")

        sym_idx = _as_int(gdxcc.gdxFindSymbol(gdx_handle, "storage_caps"))
        if sym_idx <= 0:
            return None

        dim = _as_int(gdxcc.gdxSymbolDim(gdx_handle, sym_idx))
        if dim != 5:
            return None

        rc = gdxcc.gdxDataReadStrStart(gdx_handle, sym_idx)
        rc = rc[0] if isinstance(rc, (list, tuple)) else rc
        if not rc:
            return None

        def _data_read_str(handle):
            res = gdxcc.gdxDataReadStr(handle)
            if not isinstance(res, (list, tuple)):
                return res, None, None
            rc0 = res[0]
            if isinstance(rc0, (list, tuple)):
                rc0 = rc0[0]
            uels = res[1] if len(res) > 1 else None
            vals = res[2] if len(res) > 2 else None
            return rc0, uels, vals

        while True:
            rc, uels, vals = _data_read_str(gdx_handle)
            if rc == 0:
                break

            acc_node, acc_year, t, comm, cap_type = uels[0], uels[1], uels[2], uels[3], uels[4]

            if str(acc_node) != node:
                continue
            if int(acc_year) != year:
                continue
            if str(t) != tech:
                continue
            if str(comm) != commodity:
                continue

            ct = str(cap_type)
            if ct not in prefer:
                continue

            rank = prefer.index(ct)
            v = float(vals[0]) if isinstance(vals, (list, tuple)) else float(vals)

            if best is None or rank < (best_rank if best_rank is not None else 10**9):
                best, best_rank = v, rank

        gdxcc.gdxDataReadDone(gdx_handle)
        return best

    finally:
        try:
            gdxcc.gdxClose(gdx_handle)
        except Exception:
            pass
        try:
            gdxcc.gdxFree(gdx_handle)
        except Exception:
            pass


def debug_series(name: str, x: np.ndarray, cap: float | None = None):
    x = np.asarray(x, dtype=float)
    finite = x[np.isfinite(x)]
    print(f"\n[{name}]")
    print(f"  finite_count={finite.size} / {x.size}")
    if finite.size:
        print(f"  min={finite.min():.6g} max={finite.max():.6g} mean={finite.mean():.6g} sum={finite.sum():.6g}")
        nz = np.count_nonzero(np.abs(finite) > 1e-12)
        print(f"  nonzero(|x|>1e-12)={nz}")
        print(f"  first10={np.round(finite[:10], 6)}")
    if cap is not None:
        print(f"  cap={cap} (used for % calc)")


def plot_hydro_operation_line_level(
    year: int = 2050,
    hydro_res_tech: str = "Hydro_reservoir",
    hydro_res_comm: str = "Water_in",
    hydro_turb_tech: str = "Hydro",
    hydro_turb_comm: str = "Elec",
    eps_white: float = 0.01,
):
    levels = np.zeros((len(SCEN_ORDER), 8760), dtype=float)
    gen_car = np.zeros((len(SCEN_ORDER), 24, 365), dtype=float)

    for i, scen in enumerate(SCEN_ORDER):
        gdx_path = gdx_path_for_scenario(scen)

        lvl = read_symbol_8760_filter(gdx_path, "storage_level_out", year, ACC_NODE, hydro_res_tech, hydro_res_comm)
        cap = read_storage_cap_scalar(gdx_path, year, ACC_NODE, hydro_res_tech, hydro_res_comm)
        cap = None if (cap is None or not np.isfinite(cap) or cap == 0) else float(cap)
        debug_series(f"{scen} hydro level raw ({hydro_res_tech},{hydro_res_comm})", lvl, cap)

        lvl_pct = (lvl / cap) * 100.0 if cap is not None else np.full_like(lvl, np.nan)
        lvl_pct = np.clip(lvl_pct, 0.0, 100.0)
        levels[i] = np.round(lvl_pct, 4)

        gen = read_symbol_8760_filter(gdx_path, "commodity_balance", year, ACC_NODE, hydro_turb_tech, hydro_turb_comm)
        debug_series(f"{scen} hydropower electricity ({hydro_turb_tech},{hydro_turb_comm})", gen)

        genM = to_carpet(gen)
        genM = np.round(genM, 4)
        genM[np.abs(genM) <= eps_white] = 0.0
        gen_car[i] = genM

    cmap, norm = pick_cmap_and_norm(gen_car, eps_white=eps_white)

    fig, axes = plt.subplots(len(SCEN_ORDER), 2, figsize=(14, 8.5), sharex=False, sharey=False)
    fig.subplots_adjust(left=0.08, right=0.92, top=0.92, bottom=0.12, wspace=0.20, hspace=0.40)

    xticks, xlabels = _month_ticks_hours()
    x = np.arange(8760)

    im = None
    for i, scen in enumerate(SCEN_ORDER):
        axL, axR = axes[i, 0], axes[i, 1]

        axL.plot(x, levels[i])
        axL.set_ylim(0, 100)
        axL.set_xlim(0, 8760)
        axL.set_xticks(xticks)
        axL.set_xticklabels(xlabels)
        axL.set_ylabel("(%)")
        if i == len(SCEN_ORDER) - 1:
            axL.set_xlabel("Month of the year")
        axL.set_title(f"Reservoir level ({scen})", fontsize=10, pad=6)

        im = axR.imshow(gen_car[i], origin="lower", aspect="auto", interpolation="nearest", cmap=cmap, norm=norm)
        axR.set_title(f"Hydropower output ({scen})", fontsize=10, pad=6)
        axR.set_yticks([0, 6, 12, 18, 23])
        axR.set_yticklabels(["0", "6", "12", "18", "23"])
        axR.set_xticks([0, 60, 120, 180, 240, 300, 364])
        axR.set_xticklabels(["1", "61", "121", "181", "241", "301", "365"])

        _safe_set_cursor_precision(axL, 4)
        _safe_set_cursor_precision(axR, 4)

    fig.text(0.74, 0.04, "Day of year", ha="center")
    fig.text(0.55, 0.5, "Hour of day", rotation="vertical", va="center", ha="center")

    cbar = fig.colorbar(im, ax=axes[:, 1], fraction=0.10, pad=0.02)
    cbar.set_label("Electricity (GW)")
    plt.show()


def read_transfer_caps_and_flh(gdx_path: Path, year: int, tech: str = "HV", commodity: str = "Elec"):
    def _read_table(symbol: str, expect_dim: int):
        rows = []
        gdx_handle = gdxcc.new_gdxHandle_tp()
        rc = _gdx_create(gdx_handle)
        rc = rc[0] if isinstance(rc, (list, tuple)) else rc
        if not rc:
            raise RuntimeError("Failed to create GDX handle")
        try:
            rc = _gdx_open_read(gdx_handle, str(gdx_path))
            rc = rc[0] if isinstance(rc, (list, tuple)) else rc
            if not rc:
                raise RuntimeError(f"Failed to open GDX: {gdx_path}")

            sym_idx = _as_int(gdxcc.gdxFindSymbol(gdx_handle, symbol))
            if sym_idx <= 0:
                raise KeyError(f"Symbol '{symbol}' not found in GDX")

            dim = _as_int(gdxcc.gdxSymbolDim(gdx_handle, sym_idx))
            if dim != expect_dim:
                raise RuntimeError(f"Expected {symbol} dim={expect_dim}, got dim={dim}")

            rc = gdxcc.gdxDataReadStrStart(gdx_handle, sym_idx)
            rc = rc[0] if isinstance(rc, (list, tuple)) else rc
            if not rc:
                raise RuntimeError(f"Failed to start reading {symbol}")

            def _data_read_str(handle):
                res = gdxcc.gdxDataReadStr(handle)
                if not isinstance(res, (list, tuple)):
                    return res, None, None
                rc0 = res[0]
                if isinstance(rc0, (list, tuple)):
                    rc0 = rc0[0]
                uels = res[1] if len(res) > 1 else None
                vals = res[2] if len(res) > 2 else None
                return rc0, uels, vals

            while True:
                rc0, uels, vals = _data_read_str(gdx_handle)
                if rc0 == 0:
                    break
                v = float(vals[0]) if isinstance(vals, (list, tuple)) else float(vals)
                rows.append((uels, v))

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
        return rows

    caps = {}
    for uels, v in _read_table("transfer_caps", expect_dim=7):
        a, b, link, y, t, comm, cap_type = map(str, uels)
        if int(y) != year:
            continue
        if t != tech or comm != commodity:
            continue
        if cap_type != "total":
            continue
        caps[(a, b, link)] = float(v)

    flh = {}
    for uels, v in _read_table("transfer_flows_annual", expect_dim=7):
        a, b, link, y, t, comm, flow_type = map(str, uels)
        if int(y) != year:
            continue
        if t != tech or comm != commodity:
            continue
        if flow_type != "flh":
            continue
        flh[(a, b, link)] = float(v)

    keys = sorted(set(caps.keys()) & set(flh.keys()), key=lambda k: k[2])
    labels = [k[2] for k in keys]
    cap_v = np.array([caps[k] for k in keys], dtype=float)
    flh_v = np.array([flh[k] for k in keys], dtype=float)

    return labels, cap_v, flh_v


def plot_transfers_links_scatter_two_years(
    year_left: int = 2025,
    year_right: int = 2050,
    tech: str = "HV",
    commodity: str = "Elec",
):
    years = [year_left, year_right]
    all_caps = {}
    all_flh = {}
    all_labels = set()

    for scen in SCEN_ORDER:
        gdx_path = gdx_path_for_scenario(scen)
        for y in years:
            labels, cap_v, flh_v = read_transfer_caps_and_flh(gdx_path, year=y, tech=tech, commodity=commodity)
            all_labels.update(labels)
            all_caps[(scen, y)] = dict(zip(labels, cap_v))
            all_flh[(scen, y)] = dict(zip(labels, flh_v))

    labels_sorted = sorted(all_labels)
    x = np.arange(len(labels_sorted))

    flh_all = []
    cap_all = []
    for scen in SCEN_ORDER:
        for y in years:
            flh_all.extend([all_flh[(scen, y)].get(lbl, np.nan) for lbl in labels_sorted])
            cap_all.extend([all_caps[(scen, y)].get(lbl, np.nan) for lbl in labels_sorted])

    flh_all = np.asarray(flh_all, dtype=float)
    cap_all = np.asarray(cap_all, dtype=float)

    flh_vmin = float(np.nanmin(flh_all)) if np.isfinite(np.nanmin(flh_all)) else 0.0
    flh_vmax = float(np.nanmax(flh_all)) if np.isfinite(np.nanmax(flh_all)) else 1.0
    if flh_vmin == flh_vmax:
        flh_vmin, flh_vmax = 0.0, 1.0

    cap_vmax = float(np.nanmax(cap_all)) if np.isfinite(np.nanmax(cap_all)) else 1.0
    if not np.isfinite(cap_vmax) or cap_vmax <= 0:
        cap_vmax = 1.0

    cmap = colormaps.get_cmap("PuRd")
    norm = plt.Normalize(vmin=flh_vmin, vmax=flh_vmax)

    fig, axes = plt.subplots(len(SCEN_ORDER), 2, figsize=(16, 9), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.06, right=0.92, top=0.92, bottom=0.22, wspace=0.12, hspace=0.35)

    sc_for_cbar = None
    for i, scen in enumerate(SCEN_ORDER):
        for j, y in enumerate(years):
            ax = axes[i, j]
            cap_vals = np.array([all_caps[(scen, y)].get(lbl, np.nan) for lbl in labels_sorted], dtype=float)
            flh_vals = np.array([all_flh[(scen, y)].get(lbl, np.nan) for lbl in labels_sorted], dtype=float)

            mask = np.isfinite(cap_vals) & np.isfinite(flh_vals)
            sc = ax.scatter(x[mask], cap_vals[mask], c=flh_vals[mask], cmap=cmap, norm=norm)
            sc_for_cbar = sc

            ax.set_title(f"{y} ({scen})", fontsize=10, pad=6)
            ax.set_ylim(0, cap_vmax * 1.05)
            if j == 0:
                ax.set_ylabel("Capacity (GW)")

            ax.set_xticks(x)
            ax.set_xticklabels(labels_sorted, rotation=75, ha="right")

    cbar = fig.colorbar(sc_for_cbar, ax=axes, fraction=0.10, pad=0.02)
    cbar.set_label("Full load hours (h)")

    plt.show()


@dataclass(frozen=True)
class ColSpec:
    title: str
    symbol: str
    commodity: str
    matchers: dict[str, list[str]]
    key: str
    sign_mult: float = 1.0


def read_symbol_8760_by_match(
    gdx_path: Path,
    symbol: str,
    year: int,
    node: str,
    commodity: str,
    tech_matchers: dict[str, list[str]],
    expect_dim: int = 5,
) -> dict[str, np.ndarray]:
    out = {k: np.zeros(8760, dtype=float) for k in tech_matchers}

    gdx_handle = gdxcc.new_gdxHandle_tp()
    rc = _gdx_create(gdx_handle)
    rc = rc[0] if isinstance(rc, (list, tuple)) else rc
    if not rc:
        raise RuntimeError("Failed to create GDX handle")

    try:
        rc = _gdx_open_read(gdx_handle, str(gdx_path))
        rc = rc[0] if isinstance(rc, (list, tuple)) else rc
        if not rc:
            raise RuntimeError(f"Failed to open GDX: {gdx_path}")

        sym_idx = _as_int(gdxcc.gdxFindSymbol(gdx_handle, symbol))
        if sym_idx <= 0:
            raise KeyError(f"Symbol '{symbol}' not found in GDX")

        dim = _as_int(gdxcc.gdxSymbolDim(gdx_handle, sym_idx))
        if dim != expect_dim:
            raise RuntimeError(f"Expected {symbol} dim={expect_dim}, got dim={dim}")

        rc = gdxcc.gdxDataReadStrStart(gdx_handle, sym_idx)
        rc = rc[0] if isinstance(rc, (list, tuple)) else rc
        if not rc:
            raise RuntimeError(f"Failed to start reading {symbol}")

        def _data_read_str(handle):
            res = gdxcc.gdxDataReadStr(handle)
            if not isinstance(res, (list, tuple)):
                return res, None, None
            rc0 = res[0]
            if isinstance(rc0, (list, tuple)):
                rc0 = rc0[0]
            uels = res[1] if len(res) > 1 else None
            vals = res[2] if len(res) > 2 else None
            return rc0, uels, vals

        def _match_key(tech: str) -> str | None:
            for k, pats in tech_matchers.items():
                for p in pats:
                    if tech == p or tech.startswith(p) or (p in tech):
                        return k
            return None

        while True:
            rc, uels, vals = _data_read_str(gdx_handle)
            if rc == 0:
                break

            tm, acc_node, acc_year, tech, comm = uels[0], uels[1], uels[2], uels[3], uels[4]

            if str(acc_node) != node:
                continue
            if int(acc_year) != year:
                continue
            if str(comm) != commodity:
                continue

            k = _match_key(str(tech))
            if k is None:
                continue

            idx = tm_to_idx(str(tm))
            if not (0 <= idx < 8760):
                continue

            v = float(vals[0]) if isinstance(vals, (list, tuple)) else float(vals)
            out[k][idx] += v

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


def plot_two_cols_five_rows_carpet(
    year: int,
    col_left: ColSpec,
    col_right: ColSpec,
    cbar_label: str,
    eps_white: float = 0.01,
):
    mats = np.zeros((len(SCEN_ORDER), 2, 24, 365), dtype=float)

    for i, scen in enumerate(SCEN_ORDER):
        gdx_path = gdx_path_for_scenario(scen)

        dl = read_symbol_8760_by_match(gdx_path, col_left.symbol, year, ACC_NODE, col_left.commodity, col_left.matchers)[
            col_left.key
        ]
        dr = read_symbol_8760_by_match(gdx_path, col_right.symbol, year, ACC_NODE, col_right.commodity, col_right.matchers)[
            col_right.key
        ]

        mats[i, 0] = to_carpet(col_left.sign_mult * dl)
        mats[i, 1] = to_carpet(col_right.sign_mult * dr)

    mats = np.round(mats, 4)
    mats[np.abs(mats) <= eps_white] = 0.0

    cmap, norm = pick_cmap_and_norm(mats, eps_white=eps_white)

    fig, axes = plt.subplots(len(SCEN_ORDER), 2, figsize=(14, 8.5), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.10, right=0.90, top=0.92, bottom=0.12, wspace=0.08, hspace=0.28)

    im = None
    for i, scen in enumerate(SCEN_ORDER):
        axL, axR = axes[i, 0], axes[i, 1]
        im = axL.imshow(mats[i, 0], origin="lower", aspect="auto", interpolation="nearest", cmap=cmap, norm=norm)
        axR.imshow(mats[i, 1], origin="lower", aspect="auto", interpolation="nearest", cmap=cmap, norm=norm)

        axL.set_title(f"{col_left.title} ({scen})", fontsize=10, pad=6)
        axR.set_title(f"{col_right.title} ({scen})", fontsize=10, pad=6)

        _safe_set_cursor_precision(axL, 4)
        _safe_set_cursor_precision(axR, 4)

    _format_axes_grid(fig, axes)

    cbar = fig.colorbar(im, ax=axes, fraction=0.10, pad=0.02)
    cbar.set_label(cbar_label)

    plt.show()


def plot_hydrogen_cb(year: int = 2050, pretty: bool = True) -> None:
    left_title = r"Electrolyser $H_2$ output" if pretty else "Electrolyser"
    right_title = r"Fuel cell $H_2$ consumption" if pretty else "Fuel cell"

    col_left = ColSpec(
        title=left_title,
        symbol="commodity_balance",
        commodity="H2",
        matchers={"electrolyser": ["Electrolyser", "electrolyser"]},
        key="electrolyser",
    )
    col_right = ColSpec(
        title=right_title,
        symbol="commodity_balance",
        commodity="H2",
        matchers={"h2_fc": ["H2_FC", "H2FC", "FuelCell", "Fuel cell"]},
        key="h2_fc",
    )

    plot_two_cols_five_rows_carpet(year, col_left, col_right, "Hydrogen (GW)", eps_white=0.01)


def plot_electricity_cb(year: int = 2050) -> None:
    col_left = ColSpec(
        title="DAC",
        symbol="commodity_balance",
        commodity="Elec",
        matchers={"dac": ["DAC", "Dac", "DirectAirCapture", "Direct Air Capture"]},
        key="dac",
    )
    col_right = ColSpec(
        title="Fischer–Tropsch",
        symbol="commodity_balance",
        commodity="Elec",
        matchers={"ftropsch": ["FTropschSyn", "Fischer", "FischerTropsch", "FTropsch"]},
        key="ftropsch",
    )

    plot_two_cols_five_rows_carpet(year, col_left, col_right, "Electricity (GW)", eps_white=0.01)


if __name__ == "__main__":
    plot_hydrogen_cb(year=2050, pretty=True)
    plot_electricity_cb(year=2050)

    plot_hydro_operation_line_level(
        year=2050,
        hydro_res_tech="Hydro_reservoir",
        hydro_res_comm="Water_in",
        hydro_turb_tech="Hydro",
        hydro_turb_comm="Elec",
        eps_white=0.01,
    )

    plot_transfers_links_scatter_two_years(
        year_left=2025,
        year_right=2050,
        tech="HV",
        commodity="Elec",
    )

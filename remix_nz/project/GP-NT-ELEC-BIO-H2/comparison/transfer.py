from __future__ import annotations

from pathlib import Path
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib as mpl
from gams.core import gdx as gdxcc


# ============================================================
# USER SETTINGS
# ============================================================

REGIONS_GEOJSON = r"C:\Local\REMix\remix_nz\input\shapefiles\11regionsNZ.geojson"
BASE_DIR = Path(r"C:\Local\REMix\remix_nz\project\GP-NT-ELEC-BIO-H2")

CASE_DIRS = [
    "nz_case_GP_2020-2050",
    "nz_case_NT_2020-2025-2030-2035-2040-2045-2050",
    "nz_case_ELEC+_2020-2025-2030-2035-2040-2045-2050",
    "nz_case_BIO+_2020-2025-2030-2035-2040-2045-2050",
    "nz_case_H2+_2020-2025-2030-2035-2040-2045-2050",
]

SCEN_ORDER = ["GP", "NT", "ELEC+", "BIO+", "H2+"]

YEAR_FIRST_PANEL = 2020
YEAR_OTHER_PANELS = 2050

TRANSFER_TECH = "HV"
TRANSFER_COMMODITY = "Elec"
CAP_TYPE = "total"
FLOW_NET_TYPE = "net"
FLOW_FLH_TYPE = "flh"

WIDTH_MODE = "cap"   # "net" or "cap"
SHOW_LABELS = True

LABEL_LEFT = {"AKL", "WTO", "TRN"}   # these labels appear left of centroid


# ============================================================
# GDX HELPERS
# ============================================================

def scenario_name(case_dir: str) -> str:
    if "_GP_" in case_dir: return "GP"
    if "_NT_" in case_dir: return "NT"
    if "_ELEC+_" in case_dir: return "ELEC+"
    if "_BIO+_" in case_dir: return "BIO+"
    if "_H2+_" in case_dir: return "H2+"
    return case_dir


def gdx_path_for_scenario(scenario: str) -> Path:
    case_dir = next((d for d in CASE_DIRS if scenario_name(d) == scenario), None)
    if case_dir is None:
        raise FileNotFoundError(f"No case_dir for scenario {scenario}")
    gdx_path = BASE_DIR / case_dir / "result" / f"{case_dir}.gdx"
    if not gdx_path.is_file():
        raise FileNotFoundError(f"GDX missing: {gdx_path}")
    return gdx_path


def _as_int(x):
    if isinstance(x, int): return x
    if isinstance(x, (list, tuple)):
        for v in reversed(x):
            if isinstance(v, int): return v
    return int(x)


def _gdx_create(handle):
    try: return gdxcc.gdxCreate(handle, None, 0)
    except TypeError:
        try: return gdxcc.gdxCreate(handle, 0)
        except TypeError: return gdxcc.gdxCreate(handle)


def _gdx_open_read(handle, path: str):
    try: return gdxcc.gdxOpenRead(handle, path)
    except TypeError: return gdxcc.gdxOpenRead(handle, path, "")


def _gdx_data_read_str(handle):
    res = gdxcc.gdxDataReadStr(handle)
    if not isinstance(res, (list, tuple)):
        return res, None, None
    rc0 = res[0]
    if isinstance(rc0, (list, tuple)): rc0 = rc0[0]
    uels = res[1] if len(res) > 1 else None
    vals = res[2] if len(res) > 2 else None
    return rc0, uels, vals


# ============================================================
# READ TRANSFER CAPS
# ============================================================

def read_transfer_caps(gdx_path: Path, year: int) -> dict:
    out = {}
    h = gdxcc.new_gdxHandle_tp()
    rc = _gdx_create(h)
    rc = rc[0] if isinstance(rc, (list, tuple)) else rc
    if not rc: raise RuntimeError("GDX create failed")

    try:
        rc = _gdx_open_read(h, str(gdx_path))
        rc = rc[0] if isinstance(rc, (list, tuple)) else rc
        if not rc: raise RuntimeError("GDX open failed")

        idx = _as_int(gdxcc.gdxFindSymbol(h, "transfer_caps"))
        if idx <= 0: return out

        gdxcc.gdxDataReadStrStart(h, idx)

        while True:
            rc0, uels, vals = _gdx_data_read_str(h)
            if rc0 == 0: break
            node_a, node_b, link_id, y, t, comm, ctype = map(str, uels)
            if int(y) != year: continue
            if t != TRANSFER_TECH: continue
            if comm != TRANSFER_COMMODITY: continue
            if ctype != CAP_TYPE: continue
            out[(node_a, node_b, link_id)] = float(vals[0])

        gdxcc.gdxDataReadDone(h)

    finally:
        try: gdxcc.gdxClose(h)
        except: pass
        try: gdxcc.gdxFree(h)
        except: pass

    return out


# ============================================================
# READ TRANSFER FLOWS (net + flh)
# ============================================================

def read_transfer_flows(gdx_path: Path, year: int) -> dict:
    out = {}
    h = gdxcc.new_gdxHandle_tp()
    rc = _gdx_create(h)
    rc = rc[0] if isinstance(rc, (list, tuple)) else rc
    if not rc: raise RuntimeError("GDX create failed")

    try:
        rc = _gdx_open_read(h, str(gdx_path))
        rc = rc[0] if isinstance(rc, (list, tuple)) else rc
        if not rc: raise RuntimeError("GDX open failed")

        idx = _as_int(gdxcc.gdxFindSymbol(h, "transfer_flows_annual"))
        if idx <= 0: return out

        gdxcc.gdxDataReadStrStart(h, idx)

        while True:
            rc0, uels, vals = _gdx_data_read_str(h)
            if rc0 == 0: break
            node_a, node_b, link_id, y, t, comm, ftype = map(str, uels)
            if int(y) != year: continue
            if t != TRANSFER_TECH: continue
            if comm != TRANSFER_COMMODITY: continue
            if ftype not in (FLOW_NET_TYPE, FLOW_FLH_TYPE): continue

            key = (node_a, node_b, link_id)
            if key not in out:
                out[key] = {"net": 0.0, "flh": 0.0}

            if ftype == FLOW_NET_TYPE:
                out[key]["net"] += float(vals[0])
            else:
                out[key]["flh"] += float(vals[0])

        gdxcc.gdxDataReadDone(h)

    finally:
        try: gdxcc.gdxClose(h)
        except: pass
        try: gdxcc.gdxFree(h)
        except: pass

    return out


# ============================================================
# REGIONS + ORIGINAL CENTROIDS (NO PROJECTION)
# ============================================================

def load_regions_and_centroids(path: str):
    gdf = gpd.read_file(path)

    if "node" in gdf.columns:
        gdf["region_id"] = gdf["node"].astype(str)
    else:
        gdf["region_id"] = gdf["id"].astype(str)

    # ORIGINAL centroid logic you liked
    gdf["centroid"] = gdf.geometry.centroid

    centroids = {
        rid: (pt.x, pt.y)
        for rid, pt in zip(gdf["region_id"], gdf["centroid"])
    }

    return gdf, centroids


# ============================================================
# BUILD LINK LIST
# ============================================================

def build_link_list(caps: dict, flows: dict):
    out = []
    keys = set(caps.keys()) | set(flows.keys())
    for k in keys:
        node_a, node_b, link_id = k
        cap = caps.get(k, 0.0)
        f = flows.get(k, {"net": 0.0, "flh": 0.0})
        out.append({
            "from": node_a,
            "to": node_b,
            "link": link_id,
            "cap": cap,
            "net": f["net"],
            "flh": f["flh"],
        })
    return out


# ============================================================
# PLOT ONE PANEL
# ============================================================

def plot_links(ax, regions_gdf, centroids, links, width_mode, show_labels):
    regions_gdf.plot(ax=ax, color="#f0f0f0", edgecolor="#999", linewidth=0.5)

    abs_nets = np.array([abs(l["net"]) for l in links])
    caps = np.array([l["cap"] for l in links])
    flh = np.array([l["flh"] for l in links])

    base = abs_nets if width_mode == "net" else caps
    wnorm = base / base.max() if base.max() > 0 else np.zeros_like(base)
    widths = 0.5 + 4.5 * wnorm

    flh_frac = np.clip(flh / 8760.0, 0, 1)
    cmap = plt.cm.YlOrRd
    norm = mpl.colors.Normalize(0, 1)

    for l, lw, frac in zip(links, widths, flh_frac):
        r0 = l["from"]
        r1 = l["to"]
        net = l["net"]

        if r0 not in centroids or r1 not in centroids:
            continue

        x0, y0 = centroids[r0]
        x1, y1 = centroids[r1]

        if net < 0:
            x0, y0, x1, y1 = x1, y1, x0, y0

        color = cmap(norm(frac))

        ax.annotate(
            "",
            xy=(x1, y1),
            xytext=(x0, y0),
            arrowprops=dict(
                arrowstyle="->",
                color=color,
                linewidth=lw,
                shrinkA=0,
                shrinkB=0,
                mutation_scale=10,
            ),
        )

    # Centroids
    xs = [xy[0] for xy in centroids.values()]
    ys = [xy[1] for xy in centroids.values()]
    ax.scatter(xs, ys, s=18, color="black", zorder=5)

    # Labels
    if show_labels:
        for r, (x, y) in centroids.items():
            if r in LABEL_LEFT:
                dx = -0.15
                ha = "right"
            else:
                dx = 0.15
                ha = "left"
            ax.text(x + dx, y, r, fontsize=9, ha=ha, va="center", zorder=6)

    ax.set_aspect("equal")
    ax.set_axis_off()


# ============================================================
# MULTI-PANEL PLOT
# ============================================================

def plot_multi_panel():
    regions_gdf, centroids = load_regions_and_centroids(REGIONS_GEOJSON)

    fig, axes = plt.subplots(2, 3, figsize=(12, 10))
    axes = axes.flatten()

    # Panel 0: GP 2025 → "All scenarios"
    scen0 = "GP"
    gdx0 = gdx_path_for_scenario(scen0)
    caps0 = read_transfer_caps(gdx0, YEAR_FIRST_PANEL)
    flows0 = read_transfer_flows(gdx0, YEAR_FIRST_PANEL)
    links0 = build_link_list(caps0, flows0)
    plot_links(axes[0], regions_gdf, centroids, links0, WIDTH_MODE, SHOW_LABELS)
    axes[0].set_title(f"All scenarios, {YEAR_FIRST_PANEL}", fontsize=11)

    # Panels 1–5: GP 2050, NT 2050, ELEC+ 2050, BIO+ 2050, H2+ 2050
    for i, scen in enumerate(SCEN_ORDER):
        ax = axes[i+1]
        gdx = gdx_path_for_scenario(scen)
        caps = read_transfer_caps(gdx, YEAR_OTHER_PANELS)
        flows = read_transfer_flows(gdx, YEAR_OTHER_PANELS)
        links = build_link_list(caps, flows)
        plot_links(ax, regions_gdf, centroids, links, WIDTH_MODE, SHOW_LABELS)
        ax.set_title(f"{scen}, {YEAR_OTHER_PANELS}", fontsize=11)

    # Colorbar
    cmap = plt.cm.YlOrRd
    norm = mpl.colors.Normalize(0, 1)
    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    cbar = fig.colorbar(
        sm,
        ax=[ax for ax in axes],
        fraction=0.03,
        pad=0.02,
    )
    cbar.set_label("Full load hours (% of year)", fontsize=10)
    cbar.set_ticks([0, 0.25, 0.5, 0.75, 1])
    cbar.set_ticklabels(["0", "25", "50", "75", "100"])

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    plot_multi_panel()
